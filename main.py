from fastapi import FastAPI, Form, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import RedirectResponse, JSONResponse
from sqlalchemy.orm import Session
import os
import pandas as pd

import models, schemas, services
from database import engine, get_db, SessionLocal
from chemont_export import attach_chemont_classifications, build_chemont_bundle
from descriptor_export import build_descriptor_export_fields
from step6_export import build_step6_export_fields
from workbook_contract import apply_workbook_contract, build_workbook_contract_frame
import security_check
from structure_classification import (
    apply_classification_contract,
    apply_local_taxonomy_evidence,
    classification_export_fields,
    classification_secondary_functional_groups,
    classify_structure,
)
from structural_evidence import calculate_structural_evidence
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import uuid
import time
import re
import requests
import subprocess
import tempfile
from threading import RLock

# Local, process-scoped state allows the progressive Phase B buttons to pass
# evidence bundles between requests without serialising pandas metadata into
# the browser payload.  This is deliberately an in-memory store for the local
# app; the run_id makes the lifecycle explicit and keeps B2 deterministic.
_PHASE_B_RUNS = {}
_PHASE_B_RUNS_LOCK = RLock()
_PHASE_B_ATTR_KEYS = (
    "endpoint_evidence_ledger",
    "cluster_evidence_matrix",
    "chemont_classifications",
    "uncertainty_register",
    "workbook_contract",
    "raw_extraction_tables",
)
_PHASE_B_STATE_TTL_SECONDS = 4 * 60 * 60


def _phase_b_run_id() -> str:
    return f"pb-{uuid.uuid4().hex}"


def _phase_b_plain_records(frame: pd.DataFrame) -> list[dict]:
    """Serialize row values without copying potentially large DataFrame attrs."""
    records = []
    for row in frame.to_dict(orient="records"):
        clean = {}
        for key, value in row.items():
            if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
                clean[key] = ""
            else:
                clean[key] = value
        records.append(clean)
    return records


def _capture_phase_b_state(run_id, frame, stage, trace=None, source_filename=None):
    """Persist compounds and pandas evidence attrs for the next Phase B step."""
    run_id = run_id or _phase_b_run_id()
    state = {
        "run_id": run_id,
        "stage": stage,
        "source_filename": source_filename or "report.pdf",
        "compounds": _phase_b_plain_records(frame) if isinstance(frame, pd.DataFrame) else [],
        "trace": list(trace or []),
        "attrs": {},
        "completed_stages": [],
        "updated_at": time.time(),
    }
    if isinstance(frame, pd.DataFrame):
        for key in _PHASE_B_ATTR_KEYS:
            value = frame.attrs.get(key)
            if value is not None:
                state["attrs"][key] = value.copy() if isinstance(value, pd.DataFrame) else value
    with _PHASE_B_RUNS_LOCK:
        previous = _PHASE_B_RUNS.get(run_id, {})
        # A later substage may not produce every bundle; retain earlier ones.
        merged_attrs = dict(previous.get("attrs", {}))
        merged_attrs.update(state["attrs"])
        state["attrs"] = merged_attrs
        state["completed_stages"] = list(dict.fromkeys(previous.get("completed_stages", []) + [stage]))
        state["updated_at"] = time.time()
        _PHASE_B_RUNS[run_id] = state
    return run_id


def _restore_phase_b_state(run_id, frame):
    """Reattach stored evidence DataFrames before workbook formatting."""
    if not run_id or not isinstance(frame, pd.DataFrame):
        return frame
    with _PHASE_B_RUNS_LOCK:
        state = _PHASE_B_RUNS.get(run_id)
        attrs = dict(state.get("attrs", {})) if state else {}
    for key, value in attrs.items():
        if isinstance(value, pd.DataFrame):
            frame.attrs[key] = value.copy()
        elif key == "raw_extraction_tables" and isinstance(value, list):
            frame.attrs[key] = [table.copy() if isinstance(table, pd.DataFrame) else table for table in value]
        else:
            frame.attrs[key] = value
    return frame


def _cleanup_phase_b_runs():
    cutoff = time.time() - _PHASE_B_STATE_TTL_SECONDS
    with _PHASE_B_RUNS_LOCK:
        expired = [key for key, value in _PHASE_B_RUNS.items() if value.get("updated_at", 0) < cutoff]
        for key in expired:
            _PHASE_B_RUNS.pop(key, None)


def _validate_phase_b_export(frame: pd.DataFrame, run_id: str | None) -> list[str]:
    """Return blocking B2 validation errors; external skips are non-blocking."""
    required = [
        "Compound Name", "Standardized SMILES", "Identity Status",
        "Classification Status", "Cluster ID", "Cluster Size", "Cluster_Status",
        "Domain", "Decision", "Domain_Status", "Structural_Evidence_Status",
    ]
    errors = [f"Missing required column: {column}" for column in required if column not in frame.columns]
    if errors:
        return errors
    def _valid(value):
        if value is None or pd.isna(value):
            return False
        text = str(value).strip()
        return bool(text) and Chem.MolFromSmiles(text) is not None
    invalid = ~frame["Standardized SMILES"].map(_valid)
    unsafe = frame[invalid & ((frame["Decision"].astype(str).str.upper() == "PASS") | frame["Domain"].astype(str).str.upper().str.startswith("INSIDE"))]
    if not unsafe.empty:
        errors.append(f"Unresolved rows have assessable PASS/INSIDE outcomes: {len(unsafe)}")
    retired = [column for column in ("Cramer Class", "TTC Limit", "TTC Value", "TTC Unit") if column in frame.columns]
    if retired:
        errors.append("Retired Cramer/TTC columns present: " + ", ".join(retired))
    with _PHASE_B_RUNS_LOCK:
        state = _PHASE_B_RUNS.get(run_id or "", {})
        stages = set(state.get("completed_stages", []))
        attrs = state.get("attrs", {}) if state else {}
    if run_id and not state:
        errors.append("Phase B run state is unavailable or expired; restart the PDF workflow and rerun B1a")
    if run_id and state and "B1A" not in stages:
        errors.append("B1a prerequisite is not complete")
    # B2 is released only after the explicit B1b final-enrichment checkpoint.
    # Keep the hard safety checks above (schema, structure validity, and
    # unresolved-row outcome protection) as the export gate as well.
    if run_id and state and not ({"B1B", "B1B-finalize", "B1B-chembl"} & stages):
        errors.append("B1b prerequisite is not complete")
    if run_id and state and "B1B-chemont" in stages and not isinstance(attrs.get("chemont_classifications"), pd.DataFrame):
        errors.append("ChemOnt stage is marked complete but its evidence bundle is missing")
    if run_id and state and "B1B-finalize" in stages and not isinstance(attrs.get("cluster_evidence_matrix"), pd.DataFrame):
        errors.append("B1b finalization is marked complete but cluster evidence matrix is missing")
    return errors

# Run security checks before initializing the app
# Temporarily disabled for local testing
# security_check.run_security_checks()

models.Base.metadata.create_all(bind=engine)

app = FastAPI(title="Chemical Assessment Engine", description="Chemical Identity, Functional Group, and Analog Similarity Engine")

# Temporary stable operating mode for local validation. Set the environment
# variable to "false" when external ChemOnt/EPA enrichment is intentionally
# re-enabled.
os.environ.setdefault("DETERMINISTIC_LOCAL_MODE", "true")
os.environ.setdefault("PHASE_B_STAGE_DEADLINE_SECONDS", "300")

# CORS Setup
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

import sys
if getattr(sys, 'frozen', False):
    current_dir = sys._MEIPASS
else:
    current_dir = os.path.abspath(os.path.dirname(__file__))

static_dir = os.path.join(current_dir, "static")
if not getattr(sys, 'frozen', False):
    os.makedirs(static_dir, exist_ok=True)

# Mount static files
app.mount("/static", StaticFiles(directory=static_dir, html=True), name="static")


# Isolated diagnostic only: this endpoint is not used by any assessment,
# Phase A/B workflow, cache, or workbook export.  It lets us assess whether
# the public ClassyFire lookup is reachable from this deployment.
_CLASSYFIRE_INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_CLASSYFIRE_ENTITY_URL = "https://classyfire.wishartlab.com/entities/{inchikey}.json"
_ONTOLOGY_SMOKE_CACHE_PATH = os.path.join(current_dir, "config", "ontology_smoke_cache.json")


@app.get("/diagnostics/classyfire-smoke")
def classyfire_smoke(inchikey: str = "LSXWFXONGKSEMY-UHFFFAOYSA-N"):
    normalized = (inchikey or "").strip().upper()
    if not _CLASSYFIRE_INCHIKEY_PATTERN.fullmatch(normalized):
        raise HTTPException(status_code=400, detail="Provide a valid 27-character InChIKey.")

    start = time.monotonic()
    try:
        response = requests.get(
            _CLASSYFIRE_ENTITY_URL.format(inchikey=normalized),
            headers={"Accept": "application/json", "User-Agent": "PRISM-ClassyFire-Smoke-Test/1.0"},
            timeout=(3, 7),
        )
        elapsed_ms = round((time.monotonic() - start) * 1000)
        payload = response.json() if response.headers.get("content-type", "").lower().startswith("application/json") else None
        if response.status_code != 200:
            return JSONResponse(status_code=200, content={
                "ok": False,
                "endpoint": "ClassyFire public InChIKey lookup",
                "inchikey": normalized,
                "http_status": response.status_code,
                "elapsed_ms": elapsed_ms,
                "content_type": response.headers.get("content-type", ""),
                "message": "Remote service did not return a successful classification payload.",
            })
        if not isinstance(payload, dict):
            return {"ok": False, "endpoint": "ClassyFire public InChIKey lookup", "inchikey": normalized, "http_status": 200, "elapsed_ms": elapsed_ms, "message": "Response was not a JSON object."}

        def node_name(key: str) -> str | None:
            node = payload.get(key)
            return node.get("name") if isinstance(node, dict) else None

        return {
            "ok": True,
            "endpoint": "ClassyFire public InChIKey lookup",
            "inchikey": normalized,
            "http_status": 200,
            "elapsed_ms": elapsed_ms,
            "classification_version": payload.get("classification_version"),
            "taxonomy": {
                "kingdom": node_name("kingdom"),
                "superclass": node_name("superclass"),
                "class": node_name("class"),
                "subclass": node_name("subclass"),
                "direct_parent": node_name("direct_parent"),
            },
        }
    except (requests.RequestException, ValueError) as exc:
        return {
            "ok": False,
            "endpoint": "ClassyFire public InChIKey lookup",
            "inchikey": normalized,
            "elapsed_ms": round((time.monotonic() - start) * 1000),
            "message": f"Remote request failed: {type(exc).__name__}: {str(exc)[:240]}",
        }


@app.get("/diagnostics/ontology-local-cache")
def ontology_local_cache_smoke(inchikey: str = "LSXWFXONGKSEMY-UHFFFAOYSA-N"):
    """Read a fixture-only exact-key cache; never used in PRISM processing."""
    normalized = (inchikey or "").strip().upper()
    if not _CLASSYFIRE_INCHIKEY_PATTERN.fullmatch(normalized):
        raise HTTPException(status_code=400, detail="Provide a valid 27-character InChIKey.")
    try:
        with open(_ONTOLOGY_SMOKE_CACHE_PATH, "r", encoding="utf-8") as handle:
            cache = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        return {"ok": False, "source": "Local cache fixture", "message": f"Fixture unavailable: {type(exc).__name__}"}
    record = cache.get("records", {}).get(normalized)
    return {
        "ok": record is not None,
        "source": "Local cache fixture (diagnostic only)",
        "fixture_only": bool(cache.get("fixture_only", True)),
        "cache_version": cache.get("cache_version"),
        "inchikey": normalized,
        "taxonomy": record,
        "message": "Exact InChIKey cache hit; no network request was made." if record else "No exact key in this three-record diagnostic fixture; no network request was made.",
    }


@app.get("/diagnostics/local-deterministic-classification")
def local_deterministic_classification_smoke(smiles: str = "CC(C)(C)OOC(C)(C)C"):
    """Read-only diagnostic wrapper around PRISM's local rule engine."""
    candidate = (smiles or "").strip()
    if not candidate:
        raise HTTPException(status_code=400, detail="Provide a SMILES value.")
    result = classify_structure(candidate)
    return {
        "ok": result.get("corrected_chemical_class") != "Unclassified" and not result.get("manual_review_flag", False),
        "source": "PRISM local deterministic classifier (offline)",
        "network_request_made": False,
        "input_smiles": candidate,
        "classification_structure_smiles": result.get("classification_structure_smiles"),
        "classification_standardization_version": result.get("classification_standardization_version"),
        "classification_stereochemistry_status": result.get("classification_stereochemistry_status"),
        "corrected_chemical_class": result.get("corrected_chemical_class"),
        "classification_rule_id": result.get("classification_rule_id"),
        "classification_rule_version": result.get("classification_rule_version"),
        "classification_status": result.get("classification_status"),
        "classification_scope": result.get("classification_scope"),
        "classification_review_recommended": result.get("classification_review_recommended"),
        "classification_review_recommendation": result.get("classification_review_recommendation"),
        "classification_suggested_action": result.get("classification_suggested_action"),
        "detected_feature_profile": result.get("detected_feature_profile"),
        "taxonomy_path": result.get("taxonomy_path"),
        "matched_categories": result.get("matched_categories"),
        "direct_parent": result.get("direct_parent"),
        "taxonomy_dictionary_version": result.get("taxonomy_dictionary_version"),
        "taxonomy_hierarchy": result.get("taxonomy_hierarchy"),
        "ontology_compatibility": result.get("ontology_compatibility"),
        "manual_review_flag": result.get("manual_review_flag"),
        "manual_review_reason": result.get("manual_review_reason"),
    }


@app.get("/diagnostics/pubchem-smoke")
def pubchem_smoke(inchikey: str = "LSXWFXONGKSEMY-UHFFFAOYSA-N"):
    """Isolated alternate-source probe. PubChem output is never labelled ChemOnt."""
    normalized = (inchikey or "").strip().upper()
    if not _CLASSYFIRE_INCHIKEY_PATTERN.fullmatch(normalized):
        raise HTTPException(status_code=400, detail="Provide a valid 27-character InChIKey.")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{normalized}/property/Title,CanonicalSMILES/JSON"
    start = time.monotonic()
    try:
        session = requests.Session()
        session.trust_env = False  # Tests a genuine direct route, not the inherited dummy proxy.
        response = session.get(url, headers={"Accept": "application/json"}, timeout=(3, 7))
        elapsed_ms = round((time.monotonic() - start) * 1000)
        if response.status_code != 200:
            return {"ok": False, "source": "PubChem PUG-REST (alternate, not ChemOnt)", "inchikey": normalized, "http_status": response.status_code, "elapsed_ms": elapsed_ms, "message": "PubChem did not return an identity payload."}
        props = response.json().get("PropertyTable", {}).get("Properties", [])
        return {"ok": bool(props), "source": "PubChem PUG-REST (alternate, not ChemOnt)", "inchikey": normalized, "http_status": 200, "elapsed_ms": elapsed_ms, "identity": props[0] if props else None}
    except (requests.RequestException, ValueError) as exc:
        return {"ok": False, "source": "PubChem PUG-REST (alternate, not ChemOnt)", "inchikey": normalized, "elapsed_ms": round((time.monotonic() - start) * 1000), "message": f"Alternate-source request failed: {type(exc).__name__}: {str(exc)[:240]}"}


@app.get("/diagnostics/classyfire-smoke-windows")
def classyfire_smoke_windows(inchikey: str = "LSXWFXONGKSEMY-UHFFFAOYSA-N"):
    """Isolated Windows curl/Schannel reachability probe; never used by PRISM workflows."""
    normalized = (inchikey or "").strip().upper()
    if not _CLASSYFIRE_INCHIKEY_PATTERN.fullmatch(normalized):
        raise HTTPException(status_code=400, detail="Provide a valid 27-character InChIKey.")

    url = _CLASSYFIRE_ENTITY_URL.format(inchikey=normalized)
    start = time.monotonic()
    try:
        completed = subprocess.run(
            [
                "curl.exe", "--noproxy", "*", "--connect-timeout", "5", "--max-time", "10",
                "--silent", "--show-error", "--write-out", "\\n__HTTP_STATUS__:%{http_code}",
                "--header", "Accept: application/json", url,
            ],
            capture_output=True,
            text=True,
            timeout=12,
            check=False,
        )
        elapsed_ms = round((time.monotonic() - start) * 1000)
        body, marker, status_text = completed.stdout.rpartition("\n__HTTP_STATUS__:")
        http_status = int(status_text.strip()) if marker and status_text.strip().isdigit() else None
        if completed.returncode != 0:
            return {"ok": False, "client": "Windows curl (Schannel)", "inchikey": normalized, "elapsed_ms": elapsed_ms, "curl_exit_code": completed.returncode, "message": completed.stderr.strip()[:240] or "curl returned no detail"}
        if http_status != 200:
            return {"ok": False, "client": "Windows curl (Schannel)", "inchikey": normalized, "http_status": http_status, "elapsed_ms": elapsed_ms, "message": "Remote service did not return a successful classification payload."}
        try:
            payload = json.loads(body)
        except json.JSONDecodeError:
            return {"ok": False, "client": "Windows curl (Schannel)", "inchikey": normalized, "http_status": http_status, "elapsed_ms": elapsed_ms, "message": "Response was not valid JSON."}

        def node_name(key: str) -> str | None:
            node = payload.get(key)
            return node.get("name") if isinstance(node, dict) else None

        return {
            "ok": True,
            "client": "Windows curl (Schannel)",
            "inchikey": normalized,
            "http_status": http_status,
            "elapsed_ms": elapsed_ms,
            "classification_version": payload.get("classification_version"),
            "taxonomy": {key: node_name(key) for key in ("kingdom", "superclass", "class", "subclass", "direct_parent")},
        }
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "client": "Windows curl (Schannel)", "inchikey": normalized, "elapsed_ms": round((time.monotonic() - start) * 1000), "message": f"Windows client failed: {type(exc).__name__}: {str(exc)[:240]}"}


@app.get('/assess/pdf/phase_b_state/{run_id}')
def get_phase_b_state(run_id: str):
    _cleanup_phase_b_runs()
    with _PHASE_B_RUNS_LOCK:
        state = _PHASE_B_RUNS.get(run_id)
        if not state:
            raise HTTPException(status_code=404, detail="Phase B run state not found or expired.")
        return {
            "run_id": run_id,
            "filename": state.get("source_filename", "report.pdf"),
            "compounds": state.get("compounds", []),
            "phase_b_stage": state.get("stage"),
            "completed_stages": state.get("completed_stages", []),
            "phase_b_trace": state.get("trace", []),
            "updated_at": state.get("updated_at"),
        }


def _persist_structure_classification(db: Session, identity: dict, standardized_smiles: str, original_smiles: str, classification: dict):
    """Write the classification audit row without interrupting the assessment."""
    try:
        identity_key = identity.get("inchikey") or identity.get("cas") or identity.get("name") or ""
        record = models.StructureClassificationRecord(
            identity_key=identity_key,
            input_smiles=original_smiles or "",
            standardized_smiles=standardized_smiles or "",
            classification_input_source=classification.get("classification_input_source", ""),
            corrected_chemical_class=classification.get("corrected_chemical_class", "Unclassified"),
            final_classification_record=classification.get("final_classification_record", ""),
            manual_review_flag=bool(classification.get("manual_review_flag", False)),
            manual_review_reason=classification.get("manual_review_reason", ""),
            classification_status=classification.get("classification_status", ""),
            classification_rule_id=classification.get("classification_rule_id", ""),
            classification_rule_version=classification.get("classification_rule_version", ""),
            detected_feature_profile=classification.get("detected_feature_profile", []),
        )
        audit_db = SessionLocal()
        try:
            audit_db.add(record)
            audit_db.commit()
        finally:
            audit_db.close()
    except Exception as exc:
        print(f"Structure classification audit persistence skipped: {exc}", flush=True)


def _classification_response_contract(response: schemas.AssessmentResponse) -> dict:
    """Convert an assessment response back to the shared classifier contract."""
    return {
        "corrected_chemical_class": response.corrected_chemical_class,
        "final_classification_record": response.final_classification_record,
        "manual_review_flag": response.manual_review_flag,
        "manual_review_reason": response.manual_review_reason,
        "classification_input_source": response.classification_input_source,
        "classification_structure_smiles": response.classification_structure_smiles,
        "classification_standardization_version": response.classification_standardization_version,
        "classification_stereochemistry_status": response.classification_stereochemistry_status,
        "classification_scope": response.classification_scope,
        "classification_review_recommended": response.classification_review_recommended,
        "classification_review_recommendation": response.classification_review_recommendation,
        "classification_suggested_action": response.classification_suggested_action,
        "classification_status": response.classification_status,
        "classification_rule_id": response.classification_rule_id,
        "classification_rule_version": response.classification_rule_version,
        "detected_feature_profile": response.detected_feature_profile,
        "taxonomy_path": response.taxonomy_path,
        "taxonomy_path_steps": response.taxonomy_path_steps,
        "matched_categories": response.matched_categories,
        "direct_parent": response.direct_parent,
        "taxonomy_dictionary_version": response.taxonomy_dictionary_version,
        "taxonomy_hierarchy": response.taxonomy_hierarchy,
        "topology_profile": response.topology_profile,
        "feature_profile": response.feature_profile,
        "ontology_compatibility": response.ontology_compatibility,
    }


def _resolve_chemont_bundle(inchikey: str, standardized_smiles: str) -> dict:
    """Return a normalized ChemOnt bundle for a resolved compound."""
    try:
        from chemont_adapter import enrich_with_chemont

        chemont_input = pd.DataFrame([{
            "InChIKey": inchikey or "",
            "Standardized SMILES": standardized_smiles or "",
        }])
        chemont_session = SessionLocal()
        try:
            chemont_df = enrich_with_chemont(chemont_input, chemont_session)
        finally:
            chemont_session.close()

        if isinstance(chemont_df, pd.DataFrame) and not chemont_df.empty:
            chemont_row = chemont_df.iloc[0].to_dict()
        else:
            chemont_row = {}

        return build_chemont_bundle(
            chemont_row.get("InChIKey", inchikey or ""),
            chemont_row.get("Standardized SMILES", standardized_smiles or ""),
            chemont_row.get("ChemOnt_Source", "Unavailable"),
            chemont_row.get("ChemOnt_Retrieval_Status", "Unavailable"),
            ChemOnt_Kingdom=chemont_row.get("ChemOnt_Kingdom", ""),
            ChemOnt_Superclass=chemont_row.get("ChemOnt_Superclass", ""),
            ChemOnt_Class=chemont_row.get("ChemOnt_Class", ""),
            ChemOnt_Subclass=chemont_row.get("ChemOnt_Subclass", ""),
            ChemOnt_Direct_Parent=chemont_row.get("ChemOnt_Direct_Parent", ""),
            ChemOnt_Molecular_Framework=chemont_row.get("ChemOnt_Molecular_Framework", ""),
            ChemOnt_Substituents=chemont_row.get("ChemOnt_Substituents", ""),
            ChemOnt_Classification_Version=chemont_row.get("ChemOnt_Classification_Version", ""),
        )
    except Exception as exc:
        return build_chemont_bundle(
            inchikey or "",
            standardized_smiles or "",
            "Unavailable",
            f"Unavailable: {type(exc).__name__}",
        )

@app.get("/")
def redirect_to_static():
    return RedirectResponse(url="/static/index.html")

def _direct_identity_from_request(req: schemas.AssessmentRequest) -> dict:
    """Build a workbook-style identity payload from direct row data."""
    original_smiles = (req.original_smiles or "").strip()
    standardized_smiles = (req.standardized_smiles or "").strip()
    if not original_smiles and standardized_smiles:
        original_smiles = standardized_smiles
    if not standardized_smiles and original_smiles:
        standardized_smiles = services.standardize_structure(original_smiles)

    identity_name = (req.resolved_name or req.compound_name or req.cas_number or "Direct Input").strip()
    return {
        "status": "Exact Match",
        "source": "Direct Input",
        "name": identity_name,
        "iupac_name": req.resolved_name or identity_name,
        "cid": "",
        "synonyms": [],
        "cas": req.cas_number or "",
        "smiles": original_smiles,
        "inchi": (req.inchi or "").strip(),
        "inchikey": (req.inchikey or "").strip(),
        "formula": req.molecular_formula or "",
        "molecular_weight": 0.0,
    }


def _merge_direct_identity_with_lookup(req: schemas.AssessmentRequest, direct_identity: dict, db: Session) -> dict:
    """
    Keep workbook-provided identity/structure fields, but enrich missing identity
    metadata from the same lookup path used by the Excel workflow.
    """
    lookup_identity = services.resolve_identity(req, db)
    merged = dict(direct_identity)

    if lookup_identity:
        preferred_fields = [
            "status", "source", "name", "iupac_name", "cid", "synonyms",
            "cas", "smiles", "inchi", "inchikey", "formula", "molecular_weight"
        ]
        for field in preferred_fields:
            lookup_value = lookup_identity.get(field)
            if lookup_value not in (None, "", [], {}):
                if field == "smiles":
                    # Preserve workbook-provided original/standardized structure if present.
                    if not merged.get("smiles"):
                        merged["smiles"] = lookup_value
                elif field == "source":
                    if merged.get("source") in (None, "", "Direct Input"):
                        merged["source"] = lookup_value
                else:
                    if merged.get(field) in (None, "", [], {}):
                        merged[field] = lookup_value

        # If the lookup succeeded, prefer its resolution status for audit display.
        if lookup_identity.get("status") not in (None, "", "Unresolved"):
            merged["status"] = lookup_identity["status"]

    # Use workbook name fields first, but fall back to the lookup result for display.
    merged["name"] = (
        (req.resolved_name or "").strip()
        or (req.compound_name or "").strip()
        or merged.get("name", "")
    )
    return merged


def _prism_single_assessment_contract(classification: dict) -> dict:
    """Adapt the shared PRISM result to the stable single-assessment API.

    This deliberately does not consult legacy functional-group heuristics.  It
    keeps the existing response field names so batch callers remain compatible,
    while every decision-facing value has one source of truth.
    """
    hierarchy = classification.get("taxonomy_hierarchy") or {}
    primary = str(hierarchy.get("Functional Group") or "")
    secondary_text, secondary_basis = classification_secondary_functional_groups(classification)
    feature_atoms = (classification.get("feature_profile") or {}).get("feature_atom_indices") or {}
    secondary_groups = []
    if secondary_text and secondary_text != "Not available":
        for name in (part.strip() for part in secondary_text.split(";")):
            if not name:
                continue
            secondary_groups.append({
                "functional_group": name,
                "matched_atoms": list(feature_atoms.get(name) or []),
                "SMARTS_pattern": "PRISM classifier-derived feature",
            })
    corrected_class = str(classification.get("corrected_chemical_class") or "Unclassified")
    topology = classification.get("topology_profile") or {}
    secondary_display = secondary_text if secondary_text else "Not available"
    rationale = (
        f"PRISM classification: {corrected_class}. "
        f"Primary functional group: {primary or 'Not available'}. "
        f"Secondary functional groups: {secondary_display}. "
        f"Secondary-group basis: {secondary_basis}. "
        f"Topology: {topology.get('Topology Class', 'Not available')}. "
        f"Rule: {classification.get('classification_rule_id', 'Not available')}."
    )
    return {
        # Backwards-compatible alias. The UI labels it as PRISM Classification.
        "chemical_class": corrected_class,
        "primary_functional_group": primary,
        "secondary_functional_groups": secondary_groups,
        "classification_rationale": rationale,
    }


def _scaffold_smiles(standardized_smiles: str) -> str:
    if not standardized_smiles:
        return "Not available"
    mol = Chem.MolFromSmiles(standardized_smiles)
    if mol is None:
        return "Not available"
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold)
        return scaffold_smiles if scaffold_smiles else "Not available"
    except Exception:
        return "Not available"


def _cluster_context_preview(analogs: list[dict[str, object]], scaffold_smiles: str) -> dict[str, object]:
    nearest_similarity = None
    if analogs:
        try:
            nearest_similarity = float(analogs[0].get("similarity_score", 0)) / 100.0
        except Exception:
            nearest_similarity = None

    nearest_display = round(nearest_similarity, 3) if nearest_similarity is not None else "Not available"
    return {
        "Cluster ID": "Not available in single assessment",
        "Cluster Size": "Not available in single assessment",
        "Cluster_Status": "Single-assessment preview only",
        "Domain": "Not assessable in single assessment",
        "Decision": "Not assessable in single assessment",
        "Domain_Status": "Not assessable",
        "Domain_Reasons": "Cluster and domain evidence are computed in the batch workflow; this tab shows only a nearest-analogue preview.",
        "Domain_Rule_Version": "Preview only",
        "Nearest_Neighbour_Tanimoto": nearest_display,
        "Cluster_Min_Pairwise_Tanimoto": "Not available",
        "Cluster_Median_Pairwise_Tanimoto": "Not available",
        "Cluster_Scaffold_Coverage": "Not available",
        "Cluster_Representative_Scaffold": scaffold_smiles or "Not available",
        "Property_Outlier_Flags": "Not assessable in single assessment",
        "Ionisation_Consistency": "Proxy not assessed",
        "Toxicophore_Consistency": "Proxy not assessed",
        "Uncertainty_Summary": "Batch clustering required for uncertainty register",
    }

@app.post("/assess", response_model=schemas.AssessmentResponse)
def assess_chemical(req: schemas.AssessmentRequest, db: Session = Depends(get_db)):
    # Step 1: Validate
    services.validate_input(req)
    
    # Step 2: Identity Resolution
    direct_mode = bool(req.original_smiles or req.standardized_smiles)
    if direct_mode:
        identity = _direct_identity_from_request(req)
        if req.compound_name or req.cas_number or req.resolved_name:
            identity = _merge_direct_identity_with_lookup(req, identity, db)
    else:
        identity = services.resolve_identity(req, db)
    
    # Step 3: Structure Standardization
    standardized_smiles = (req.standardized_smiles or "").strip() if direct_mode else services.standardize_structure(identity["smiles"])
    if direct_mode and not standardized_smiles and identity["smiles"]:
        standardized_smiles = services.standardize_structure(identity["smiles"])
    structure_info = {
        "original_smiles": identity["smiles"],
        "standardized_smiles": standardized_smiles,
        "image_base64": services.get_structure_image_base64(standardized_smiles)
    }
    
    # Step 4: Functional Group Detection
    functional_groups = services.detect_functional_groups(standardized_smiles)
    
    # Step 5: Parent Moiety Detection
    parent_moiety = services.detect_parent_moiety(standardized_smiles)

    # Step 5b: Structural evidence and descriptor generation
    structural_evidence = calculate_structural_evidence(standardized_smiles)
    scaffold_smiles = _scaffold_smiles(standardized_smiles)
    # Single assessment does not make external ontology calls.  The classifier
    # includes the approved local compatibility mapping in its own result.
    chemont_classification = {}

    # Deterministic local structure classification is the sole source for all
    # displayed class, hierarchy, functional-group, and review fields.
    structure_classification = classify_structure(standardized_smiles, identity["smiles"])
    structure_classification = apply_local_taxonomy_evidence(
        structure_classification, identity.get("inchikey", "")
    )
    prism_contract = _prism_single_assessment_contract(structure_classification)

    # Step 7: Molecular Fingerprint Generation
    fingerprint = services.generate_fingerprint(standardized_smiles)
    
    # Step 8, 9, 10: Analog Search & Interpretation
    analogs = services.analog_search_and_score(standardized_smiles, fingerprint, parent_moiety, functional_groups, db)
    
    # Step 11: Confidence Scoring
    confidence = services.determine_confidence(identity["status"], standardized_smiles, parent_moiety, functional_groups)
    
    # Classification rationale is derived from the same PRISM result as the
    # hierarchy. Legacy descriptors and analogue results remain supporting
    # evidence only and cannot alter this decision.
    chemical_class = prism_contract["chemical_class"]
    primary_functional_group = prism_contract["primary_functional_group"]
    secondary_functional_groups = prism_contract["secondary_functional_groups"]
    cluster_context = _cluster_context_preview(analogs, scaffold_smiles)
    classification_rationale = prism_contract["classification_rationale"]
    rationale = classification_rationale
    
    # Store Audit Record
    db_record = models.AssessmentRecord(
        input_compound_name=req.compound_name or req.resolved_name,
        input_cas_number=req.cas_number,
        input_molecular_formula=req.molecular_formula,
        input_concentration=req.concentration,
        input_concentration_unit=req.concentration_unit,
        resolved_name=identity["name"],
        resolved_cas=identity["cas"],
        resolved_smiles=identity["smiles"],
        resolved_inchi=identity["inchi"],
        resolved_inchikey=identity["inchikey"],
        resolved_formula=identity["formula"],
        resolved_mw=identity["molecular_weight"],
        resolution_status=identity["status"],
        source_database=identity["source"],
        standardized_smiles=standardized_smiles,
        functional_groups=[fg for fg in functional_groups],
        parent_moiety=parent_moiety,
        chemical_class=chemical_class,
        morgan_fingerprint="generated" if fingerprint else None,
        analogs=[a for a in analogs],
        grouping_rationale=rationale,
        confidence_score=confidence
    )
    db.add(db_record)
    _persist_structure_classification(db, identity, standardized_smiles, identity["smiles"], structure_classification)
    db.commit()
    db.refresh(db_record)
    
    return schemas.AssessmentResponse(
        identity=identity,
        structure=structure_info,
        functional_groups=functional_groups,
        parent_moiety=parent_moiety,
        chemical_class=chemical_class,
        corrected_chemical_class=structure_classification["corrected_chemical_class"],
        final_classification_record=structure_classification["final_classification_record"],
        manual_review_flag=structure_classification["manual_review_flag"],
        manual_review_reason=structure_classification["manual_review_reason"],
        classification_input_source=structure_classification["classification_input_source"],
        classification_structure_smiles=structure_classification["classification_structure_smiles"],
        classification_standardization_version=structure_classification["classification_standardization_version"],
        classification_stereochemistry_status=structure_classification["classification_stereochemistry_status"],
        classification_status=structure_classification["classification_status"],
        classification_rule_id=structure_classification["classification_rule_id"],
        classification_rule_version=structure_classification["classification_rule_version"],
        classification_scope=structure_classification["classification_scope"],
        classification_review_recommended=structure_classification["classification_review_recommended"],
        classification_review_recommendation=structure_classification["classification_review_recommendation"],
        classification_suggested_action=structure_classification["classification_suggested_action"],
        detected_feature_profile=structure_classification["detected_feature_profile"],
        primary_functional_group=primary_functional_group,
        secondary_functional_groups=secondary_functional_groups,
        taxonomy_path=structure_classification["taxonomy_path"],
        taxonomy_path_steps=structure_classification["taxonomy_path_steps"],
        matched_categories=structure_classification["matched_categories"],
        direct_parent=structure_classification["direct_parent"],
        taxonomy_dictionary_version=structure_classification["taxonomy_dictionary_version"],
        taxonomy_hierarchy=structure_classification["taxonomy_hierarchy"],
        ontology_compatibility=structure_classification["ontology_compatibility"],
        topology_profile=structure_classification["topology_profile"],
        feature_profile=structure_classification["feature_profile"],
        structural_evidence=structural_evidence,
        scaffold=scaffold_smiles,
        cramer_class="",
        ttc_limit="",
        ttc_value=0.0,
        ttc_unit="",
        cramer_rule_version="",
        cramer_path="",
        cramer_evidence="",
        toxtree_results={},
        chemont_classification=chemont_classification,
        cluster_context=cluster_context,
        final_classification_details=structure_classification["final_classification_details"],
        classification_rationale=classification_rationale,
        analogs=analogs,
        confidence=confidence,
        rationale=rationale
    )


from fastapi import UploadFile, File
import pandas as pd
from io import BytesIO
import time


@app.post('/assess/pdf/cluster_reasoning')
async def generate_cluster_reasoning(data: dict):
    """Generate a grounded explanation for one reviewer-selected compound.

    This endpoint is intentionally on-demand.  It never recalculates or
    changes clustering; it only turns the deterministic row evidence into a
    technical narrative and returns the deterministic fallback when AI access
    is unavailable.
    """
    compound = data.get("compound") if isinstance(data, dict) else None
    if not isinstance(compound, dict):
        raise HTTPException(status_code=400, detail="A single compound evidence record is required.")
    try:
        from xai_engine import generate_compound_cluster_reasoning
        reasoning = generate_compound_cluster_reasoning(compound)
        return JSONResponse(content={
            "compound_name": str(compound.get("Compound Name") or ""),
            "cluster_id": str(compound.get("Cluster ID") or ""),
            "reasoning": reasoning,
        })
    except Exception:
        raise HTTPException(status_code=500, detail="Unable to prepare the grounded cluster reasoning response.")

@app.post('/assess/excel')
async def assess_excel(file: UploadFile = File(...), db: Session = Depends(get_db)):
    content = await file.read()
    df = pd.read_excel(BytesIO(content))
    
    results = []
    chemont_columns = [
        "InChIKey",
        "Standardized SMILES",
        "ChemOnt_Kingdom",
        "ChemOnt_Superclass",
        "ChemOnt_Class",
        "ChemOnt_Subclass",
        "ChemOnt_Direct_Parent",
        "ChemOnt_Molecular_Framework",
        "ChemOnt_Substituents",
        "ChemOnt_Classification_Version",
        "ChemOnt_Source",
        "ChemOnt_Retrieval_Status",
    ]

    def _attach_chemont_fields(row_dict: dict, chemont_payload: dict | None, fallback_inchikey: str = "", fallback_smiles: str = "") -> None:
        payload = chemont_payload or {}
        row_dict["InChIKey"] = payload.get("InChIKey", fallback_inchikey)
        row_dict["Standardized SMILES"] = payload.get("Standardized SMILES", fallback_smiles)
        row_dict["ChemOnt_Kingdom"] = payload.get("ChemOnt_Kingdom", "")
        row_dict["ChemOnt_Superclass"] = payload.get("ChemOnt_Superclass", "")
        row_dict["ChemOnt_Class"] = payload.get("ChemOnt_Class", "")
        row_dict["ChemOnt_Subclass"] = payload.get("ChemOnt_Subclass", "")
        row_dict["ChemOnt_Direct_Parent"] = payload.get("ChemOnt_Direct_Parent", "")
        row_dict["ChemOnt_Molecular_Framework"] = payload.get("ChemOnt_Molecular_Framework", "")
        row_dict["ChemOnt_Substituents"] = payload.get("ChemOnt_Substituents", "")
        row_dict["ChemOnt_Classification_Version"] = payload.get("ChemOnt_Classification_Version", "")
        row_dict["ChemOnt_Source"] = payload.get("ChemOnt_Source", "")
        row_dict["ChemOnt_Retrieval_Status"] = payload.get("ChemOnt_Retrieval_Status", "")

    for _, row in df.iterrows():
        # Heuristically find columns
        compound_col = next((c for c in df.columns if 'compound' in str(c).lower() or 'chemical' in str(c).lower()), None)
        cas_col = next((c for c in df.columns if 'cas' in str(c).lower()), None)
        
        c_name = str(row[compound_col]).strip() if compound_col else ''
        if c_name == 'nan': c_name = ''
        
        c_cas = str(row[cas_col]).strip() if cas_col else ''
        import re
        if c_cas == 'nan' or not re.match(r'^\d{2,7}-\d{2}-\d$', c_cas): 
            c_cas = ''
        
        if not c_name and not c_cas:
            continue
            
        req = schemas.AssessmentRequest(
            compound_name=c_name,
            cas_number=c_cas
        )
        try:
            res = assess_chemical(req, db)
            row_dict = row.to_dict()
            row_dict['Resolved Name'] = res.identity.get('name', '')
            row_dict['Cleaned Compound Name'] = row_dict.get('Cleaned Compound Name') or row_dict.get('Compound Name', '') or c_name
            row_dict['Conc_Float'] = row_dict.get('Conc_Float', row.get('Conc_Float', ''))
            row_dict['Identity Status'] = res.identity.get('status', '')
            row_dict['Original SMILES'] = res.structure.get('original_smiles', '')
            row_dict['Standardized SMILES'] = res.structure.get('standardized_smiles', '')
            row_dict['InChIKey'] = res.identity.get('inchikey', '')
            row_dict['InChI'] = res.identity.get('inchi', '')
            row_dict['Chemical Class'] = res.chemical_class
            row_dict['Corrected Chemical Class'] = res.corrected_chemical_class
            row_dict['Final Classification Record'] = res.final_classification_record
            row_dict['Classification Input Source'] = res.classification_input_source
            row_dict['Classification Status'] = res.classification_status
            row_dict['Classification Rule ID'] = res.classification_rule_id
            row_dict['Classification Rule Version'] = res.classification_rule_version
            row_dict['Classification Scope'] = res.classification_scope
            row_dict['Classification Review Recommended'] = res.classification_review_recommended
            row_dict['Classification Review Recommendation'] = res.classification_review_recommendation
            row_dict['Manual Review Flag'] = res.manual_review_flag
            row_dict['Manual Review Reason'] = res.manual_review_reason
            row_dict['Detected Feature Profile'] = "; ".join(res.detected_feature_profile)
            row_dict['Primary Functional Group'] = res.primary_functional_group
            row_dict['Secondary Functional Groups'] = "; ".join(
                fg.functional_group for fg in res.secondary_functional_groups
            )
            row_dict.update(classification_export_fields(_classification_response_contract(res)))
            row_dict.update(build_descriptor_export_fields(res.structural_evidence, res.identity))
            row_dict.update(build_step6_export_fields(res.structural_evidence, res.scaffold))
            _attach_chemont_fields(
                row_dict,
                res.chemont_classification,
                res.identity.get('inchikey', ''),
                res.structure.get('standardized_smiles', '')
            )
            row_dict['Taxonomy Path'] = res.taxonomy_path
            row_dict['Taxonomy Path Steps'] = " → ".join(res.taxonomy_path_steps)
            row_dict['Taxonomy Hierarchy'] = " | ".join(f"{k}: {v}" for k, v in (res.taxonomy_hierarchy or {}).items())
            row_dict['Topology Profile'] = " | ".join(f"{k}: {v}" for k, v in (res.topology_profile or {}).items())
            row_dict['Classification Rationale'] = res.classification_rationale
            row_dict['Confidence'] = res.confidence
            row_dict['Tanimoto Analogs Found'] = len(res.analogs)
            results.append(row_dict)
        except Exception as e:
            row_dict = row.to_dict()
            row_dict['Resolved Name'] = 'Error'
            row_dict['Original SMILES'] = ''
            row_dict['Standardized SMILES'] = ''
            row_dict['InChIKey'] = ''
            row_dict['InChI'] = ''
            row_dict['Chemical Class'] = str(e)
            row_dict['Corrected Chemical Class'] = 'Unclassified'
            row_dict['Final Classification Record'] = ''
            row_dict['Topology Profile'] = 'Topology Class: Unclassified | Topology Modifiers: Unclassified | Ring System: Unclassified | Ring Profile: Not assessable'
            row_dict.update(build_descriptor_export_fields({
                "Structural_Evidence_Status": "Not assessable",
                "Structural_Evidence_Version": "",
                "MW": None,
                "LogP": None,
                "TPSA": None,
                "Min_EState": None,
                "Max_EState": None,
                "3D_Asphericity": None,
                "3D_PMI1": None,
                "3D_PMI2": None,
                "3D_PMI3": None,
                "3D_RadiusOfGyration": None,
                "HBD": None,
                "HBA": None,
                "Rotatable_Bonds": None,
                "Formal_Charge": None,
                "Ring_Count": None,
                "Aromatic_Ring_Count": None,
                "Fraction_CSP3": None,
                "Heavy_Atom_Count": None,
                "Molecular_Refractivity": None,
                "Ionisation_Indicator": "Unknown",
                "Toxicophore_Profile": "Not assessable",
                "All_Structural_Alerts": "Not assessable",
                "Alerts": "Not assessable",
            }, {"molecular_weight": None}))
            row_dict.update(build_step6_export_fields({
                "Structural_Evidence_Status": "Not assessable",
                "Structural_Evidence_Version": "",
                "Toxicophore_Profile": "Not assessable",
                "All_Structural_Alerts": "Not assessable",
                "Alerts": "Not assessable",
            }, "Not available"))
            _attach_chemont_fields(
                row_dict,
                {
                    "InChIKey": '',
                    "Standardized SMILES": '',
                    "ChemOnt_Kingdom": '',
                    "ChemOnt_Superclass": '',
                    "ChemOnt_Class": '',
                    "ChemOnt_Subclass": '',
                    "ChemOnt_Direct_Parent": '',
                    "ChemOnt_Molecular_Framework": '',
                    "ChemOnt_Substituents": '',
                    "ChemOnt_Classification_Version": '',
                    "ChemOnt_Source": '',
                    "ChemOnt_Retrieval_Status": 'Unavailable: processing error',
                }
            )
            results.append(row_dict)
        time.sleep(0.2) # Rate limit
        
    out_df = pd.DataFrame(results)
    if not out_df.empty:
        out_df = attach_chemont_classifications(out_df)
        out_df = apply_workbook_contract(out_df)
        out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
    out_buffer = BytesIO()
    import excel_formatter
    excel_formatter.format_excel_output(out_df, out_buffer)
    out_buffer.seek(0)
    
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        out_buffer, 
        media_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        headers={'Content-Disposition': f'attachment; filename=processed_{file.filename}'}
    )

import sys
import os

if getattr(sys, 'frozen', False):
    sys.path.insert(0, sys._MEIPASS)
else:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from extractor import UniversalPDFExtractor


def _pdf_temp_path(filename: str) -> str:
    """Return a writable, traversal-safe path for one staged PDF workflow."""
    safe_filename = os.path.basename(str(filename or "report.pdf"))
    safe_filename = re.sub(r'[<>:"/\\\\|?*]', '_', safe_filename).strip(' .') or "report.pdf"
    upload_directory = os.path.join(tempfile.gettempdir(), "PRISM", "pdf_uploads")
    os.makedirs(upload_directory, exist_ok=True)
    return os.path.join(upload_directory, f"temp_{safe_filename}")


def _load_raw_extraction_tables(filename: str) -> list[pd.DataFrame]:
    """Restore successful PDF table extractions before deduplication.

    These rows are a read-only audit view.  They retain page/table provenance
    and do not feed back into identity resolution, classification, or clustering.
    """
    raw_path = f"{_pdf_temp_path(filename)}_raw.json"
    if not os.path.exists(raw_path):
        return []
    try:
        import json
        with open(raw_path, encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, list):
            return []
        return [pd.DataFrame(table) for table in payload if isinstance(table, list) and table]
    except (OSError, ValueError, TypeError):
        return []

@app.post('/assess/pdf/step1')
async def assess_pdf_step1(file: UploadFile = File(...)):
    try:
        # Do not use the browser-provided filename directly as a filesystem
        # path. Apart from preventing path traversal, this handles filenames
        # copied from SharePoint/email that contain Windows-reserved characters
        # and used to fail here as an unhelpful HTTP 500 before extraction.
        original_filename = os.path.basename(file.filename or "report.pdf")
        safe_filename = re.sub(r'[<>:"/\\\\|?*]', '_', original_filename).strip(' .') or "report.pdf"
        if not safe_filename.lower().endswith('.pdf'):
            raise HTTPException(status_code=400, detail="Please select a PDF file.")

        content = await file.read()
        if not content:
            raise HTTPException(status_code=400, detail="The selected PDF is empty.")

        # The application folder can be read-only under beta distribution or
        # corporate endpoint controls.  Stage uploads in the user's temp area,
        # which is writable and isolated from the installed application.
        temp_pdf = _pdf_temp_path(safe_filename)
        with open(temp_pdf, 'wb') as temporary_file:
            temporary_file.write(content)

        extractor = UniversalPDFExtractor(temp_pdf)
        agency = getattr(extractor, 'agency', 'Unknown')
        analysis_type = getattr(extractor, 'analysis_type', 'Unknown')
        extractor.close()
        return {
            # The client reuses this value in steps 2–4 to locate the same
            # temporary file, so return the validated safe filename.
            "filename": safe_filename,
            "agency": agency,
            "analysis_type": analysis_type
        }
    except HTTPException:
        raise
    except Exception as e:
        if 'temp_pdf' in locals() and os.path.exists(temp_pdf):
            os.remove(temp_pdf)
        return {"error": str(e)}



from extractor import _extract_table_worker

@app.post('/assess/pdf/step234')
async def assess_pdf_step234(filename: str = Form(...), use_parallel: bool = Form(True)):
    temp_pdf = _pdf_temp_path(filename)
    if not os.path.exists(temp_pdf):
        return {"error": "Session expired or file not found."}
        
    def event_generator():
        try:
            yield json.dumps({"type": "progress", "message": "Pre-scanning PDF for tables..."}) + "\n"
            
            extractor = UniversalPDFExtractor(temp_pdf)
            pages_subset = None
            
            if extractor.agency in ['Namsa', 'Wuxi']:
                import fitz
                try:
                    doc = fitz.open(temp_pdf)
                    pages_subset = set()
                    for page_num in range(min(150, len(doc))):
                        text = doc.load_page(page_num).get_text("text").lower()
                        if 'table' in text and 'results' in text and any(m in text for m in ['gc', 'lc', 'icp', 'qtof', 'direct inject']):
                            pages_subset.add(page_num)
                    doc.close()
                except Exception as e:
                    print(f"Fitz prescan error: {e}")
                    
            import threading
            class CaptionThread(threading.Thread):
                def __init__(self, extractor, pages_subset):
                    super().__init__()
                    self.extractor = extractor
                    self.pages_subset = pages_subset
                    self.captions = None
                    self.exc = None
                def run(self):
                    try:
                        self.captions = self.extractor.find_captions(pages_subset=self.pages_subset)
                    except Exception as e:
                        self.exc = e

            cap_thread = CaptionThread(extractor, pages_subset)
            cap_thread.start()
            
            while cap_thread.is_alive():
                yield json.dumps({"type": "progress", "message": "Scanning PDF to find tables... (this can take up to a minute)"}) + "\n"
                cap_thread.join(timeout=2.0)
                
            if cap_thread.exc:
                raise cap_thread.exc
                
            captions = cap_thread.captions
            extractor_agency = extractor.agency
            extractor_analysis_type = extractor.analysis_type
            extractor.close()
            extractor = None
            
            if not captions:
                captions = [{'page': p, 'caption': f'Table Page {p}', 'is_valid': True} for p in range(1, min(30, 100))]
            
            valid_count = sum(1 for c in captions if c.get('is_valid'))
            invalid_count = len(captions) - valid_count
            
            records = []
            raw_tables = []
            
            if use_parallel:
                yield json.dumps({"type": "progress", "message": f"Extracting {valid_count} valid tables in parallel..."}) + "\n"
                import concurrent.futures
                valid_captions = [c for c in captions if c.get('is_valid')]
                processed = 0
                
                # Threads avoid Windows process-spawn/access-denied failures;
                # each worker owns and closes its own PDF extractor instance.
                with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
                    futures = {executor.submit(_extract_table_worker, temp_pdf, c['page'], c['caption'], c.get('bbox'), extractor_agency, extractor_analysis_type): c for c in valid_captions}
                    not_done = set(futures.keys())
                    while not_done:
                        done, not_done = concurrent.futures.wait(not_done, timeout=2.0, return_when=concurrent.futures.FIRST_COMPLETED)
                        if not done:
                            yield json.dumps({"type": "progress", "message": f"Still extracting tables... ({processed} of {valid_count} done)"}) + "\n"
                        else:
                            for future in done:
                                c = futures[future]
                                df, err = future.result()
                                processed += 1
                                yield json.dumps({"type": "progress", "message": f"Extracted {processed} of {valid_count} tables..."}) + "\n"
                                
                                if err:
                                    print(f"Parallel err: {err}")
                                elif df is not None and not df.empty:
                                    raw_df = df.copy()
                                    raw_df.insert(0, 'Source', f"Page {c['page']}: {c['caption']}")
                                    raw_df = raw_df.fillna("")
                                    raw_tables.append(raw_df.to_dict('records'))
                                    for _, row in df.iterrows():
                                        c_name = row.get('Compound', row.get('Compound / Identification', row.get('Compound Name', row.get('CHEMICAL', row.get('Element', '')))))
                                        cas = row.get('CASRN', row.get('CAS', ''))
                                        conc = row.get("Concentration (ug/device)", row.get("Est Conc (ug/device)", row.get("Est Conc (ug/mL)", row.get("ug_device", row.get("Amount (ug/M^3)", row.get("AREA", ""))))))
                                        if c_name or cas:
                                            records.append({'Compound Name': c_name, 'CASRN': cas, 'Concentration': conc, 'Source': f"Page {c['page']}: {c['caption']}"})
            else:
                try:
                    extractor = UniversalPDFExtractor(temp_pdf)
                    processed = 0
                    for c in captions:
                        if not c.get('is_valid'): continue
                        
                        class SequentialThread(threading.Thread):
                            def __init__(self, ext, page, caption, bbox):
                                super().__init__()
                                self.ext = ext
                                self.page = page
                                self.caption = caption
                                self.bbox = bbox
                                self.df = None
                                self.exc = None
                            def run(self):
                                try:
                                    self.df = self.ext.extract_table(self.page, self.caption, caption_bbox=self.bbox)
                                except Exception as e:
                                    self.exc = e

                        seq_thread = SequentialThread(extractor, c['page'], c['caption'], c.get('bbox'))
                        seq_thread.start()
                        
                        while seq_thread.is_alive():
                            yield json.dumps({"type": "progress", "message": f"Still extracting tables sequentially... ({processed} of {valid_count} done)"}) + "\n"
                            seq_thread.join(timeout=2.0)
                            
                        if seq_thread.exc:
                            print(f"Sequential err: {seq_thread.exc}")
                            df = None
                        else:
                            df = seq_thread.df
                            
                        processed += 1
                        yield json.dumps({"type": "progress", "message": f"Extracted {processed} of {valid_count} tables sequentially..."}) + "\n"
                        
                        if df is not None and not df.empty:
                            raw_df = df.copy()
                            raw_df.insert(0, 'Source', f"Page {c['page']}: {c['caption']}")
                            raw_df = raw_df.fillna("")
                            raw_tables.append(raw_df.to_dict('records'))
                            for _, row in df.iterrows():
                                c_name = row.get('Compound', row.get('Compound / Identification', row.get('Compound Name', row.get('CHEMICAL', row.get('Element', '')))))
                                cas = row.get('CASRN', row.get('CAS', ''))
                                conc = row.get("Concentration (ug/device)", row.get("Est Conc (ug/device)", row.get("Est Conc (ug/mL)", row.get("ug_device", row.get("Amount (ug/M^3)", row.get("AREA", ""))))))
                                if c_name or cas:
                                    records.append({'Compound Name': c_name, 'CASRN': cas, 'Concentration': conc, 'Source': f"Page {c['page']}: {c['caption']}"})
                finally:
                    if extractor is not None:
                        extractor.close()
            
            with open(f'{temp_pdf}_records.json', 'w') as f:
                json.dump(records, f)
                
            with open(f'{temp_pdf}_raw.json', 'w') as f:
                json.dump(raw_tables, f)
                
            result = {
                "total_tables": len(captions),
                "valid_tables": valid_count,
                "invalid_tables": invalid_count,
                "captions": captions,
                "raw_compounds_count": len(records)
            }
            yield json.dumps({"type": "result", "data": result}) + "\n"
            
        except Exception as e:
            if 'extractor' in locals() and extractor is not None:
                try:
                    extractor.close()
                except Exception:
                    pass
            yield json.dumps({"type": "error", "message": str(e)}) + "\n"

    from fastapi.responses import StreamingResponse
    return StreamingResponse(event_generator(), media_type="text/event-stream")

@app.post('/assess/pdf/step5')
async def assess_pdf_step5(filename: str = Form(...), use_parallel: bool = Form(True)):
    temp_pdf = _pdf_temp_path(filename)
    records_file = f'{temp_pdf}_records.json'
    
    if not os.path.exists(records_file):
        return {"error": "Session expired or records not found."}
        
    try:
        import json
        import pandas as pd
        with open(records_file, 'r') as f:
            records = json.load(f)
            
        raw_file = f'{temp_pdf}_raw.json'
        raw_tables = []
        if os.path.exists(raw_file):
            with open(raw_file, 'r') as f:
                raw_tables_data = json.load(f)
                for tbl_data in raw_tables_data:
                    raw_tables.append(pd.DataFrame(tbl_data))
                    
        df = pd.DataFrame(records)
        if df.empty:
            raise Exception("No records to process.")
            
        import re
        def clean_compound(name):
            if not isinstance(name, str): return name
            name = re.sub(r'(?i)\brelated\s+component\b', '', name)
            name = re.sub(r'(?i)\brelated\s+compou?nd\s*\d*\b', '', name)
            name = re.sub(r'(?i)\s*\(n=\d+\)', '', name)
            name = re.sub(r'(?i),\s*loss\s+of\s+.*$', '', name)
            name = re.sub(r'(?i),\s*(hydroxy|methyl|oxy)$', '', name)
            name = re.sub(r'(?i)\s+dimer$', '', name)
            name = re.sub(r',\s*$', '', name)
            return re.sub(r'\s+', ' ', name).strip()
        
        def parse_conc(val):
            if pd.isna(val) or val is None: return -1.0
            s = str(val).replace(',', '').strip()
            match = re.search(r'[\d.]+', s)
            if match:
                try: return float(match.group(0))
                except: return -1.0
            return -1.0
            
        df['Cleaned Compound Name'] = df['Compound Name'].apply(clean_compound)
        df['Conc_Float'] = df['Concentration'].apply(parse_conc)
        
        valid_df = df[df['Cleaned Compound Name'] != ''].copy()
        if not valid_df.empty:
            idx = valid_df.groupby('Compound Name')['Conc_Float'].idxmax()
            df = valid_df.loc[idx].copy()
            
        df = df.fillna('')
        
        return {
            "unique_compounds_count": len(df),
            "compounds": df.to_dict(orient='records')
        }
    except Exception as e:
        return {"error": str(e)}
    finally:
        if os.path.exists(temp_pdf): os.remove(temp_pdf)
        if os.path.exists(records_file): os.remove(records_file)

from fastapi import Request
from fastapi.responses import JSONResponse

from fastapi.responses import StreamingResponse
import json

@app.post('/assess/pdf/resolve_identity')
async def assess_pdf_resolve_identity(request: Request):
    try:
        data = await request.json()
    except Exception as e:
        return JSONResponse(status_code=400, content={"error": f"Invalid JSON body: {str(e)}"})

    async def event_generator(data_dict):
        try:
            filename = data_dict.get('filename', 'report.pdf')
            compounds = data_dict.get('compounds', [])
            use_parallel = data_dict.get('use_parallel', False)
            
            if not compounds:
                yield json.dumps({"type": "error", "message": "No compounds to process"}) + "\n"
                return

            import re
            df = pd.DataFrame(compounds)
            print("Yielding progress", flush=True)
            msg = json.dumps({"type": "progress", "message": f"[{filename}] Received {len(df)} compounds. Running Phase A..."}) + "\n"
            yield msg + (" " * 2048) + "\n"
            print("Yielded progress", flush=True)

            results = []
            total = len(df)

            def _apply_pdf_identity_result(row: pd.Series, ident_res: dict, smiles: str, std_smiles: str, fgs: list, parent, c_class: str, conf, fp, struct_classification: dict, analogs: list) -> dict:
                row_dict = row.to_dict()
                structural_evidence = calculate_structural_evidence(std_smiles)
                row_dict['Resolved Name'] = ident_res.get('name', '')
                row_dict['Original SMILES'] = smiles
                row_dict['Standardized SMILES'] = std_smiles
                row_dict['InChIKey'] = ident_res.get('inchikey', '')
                row_dict['InChI'] = ident_res.get('inchi', '')
                row_dict['Chemical Class'] = c_class
                row_dict['Corrected Chemical Class'] = struct_classification['corrected_chemical_class']
                row_dict['Final Classification Record'] = struct_classification['final_classification_record']
                row_dict['Classification Input Source'] = struct_classification['classification_input_source']
                row_dict['Classification Status'] = struct_classification['classification_status']
                row_dict['Classification Rule ID'] = struct_classification['classification_rule_id']
                row_dict['Classification Rule Version'] = struct_classification['classification_rule_version']
                row_dict['Manual Review Flag'] = struct_classification['manual_review_flag']
                row_dict['Manual Review Reason'] = struct_classification['manual_review_reason']
                row_dict['Detected Feature Profile'] = "; ".join(struct_classification['detected_feature_profile'])
                exported_classification = classification_export_fields(struct_classification)
                row_dict['Primary Functional Group'] = exported_classification['Primary Functional Group']
                row_dict['Secondary Functional Groups'] = exported_classification['Secondary Functional Groups']
                row_dict['Taxonomy Path'] = struct_classification['taxonomy_path']
                row_dict['Taxonomy Path Steps'] = " → ".join(struct_classification['taxonomy_path_steps'])
                row_dict['Taxonomy Hierarchy'] = " | ".join(f"{k}: {v}" for k, v in (struct_classification.get('taxonomy_hierarchy') or {}).items())
                row_dict['Topology Profile'] = " | ".join(f"{k}: {v}" for k, v in (struct_classification.get('topology_profile') or {}).items())
                row_dict.update(structural_evidence)
                row_dict['Scaffold'] = _scaffold_smiles(std_smiles)
                row_dict.update(_resolve_chemont_bundle(ident_res.get('inchikey', ''), std_smiles))
                row_dict['Classification Rationale'] = _prism_single_assessment_contract(struct_classification)['classification_rationale']
                row_dict.update(exported_classification)
                row_dict['Confidence'] = conf
                row_dict['Tanimoto Analogs Found'] = len(analogs)
                return row_dict

            if use_parallel:
                import asyncio
                import aiohttp
                from services_async import resolve_identity_async
                from services import standardize_structure, detect_functional_groups, detect_parent_moiety, assign_chemical_class, generate_fingerprint, analog_search_and_score, determine_confidence
                from database import SessionLocal

                print("Setting up async paced queue", flush=True)
                rows = list(df.iterrows())
                connector = aiohttp.TCPConnector(limit=10)

                async with aiohttp.ClientSession(connector=connector) as session:
                    print("Entered session", flush=True)

                    semaphore = asyncio.Semaphore(15)
                    from services_async import analog_search_and_score_async, RateLimiter

                    # Maximum 4.5 requests per second to stay safely under 5 req/sec limit
                    pubchem_rate_limiter = RateLimiter(4.5)

                    async def process_compound(idx, row):
                        async with semaphore:
                            c_name = str(row.get('Cleaned Compound Name', '')).strip()
                            if c_name == 'nan':
                                c_name = ''
                            c_cas = str(row.get('CASRN', '')).strip()
                            if c_cas == 'nan' or not re.match(r'^\d{2,7}-\d{2}-\d$', c_cas):
                                c_cas = ''

                            req = schemas.AssessmentRequest(compound_name=c_name, cas_number=c_cas)
                            if not c_name and not c_cas:
                                return (idx, None)

                            try:
                                ident_res = await resolve_identity_async(req, session, rate_limiter=pubchem_rate_limiter)
                                if ident_res is None or isinstance(ident_res, Exception):
                                    row_dict = row.to_dict()
                                    row_dict['Resolved Name'] = 'Error'
                                    row_dict['Chemical Class'] = str(ident_res) if isinstance(ident_res, Exception) else 'Empty'
                                    row_dict['Corrected Chemical Class'] = 'Unclassified'
                                    row_dict['Final Classification Record'] = 'Structure unavailable -> Manual Review'
                                    row_dict['Classification Input Source'] = 'Unavailable'
                                    row_dict['Classification Status'] = 'Manual review'
                                    row_dict['Classification Rule ID'] = 'STC-UNK-000'
                                    row_dict['Classification Rule Version'] = 'unknown'
                                    row_dict['Manual Review Flag'] = True
                                    row_dict['Manual Review Reason'] = 'Identity resolution did not return a valid structure.'
                                    row_dict['Detected Feature Profile'] = ''
                                    row_dict['Primary Functional Group'] = ''
                                    row_dict['Secondary Functional Groups'] = ''
                                    row_dict['Taxonomy Path'] = 'Unclassified'
                                    row_dict['Taxonomy Path Steps'] = 'Unclassified'
                                    row_dict['Topology Profile'] = 'Topology Class: Unclassified | Topology Modifiers: Unclassified | Ring System: Unclassified | Ring Profile: Not assessable'
                                    row_dict['Classification Rationale'] = 'Not available'
                                    row_dict.update(_resolve_chemont_bundle('', ''))
                                    return (idx, row_dict)

                                smiles = ident_res.get('smiles', '')

                                db_local = SessionLocal()
                                try:
                                    def cpu_sync():
                                        std_smiles = standardize_structure(smiles)
                                        fgs = detect_functional_groups(std_smiles)
                                        parent = detect_parent_moiety(std_smiles)
                                        c_class = assign_chemical_class(parent, fgs)
                                        conf = determine_confidence(ident_res.get('status'), std_smiles, parent, fgs)
                                        fp = generate_fingerprint(std_smiles)
                                        struct_class = classify_structure(std_smiles, smiles)
                                        return std_smiles, fgs, parent, c_class, conf, fp, struct_class

                                    std_smiles, fgs, parent, c_class, conf, fp, struct_classification = await asyncio.to_thread(cpu_sync)

                                    analogs = await analog_search_and_score_async(std_smiles, fp, parent, fgs, db_local, session, rate_limiter=pubchem_rate_limiter)

                                    row_dict = _apply_pdf_identity_result(row, ident_res, smiles, std_smiles, fgs, parent, c_class, conf, fp, struct_classification, analogs)
                                    return (idx, row_dict)
                                finally:
                                    db_local.close()

                            except Exception as e:
                                row_dict = row.to_dict()
                                row_dict['Resolved Name'] = 'Error'
                                row_dict['Chemical Class'] = str(e)
                                row_dict['Corrected Chemical Class'] = 'Unclassified'
                                row_dict['Final Classification Record'] = 'Structure unavailable -> Manual Review'
                                row_dict['Classification Input Source'] = 'Unavailable'
                                row_dict['Classification Status'] = 'Manual review'
                                row_dict['Classification Rule ID'] = 'STC-UNK-000'
                                row_dict['Classification Rule Version'] = 'unknown'
                                row_dict['Manual Review Flag'] = True
                                row_dict['Manual Review Reason'] = str(e)
                                row_dict['Detected Feature Profile'] = ''
                                row_dict['Primary Functional Group'] = ''
                                row_dict['Secondary Functional Groups'] = ''
                                row_dict['Taxonomy Path'] = 'Unclassified'
                                row_dict['Taxonomy Path Steps'] = 'Unclassified'
                                row_dict['Taxonomy Hierarchy'] = 'Parent Class: Unclassified | Functional Group: Unclassified | Subclass: Unclassified | Structural Type: Unclassified'
                                row_dict['Topology Profile'] = 'Topology Class: Unclassified | Topology Modifiers: Unclassified | Ring System: Unclassified | Ring Profile: Not assessable'
                                row_dict.update(_resolve_chemont_bundle('', ''))
                                row_dict['Classification Rationale'] = 'Not available'
                                return (idx, row_dict)

                    tasks = []
                    for idx, row in enumerate(rows):
                        tasks.append(asyncio.create_task(process_compound(idx, row[1])))

                    results = [None] * total
                    processed_count = 0

                    for f in asyncio.as_completed(tasks):
                        idx, row_result = await f
                        results[idx] = row_result

                        processed_count += 1
                        if processed_count % 3 == 0 or processed_count == total:
                            yield json.dumps({"type": "progress", "message": f"[{filename}] Processed {processed_count} of {total} compounds..."}) + "\n"

                    results = [r for r in results if r is not None]

            else:
                from concurrent.futures import ThreadPoolExecutor
                from database import SessionLocal

                def process_row(args):
                    i, row_tup = args
                    _, row = row_tup
                    c_name = str(row.get('Cleaned Compound Name', '')).strip()
                    if c_name == 'nan':
                        c_name = ''
                    c_cas = str(row.get('CASRN', '')).strip()
                    if c_cas == 'nan' or not re.match(r'^\d{2,7}-\d{2}-\d$', c_cas):
                        c_cas = ''
                    if not c_name and not c_cas:
                        return None
                    req = schemas.AssessmentRequest(compound_name=c_name, cas_number=c_cas)
                    db_local = SessionLocal()
                    try:
                        res = assess_chemical(req, db_local)
                        row_dict = row.to_dict()
                        row_dict.update(res.structural_evidence)
                        row_dict['Resolved Name'] = res.identity.get('name', '')
                        row_dict['Original SMILES'] = res.structure.get('original_smiles', '')
                        row_dict['Standardized SMILES'] = res.structure.get('standardized_smiles', '')
                        row_dict['InChIKey'] = res.identity.get('inchikey', '')
                        row_dict['InChI'] = res.identity.get('inchi', '')
                        row_dict['Chemical Class'] = res.chemical_class
                        row_dict['Corrected Chemical Class'] = res.corrected_chemical_class
                        row_dict['Final Classification Record'] = res.final_classification_record
                        row_dict['Classification Input Source'] = res.classification_input_source
                        row_dict['Classification Status'] = res.classification_status
                        row_dict['Classification Rule ID'] = res.classification_rule_id
                        row_dict['Classification Rule Version'] = res.classification_rule_version
                        row_dict['Classification Scope'] = res.classification_scope
                        row_dict['Classification Review Recommended'] = res.classification_review_recommended
                        row_dict['Classification Review Recommendation'] = res.classification_review_recommendation
                        row_dict['Manual Review Flag'] = res.manual_review_flag
                        row_dict['Manual Review Reason'] = res.manual_review_reason
                        row_dict['Detected Feature Profile'] = "; ".join(res.detected_feature_profile)
                        row_dict['Primary Functional Group'] = res.primary_functional_group
                        row_dict['Secondary Functional Groups'] = "; ".join(
                            fg.functional_group for fg in res.secondary_functional_groups
                        )
                        row_dict.update(classification_export_fields(_classification_response_contract(res)))
                        row_dict['Taxonomy Path'] = res.taxonomy_path
                        row_dict['Taxonomy Path Steps'] = " → ".join(res.taxonomy_path_steps)
                        row_dict['Taxonomy Hierarchy'] = " | ".join(f"{k}: {v}" for k, v in (res.taxonomy_hierarchy or {}).items())
                        row_dict['Topology Profile'] = " | ".join(f"{k}: {v}" for k, v in (res.topology_profile or {}).items())
                        row_dict['Scaffold'] = res.scaffold
                        row_dict.update(_resolve_chemont_bundle(res.identity.get('inchikey', ''), res.structure.get('standardized_smiles', '')))
                        row_dict['Classification Rationale'] = res.classification_rationale
                        row_dict['Confidence'] = res.confidence
                        row_dict['Tanimoto Analogs Found'] = len(res.analogs)
                        return row_dict
                    except Exception as e:
                        row_dict = row.to_dict()
                        row_dict['Resolved Name'] = 'Error'
                        row_dict['Chemical Class'] = str(e)
                        row_dict['Corrected Chemical Class'] = 'Unclassified'
                        row_dict['Final Classification Record'] = 'Structure unavailable -> Manual Review'
                        row_dict['Classification Input Source'] = 'Unavailable'
                        row_dict['Classification Status'] = 'Manual review'
                        row_dict['Classification Rule ID'] = 'STC-UNK-000'
                        row_dict['Classification Rule Version'] = 'unknown'
                        row_dict['Manual Review Flag'] = True
                        row_dict['Manual Review Reason'] = str(e)
                        row_dict['Detected Feature Profile'] = ''
                        row_dict['Primary Functional Group'] = ''
                        row_dict['Secondary Functional Groups'] = ''
                        row_dict['Taxonomy Path'] = 'Unclassified'
                        row_dict['Taxonomy Path Steps'] = 'Unclassified'
                        row_dict['Taxonomy Hierarchy'] = 'Parent Class: Unclassified | Functional Group: Unclassified | Subclass: Unclassified | Structural Type: Unclassified'
                        row_dict['Topology Profile'] = 'Topology Class: Unclassified | Topology Modifiers: Unclassified | Ring System: Unclassified | Ring Profile: Not assessable'
                        row_dict.update(_resolve_chemont_bundle('', ''))
                        row_dict['Classification Rationale'] = 'Not available'
                        return row_dict
                    finally:
                        db_local.close()

                with ThreadPoolExecutor(max_workers=10) as executor:
                    for idx, res_dict in enumerate(executor.map(process_row, enumerate(df.iterrows(), 1)), 1):
                        if res_dict:
                            results.append(res_dict)
                        if idx % 10 == 0 or idx == total:
                            yield json.dumps({"type": "progress", "message": f"[{filename}] Processed {idx} of {total} compounds..."}) + "\n"

            out_df = pd.DataFrame(results)
            import identity_rescue
            out_df = identity_rescue.rescue_unresolved_compounds(out_df)
            
            import math
            def clean_nan(obj):
                if isinstance(obj, dict):
                    return {k: clean_nan(v) for k, v in obj.items()}
                elif isinstance(obj, list):
                    return [clean_nan(v) for v in obj]
                elif isinstance(obj, float) and math.isnan(obj):
                    return None
                return obj
                
            out_df_json = out_df.fillna('').to_dict(orient='records')
            out_df_json = clean_nan(out_df_json)
            
            yield json.dumps({"type": "result", "filename": filename, "compounds": out_df_json}, allow_nan=False) + "\n"
            
        except Exception as e:
            import traceback
            trace = traceback.format_exc()
            print(f"Resolve identity error: {e}\n{trace}", flush=True)
            yield json.dumps({"error": str(e), "trace": trace}) + "\n"
            
    return StreamingResponse(event_generator(data), media_type="application/x-ndjson")


@app.post('/evidence/import')
async def import_endpoint_evidence(request: Request):
    """Controlled manual JSON ingestion for Stage 4 local evidence records."""
    from database import SessionLocal
    from evidence_ledger import EvidenceValidationError, ingest_evidence_records
    payload = await request.json()
    records = payload.get("records", [])
    imported_by = payload.get("imported_by", "")
    if not isinstance(records, list) or not records:
        raise HTTPException(status_code=400, detail="records must be a non-empty list.")
    session = SessionLocal()
    try:
        return {"outcomes": ingest_evidence_records(session, records, imported_by)}
    except EvidenceValidationError as exc:
        session.rollback()
        raise HTTPException(status_code=400, detail=str(exc))
    finally:
        session.close()


@app.post('/evidence/import/csv')
async def import_endpoint_evidence_csv(file: UploadFile = File(...), imported_by: str = Form(...)):
    """Controlled CSV ingestion; records are validated before local retention."""
    import csv
    from io import StringIO
    from database import SessionLocal
    from evidence_ledger import EvidenceValidationError, ingest_evidence_records
    try:
        content = (await file.read()).decode("utf-8-sig")
        records = list(csv.DictReader(StringIO(content)))
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Evidence CSV must be UTF-8 encoded.")
    session = SessionLocal()
    try:
        return {"outcomes": ingest_evidence_records(session, records, imported_by)}
    except EvidenceValidationError as exc:
        session.rollback()
        raise HTTPException(status_code=400, detail=str(exc))
    finally:
        session.close()


@app.get('/evidence/export/csv')
async def export_endpoint_evidence_csv(identity_key: str | None = None):
    """Export immutable Stage 4 evidence/provenance records for audit review."""
    from database import SessionLocal
    from evidence_ledger import export_evidence_csv
    session = SessionLocal()
    try:
        content = export_evidence_csv(session, [identity_key] if identity_key else None)
    finally:
        session.close()
    return StreamingResponse(iter([content]), media_type="text/csv", headers={"Content-Disposition": "attachment; filename=endpoint_evidence_ledger.csv"})


@app.post('/assess/pdf/run_ai_clustering')
async def assess_pdf_run_ai_clustering(request: Request):
    try:
        data = await request.json()
        filename = data.get('filename', 'report.pdf')
        compounds = data.get('compounds', [])
        
        if not compounds:
            raise Exception('No compounds to process')
            
        out_df = pd.DataFrame(compounds)
        print(f"[{filename}] Running Phase B: AI Clustering on {len(out_df)} compounds...", flush=True)
        out_df = attach_chemont_classifications(out_df)
        
        import clustering_engine
        out_df = clustering_engine.process_clustering_and_adme(out_df)
        
        from io import BytesIO
        out_buffer = BytesIO()
        import excel_formatter

        raw_file = f'{_pdf_temp_path(filename)}_raw.json'
        raw_tables = []
        if os.path.exists(raw_file):
            import json
            with open(raw_file, 'r') as f:
                raw_tables_data = json.load(f)
                for tbl_data in raw_tables_data:
                    raw_tables.append(pd.DataFrame(tbl_data))
            os.remove(raw_file)
            
        excel_formatter.format_excel_output(out_df, out_buffer, raw_tables=raw_tables)
        out_buffer.seek(0)
        
        import base64
        excel_base64 = base64.b64encode(out_buffer.getvalue()).decode('utf-8')
        
        out_df_json = out_df.fillna('').to_dict(orient='records')
        
        return JSONResponse(content={
            "filename": f'processed_{filename}',
            "excel_base64": excel_base64,
            "compounds": out_df_json
        })
    except Exception as e:
        import traceback
        trace = traceback.format_exc()
        print(f"Phase B error: {trace}", flush=True)
        raise HTTPException(status_code=500, detail=str(e))


@app.post('/assess/pdf/run_deterministic_enrichment')
async def assess_pdf_run_deterministic_enrichment(request: Request):
    """Run deterministic Phase B enrichment and export the workbook."""
    phase_b_trace = []
    try:
        data = await request.json()
        filename = data.get('filename', 'report.pdf')
        compounds = data.get('compounds', [])
        run_id = data.get('run_id') or _phase_b_run_id()
        stage = str(data.get('phase_b_stage', 'all')).strip().lower()
        # A refresh/retry may send only the run_id. Recover the last committed
        # compound snapshot for every post-Phase-A stage.
        if not compounds and stage not in {"a", "phase_a", "phasea", "snapshot"}:
            with _PHASE_B_RUNS_LOCK:
                saved_run = _PHASE_B_RUNS.get(run_id)
                compounds = list(saved_run.get("compounds", [])) if saved_run else []
        if not compounds:
            raise HTTPException(status_code=400, detail="No compounds available for workbook export.")

        def _normalise_workbook_name(source_filename: str) -> str:
            base_name = os.path.basename(str(source_filename or 'report.pdf'))
            stem = os.path.splitext(base_name)[0].strip() or 'report'
            for prefix in ('phase_a_', 'processed_'):
                while stem.lower().startswith(prefix):
                    stem = stem[len(prefix):].strip() or 'report'
            return f"{stem}.xlsx"

        workbook_name = _normalise_workbook_name(filename)

        out_df = pd.DataFrame(compounds)
        # Reapply the current deterministic contract from structure before any
        # B-stage clustering or export. This makes resume/retry independent of
        # stale fields carried by a prior Phase A workbook.
        out_df = apply_classification_contract(out_df)
        # Restore the extraction snapshot once and persist it with this local
        # run so both Phase A snapshots and the final B2 workbook include the
        # complete pre-deduplication audit sheet.
        raw_extraction_tables = _load_raw_extraction_tables(filename)
        if raw_extraction_tables:
            out_df.attrs["raw_extraction_tables"] = raw_extraction_tables
        phase_b_trace.append({
            "step": "input_loaded",
            "status": "ok",
            "rows": int(len(out_df)),
            "columns": int(len(out_df.columns)),
            "note": "Request payload converted to DataFrame and classification contract refreshed"
        })
        if stage in {"a", "phase_a", "phasea", "snapshot"}:
            if not out_df.empty:
                phase_b_trace.append({"step": "workbook_contract", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Applying workbook contract for Phase A snapshot"})
                out_df = apply_workbook_contract(out_df)
                out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
                phase_b_trace.append({"step": "workbook_contract", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Workbook contract attached"})
            out_buffer = BytesIO()
            import excel_formatter
            phase_b_trace.append({"step": "excel_format", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Formatting Phase A snapshot workbook"})
            excel_formatter.format_excel_output(
                out_df, out_buffer, raw_tables=out_df.attrs.get("raw_extraction_tables", [])
            )
            phase_b_trace.append({"step": "excel_format", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Phase A snapshot workbook formatted"})
            out_buffer.seek(0)

            import base64
            phase_b_trace.append({"step": "base64_encode", "status": "start", "note": "Encoding Phase A snapshot workbook bytes"})
            excel_base64 = base64.b64encode(out_buffer.getvalue()).decode('utf-8')
            phase_b_trace.append({"step": "base64_encode", "status": "ok", "note": "Phase A snapshot workbook encoded"})
            out_df_json = _phase_b_plain_records(out_df)
            run_id = _capture_phase_b_state(run_id, out_df, "A", phase_b_trace, filename)
            return JSONResponse(content={
                "filename": f'phase_a_{workbook_name}',
                "source_filename": str(filename or 'report.pdf'),
                "excel_base64": excel_base64,
                "compounds": out_df_json,
                "phase_b_trace": phase_b_trace,
                "phase_b_stage": "A",
                "run_id": run_id,
            })

        if stage in {"b1", "b1a", "core", "enrich", "enrichment"}:
            if not out_df.empty:
                previous_xai_enabled = os.environ.get("XAI_ENABLED")
                previous_xai_stage6_enabled = os.environ.get("XAI_STAGE6_ENABLED")
                previous_loky_max_cpu = os.environ.get("LOKY_MAX_CPU_COUNT")
                os.environ["XAI_ENABLED"] = "false"
                os.environ["XAI_STAGE6_ENABLED"] = "false"
                os.environ["LOKY_MAX_CPU_COUNT"] = "1"
                import clustering_engine
                try:
                    phase_b_trace.append({"step": "phase_b_engine", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Entering Phase B1a deterministic core engine"})
                    out_df = clustering_engine.process_clustering_and_adme_core(out_df, debug_trace=phase_b_trace)
                    phase_b_trace.append({"step": "phase_b_engine", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Phase B1a core engine returned enriched DataFrame"})
                finally:
                    if previous_xai_enabled is None:
                        os.environ.pop("XAI_ENABLED", None)
                    else:
                        os.environ["XAI_ENABLED"] = previous_xai_enabled
                    if previous_xai_stage6_enabled is None:
                        os.environ.pop("XAI_STAGE6_ENABLED", None)
                    else:
                        os.environ["XAI_STAGE6_ENABLED"] = previous_xai_stage6_enabled
                    if previous_loky_max_cpu is None:
                        os.environ.pop("LOKY_MAX_CPU_COUNT", None)
                    else:
                        os.environ["LOKY_MAX_CPU_COUNT"] = previous_loky_max_cpu
            if not out_df.empty:
                phase_b_trace.append({"step": "workbook_contract", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Applying workbook contract"})
                # Capture engine attrs before the contract transformation.
                if raw_extraction_tables:
                    out_df.attrs["raw_extraction_tables"] = raw_extraction_tables
                run_id = _capture_phase_b_state(run_id, out_df, "B1A", phase_b_trace, filename)
                out_df = apply_workbook_contract(out_df)
                out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
                phase_b_trace.append({"step": "workbook_contract", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Workbook contract attached"})
            out_df_json = _phase_b_plain_records(out_df)
            phase_b_trace.append({"step": "phase_b_stage", "status": "ok", "note": "B1a core enrichment complete"})
            return JSONResponse(content={
                "filename": f'processed_{workbook_name}',
                "source_filename": str(filename or 'report.pdf'),
                "compounds": out_df_json,
                "phase_b_trace": phase_b_trace,
                "phase_b_stage": "B1A",
                "next_stage": "B1B",
                "run_id": run_id,
            })

        if stage in {"b1b1", "b1b_1", "chemont", "b1b2", "b1b_2", "epa", "epa_ctx", "b1b3", "b1b_3", "chembl", "b1b4", "b1b_4"}:
            substage = {
                "b1b1": "chemont", "b1b_1": "chemont", "chemont": "chemont",
                "b1b2": "epa_ctx", "b1b_2": "epa_ctx", "epa": "epa_ctx", "epa_ctx": "epa_ctx",
                "b1b3": "chembl", "b1b_3": "chembl", "chembl": "chembl",
                "b1b4": "finalize", "b1b_4": "finalize",
            }[stage]
            next_stage = {"chemont": "B1B-2", "epa_ctx": "B1B-3", "chembl": "B1B-4", "finalize": "B2"}[substage]
            if not out_df.empty:
                previous_xai_enabled = os.environ.get("XAI_ENABLED")
                previous_xai_stage6_enabled = os.environ.get("XAI_STAGE6_ENABLED")
                previous_loky_max_cpu = os.environ.get("LOKY_MAX_CPU_COUNT")
                os.environ["XAI_ENABLED"] = "false"
                os.environ["XAI_STAGE6_ENABLED"] = "false"
                os.environ["LOKY_MAX_CPU_COUNT"] = "1"
                import clustering_engine
                try:
                    phase_b_trace.append({"step": "phase_b_engine", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": f"Entering Phase B1b substage {substage}"})
                    out_df = clustering_engine.process_clustering_and_adme_optional_stage(out_df, substage, debug_trace=phase_b_trace)
                    phase_b_trace.append({"step": "phase_b_engine", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": f"Phase B1b substage {substage} returned enriched DataFrame"})
                finally:
                    if previous_xai_enabled is None:
                        os.environ.pop("XAI_ENABLED", None)
                    else:
                        os.environ["XAI_ENABLED"] = previous_xai_enabled
                    if previous_xai_stage6_enabled is None:
                        os.environ.pop("XAI_STAGE6_ENABLED", None)
                    else:
                        os.environ["XAI_STAGE6_ENABLED"] = previous_xai_stage6_enabled
                    if previous_loky_max_cpu is None:
                        os.environ.pop("LOKY_MAX_CPU_COUNT", None)
                    else:
                        os.environ["LOKY_MAX_CPU_COUNT"] = previous_loky_max_cpu
            if not out_df.empty:
                if raw_extraction_tables:
                    out_df.attrs["raw_extraction_tables"] = raw_extraction_tables
                run_id = _capture_phase_b_state(run_id, out_df, f"B1B-{substage}", phase_b_trace, filename)
                out_df = apply_workbook_contract(out_df)
                out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
            out_df_json = _phase_b_plain_records(out_df)
            phase_b_trace.append({"step": "phase_b_stage", "status": "ok", "note": f"B1b substage {substage} complete"})
            return JSONResponse(content={
                "filename": f"processed_{workbook_name}",
                "source_filename": str(filename or "report.pdf"),
                "compounds": out_df_json,
                "phase_b_trace": phase_b_trace,
                "phase_b_stage": f"B1B-{ {'chemont':'1','epa_ctx':'2','chembl':'3','finalize':'4'}[substage] }",
                "next_stage": next_stage,
                "run_id": run_id,
            })

        if stage in {"b1b", "optional", "finalize"}:
            if not out_df.empty:
                previous_xai_enabled = os.environ.get("XAI_ENABLED")
                previous_xai_stage6_enabled = os.environ.get("XAI_STAGE6_ENABLED")
                previous_loky_max_cpu = os.environ.get("LOKY_MAX_CPU_COUNT")
                os.environ["XAI_ENABLED"] = "false"
                os.environ["XAI_STAGE6_ENABLED"] = "false"
                os.environ["LOKY_MAX_CPU_COUNT"] = "1"
                import clustering_engine
                try:
                    phase_b_trace.append({"step": "phase_b_engine", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Entering Phase B1b optional enrichment engine"})
                    out_df = clustering_engine.process_clustering_and_adme_optional(out_df, debug_trace=phase_b_trace)
                    phase_b_trace.append({"step": "phase_b_engine", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Phase B1b optional engine returned enriched DataFrame"})
                finally:
                    if previous_xai_enabled is None:
                        os.environ.pop("XAI_ENABLED", None)
                    else:
                        os.environ["XAI_ENABLED"] = previous_xai_enabled
                    if previous_xai_stage6_enabled is None:
                        os.environ.pop("XAI_STAGE6_ENABLED", None)
                    else:
                        os.environ["XAI_STAGE6_ENABLED"] = previous_xai_stage6_enabled
                    if previous_loky_max_cpu is None:
                        os.environ.pop("LOKY_MAX_CPU_COUNT", None)
                    else:
                        os.environ["LOKY_MAX_CPU_COUNT"] = previous_loky_max_cpu
            if not out_df.empty:
                if raw_extraction_tables:
                    out_df.attrs["raw_extraction_tables"] = raw_extraction_tables
                run_id = _capture_phase_b_state(run_id, out_df, "B1B", phase_b_trace, filename)
                phase_b_trace.append({"step": "workbook_contract", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Applying workbook contract"})
                out_df = apply_workbook_contract(out_df)
                out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
                phase_b_trace.append({"step": "workbook_contract", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Workbook contract attached"})
            out_df_json = _phase_b_plain_records(out_df)
            phase_b_trace.append({"step": "phase_b_stage", "status": "ok", "note": "B1b optional enrichment complete"})
            return JSONResponse(content={
                "filename": f'processed_{workbook_name}',
                "source_filename": str(filename or 'report.pdf'),
                "compounds": out_df_json,
                "phase_b_trace": phase_b_trace,
                "phase_b_stage": "B1B",
                "next_stage": "B2",
                "run_id": run_id,
            })

        if stage in {"b2", "export"}:
            # Reattach all evidence bundles accumulated by B1a/B1b stages.
            _restore_phase_b_state(run_id, out_df)
            validation_errors = _validate_phase_b_export(out_df, run_id)
            phase_b_trace.append({"step": "b2_validation", "status": "error" if validation_errors else "ok", "errors": validation_errors, "note": "Pre-export structural and prerequisite validation"})
            if validation_errors:
                return JSONResponse(status_code=409, content={
                    "detail": "B2 export blocked by validation gate.",
                    "validation_errors": validation_errors,
                    "phase_b_trace": phase_b_trace,
                    "phase_b_stage": "B2",
                    "run_id": run_id,
                })
            if not out_df.empty:
                out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
            out_buffer = BytesIO()
            import excel_formatter
            phase_b_trace.append({"step": "excel_format", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Formatting workbook"})
            excel_formatter.format_excel_output(
                out_df, out_buffer, raw_tables=out_df.attrs.get("raw_extraction_tables", [])
            )
            phase_b_trace.append({"step": "excel_format", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Workbook formatting complete"})
            out_buffer.seek(0)

            import base64
            phase_b_trace.append({"step": "base64_encode", "status": "start", "note": "Encoding workbook bytes"})
            excel_base64 = base64.b64encode(out_buffer.getvalue()).decode('utf-8')
            phase_b_trace.append({"step": "base64_encode", "status": "ok", "note": "Workbook encoded"})
            out_df_json = _phase_b_plain_records(out_df)
            return JSONResponse(content={
                "filename": f'processed_{workbook_name}',
                "source_filename": str(filename or 'report.pdf'),
                "excel_base64": excel_base64,
                "compounds": out_df_json,
                "phase_b_trace": phase_b_trace,
                "phase_b_stage": "B2",
                "run_id": run_id,
            })

        if not out_df.empty:
            previous_xai_enabled = os.environ.get("XAI_ENABLED")
            previous_xai_stage6_enabled = os.environ.get("XAI_STAGE6_ENABLED")
            previous_loky_max_cpu = os.environ.get("LOKY_MAX_CPU_COUNT")
            os.environ["XAI_ENABLED"] = "false"
            os.environ["XAI_STAGE6_ENABLED"] = "false"
            os.environ["LOKY_MAX_CPU_COUNT"] = "1"
            import clustering_engine
            try:
                phase_b_trace.append({"step": "phase_b_engine", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Entering deterministic Phase B engine"})
                out_df = clustering_engine.process_clustering_and_adme(out_df, debug_trace=phase_b_trace)
                phase_b_trace.append({"step": "phase_b_engine", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Phase B engine returned enriched DataFrame"})
            finally:
                if previous_xai_enabled is None:
                    os.environ.pop("XAI_ENABLED", None)
                else:
                    os.environ["XAI_ENABLED"] = previous_xai_enabled
                if previous_xai_stage6_enabled is None:
                    os.environ.pop("XAI_STAGE6_ENABLED", None)
                else:
                    os.environ["XAI_STAGE6_ENABLED"] = previous_xai_stage6_enabled
                if previous_loky_max_cpu is None:
                    os.environ.pop("LOKY_MAX_CPU_COUNT", None)
                else:
                    os.environ["LOKY_MAX_CPU_COUNT"] = previous_loky_max_cpu
        if not out_df.empty:
            phase_b_trace.append({"step": "workbook_contract", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Applying workbook contract"})
            out_df = apply_workbook_contract(out_df)
            out_df.attrs["workbook_contract"] = build_workbook_contract_frame()
            phase_b_trace.append({"step": "workbook_contract", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Workbook contract attached"})

        out_buffer = BytesIO()
        import excel_formatter
        if raw_extraction_tables:
            out_df.attrs["raw_extraction_tables"] = raw_extraction_tables
        phase_b_trace.append({"step": "excel_format", "status": "start", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Formatting workbook"})
        excel_formatter.format_excel_output(
            out_df, out_buffer, raw_tables=out_df.attrs.get("raw_extraction_tables", [])
        )
        phase_b_trace.append({"step": "excel_format", "status": "ok", "rows": int(len(out_df)), "columns": int(len(out_df.columns)), "note": "Workbook formatting complete"})
        out_buffer.seek(0)

        import base64
        phase_b_trace.append({"step": "base64_encode", "status": "start", "note": "Encoding workbook bytes"})
        excel_base64 = base64.b64encode(out_buffer.getvalue()).decode('utf-8')
        phase_b_trace.append({"step": "base64_encode", "status": "ok", "note": "Workbook encoded"})
        out_df_json = out_df.fillna('').to_dict(orient='records')

        return JSONResponse(content={
            "filename": f'processed_{workbook_name}',
            "source_filename": str(filename or 'report.pdf'),
            "excel_base64": excel_base64,
            "compounds": out_df_json,
            "phase_b_trace": phase_b_trace,
            "phase_b_stage": "ALL"
        })
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        trace = traceback.format_exc()
        phase_b_trace.append({"step": "error", "status": "error", "note": str(e)})
        print(f"Phase A workbook export error: {trace}", flush=True)
        return JSONResponse(status_code=500, content={
            "detail": str(e),
            "phase_b_trace": phase_b_trace,
            "traceback": trace,
        })

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8055)
