"""Stage 6 deterministic cluster evidence matrix and grounded-output validation."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_RULES_PATH = Path(__file__).with_name("config") / "stage6_xai_rules.json"
MATRIX_COLUMNS = ("Cluster_ID", "Evidence_ID", "Uncertainty_ID", "Evidence_Domain", "Evidence_Type", "Evidence_Text", "Source_Fields", "Matrix_Version")


class GroundedOutputError(ValueError):
    """Raised when an LLM response cannot be tied to the supplied matrix."""


def load_stage6_rules(path: str | Path | None = None) -> dict[str, Any]:
    with open(path or DEFAULT_RULES_PATH, encoding="utf-8") as handle:
        rules = json.load(handle)
    required = {"matrix_version", "prompt_version", "max_source_text_length", "max_rendered_characters", "required_output_fields", "banned_claim_patterns", "endpoint_claim_terms"}
    if missing := required.difference(rules):
        raise GroundedOutputError("Stage 6 rules missing: " + ", ".join(sorted(missing)))
    return rules


def _canonical_text(value: Any) -> str:
    if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
        return ""
    return " ".join(str(value).replace("\x00", " ").split())


def _text(value: Any, limit: int) -> str:
    return _canonical_text(value)[:limit]


def _column_series(frame: pd.DataFrame, column: str) -> pd.Series | None:
    if column not in frame.columns:
        return None
    subset = frame.loc[:, frame.columns == column]
    if isinstance(subset, pd.Series):
        return subset
    if subset.shape[1] == 1:
        return subset.iloc[:, 0]
    return subset.apply(
        lambda row: next((value for value in row.tolist() if not pd.isna(value) and str(value).strip()),
                         next((value for value in row.tolist() if not pd.isna(value)), "")),
        axis=1,
    )


def _stable_id(prefix: str, cluster_id: str, domain: str, text: str) -> str:
    digest = hashlib.sha256(f"{cluster_id}|{domain}|{text}".encode("utf-8")).hexdigest()[:12].upper()
    return f"{prefix}-{digest}"


def _unique_values(frame: pd.DataFrame, column: str, limit: int) -> str:
    series = _column_series(frame, column)
    if series is None:
        return "Not supplied"
    # Preserve full canonical source text for the evidence-ID hash. Display
    # truncation occurs only in add(), so distinct long source statements do
    # not collapse into the same trace identifier.
    values = sorted({_canonical_text(value) for value in series.tolist() if _canonical_text(value)})
    return "; ".join(values) if values else "Not supplied"


def _numeric_ranges(frame: pd.DataFrame) -> str:
    pieces = []
    for column in ("MW", "LogP", "TPSA", "HBD", "HBA", "Rotatable_Bonds", "Formal_Charge", "Ring_Count", "Aromatic_Ring_Count", "Fraction_CSP3"):
        series = _column_series(frame, column)
        if series is None:
            continue
        values = pd.to_numeric(series, errors="coerce").dropna()
        if not values.empty:
            pieces.append(f"{column}={values.min():.2f}-{values.max():.2f}")
    return "; ".join(pieces) or "Not supplied"


def build_cluster_evidence_matrix(df: pd.DataFrame, uncertainty_register: pd.DataFrame | None = None, endpoint_evidence_ledger: pd.DataFrame | None = None, cluster_col: str = "Cluster ID", rules_path: str | Path | None = None) -> pd.DataFrame:
    """Create deterministic, review-only evidence rows; no chemical decision is changed."""
    rules = load_stage6_rules(rules_path)
    limit = int(rules["max_source_text_length"])
    records: list[dict[str, str]] = []
    if cluster_col not in df.columns:
        return pd.DataFrame(columns=MATRIX_COLUMNS)
    for cluster_value, frame in df.groupby(cluster_col, dropna=False, sort=True):
        cluster_id = _text(cluster_value, limit) or "Unassigned"
        if cluster_id == "Unclustered":
            continue
        def add(domain: str, evidence_type: str, text: str, source_fields: str, uncertainty: bool = False) -> None:
            canonical_text = _canonical_text(text)
            text = canonical_text[:limit]
            if not text or text == "Not supplied":
                return
            evidence_id = _stable_id("E", cluster_id, domain + "|" + evidence_type, canonical_text)
            records.append({"Cluster_ID": cluster_id, "Evidence_ID": evidence_id, "Uncertainty_ID": _stable_id("U", cluster_id, domain, canonical_text) if uncertainty else "", "Evidence_Domain": domain, "Evidence_Type": evidence_type, "Evidence_Text": text, "Source_Fields": source_fields, "Matrix_Version": rules["matrix_version"]})
        add("Grouping", "Deterministic cluster", f"Members={len(frame)}; Cluster_Status={_unique_values(frame, 'Cluster_Status', limit)}; Decision={_unique_values(frame, 'Decision', limit)}", "Cluster ID; Cluster_Status; Decision")
        add("Structure", "Scaffold", _unique_values(frame, "Scaffold", limit), "Scaffold")
        add("Structure", "Toxicophore profile", _unique_values(frame, "Toxicophore_Profile", limit), "Toxicophore_Profile")
        add("Structure", "ChemOnt hierarchy", "; ".join(
            f"{label}={_unique_values(frame, column, limit)}"
            for label, column in (("Superclass", "ChemOnt_Superclass"), ("Class", "ChemOnt_Class"), ("Subclass", "ChemOnt_Subclass"), ("Direct Parent", "ChemOnt_Direct_Parent"))
            if column in frame.columns
        ) or "Not supplied", "ChemOnt_Superclass; ChemOnt_Class; ChemOnt_Subclass; ChemOnt_Direct_Parent", True)
        add("Reactivity", "Structural alerts", _unique_values(frame, "All_Structural_Alerts", limit), "All_Structural_Alerts")
        add("Physicochemical", "Calculated property ranges", _numeric_ranges(frame), "MW; LogP; TPSA; HBD; HBA; Rotatable_Bonds; Formal_Charge; Ring_Count; Aromatic_Ring_Count; Fraction_CSP3")
        add("Domain", "Domain status", _unique_values(frame, "Domain_Status", limit), "Domain_Status", True)
        add("Domain", "Domain reasons", _unique_values(frame, "Domain_Reasons", limit), "Domain_Reasons", True)
        add("Identity", "Identity status", _unique_values(frame, "Identity Status", limit), "Identity Status", True)
        add("Physicochemical", "Member property-outlier flags", _unique_values(frame, "Property_Outlier_Flags", limit), "Property_Outlier_Flags", True)
        add("Biological", "Available targets", _unique_values(frame, "ChEMBL_Targets", limit), "ChEMBL_Targets", True)
        if isinstance(uncertainty_register, pd.DataFrame) and not uncertainty_register.empty:
            scope_series = _column_series(uncertainty_register, "Scope")
            scope_id_series = _column_series(uncertainty_register, "Scope_ID")
            if scope_series is None or scope_id_series is None:
                rows = pd.DataFrame()
            else:
                rows = uncertainty_register[(scope_series == "Cluster") & (scope_id_series.astype(str) == cluster_id)]
            for _, row in rows.sort_values(["Category", "Level", "Reason"], kind="stable").iterrows():
                add("Uncertainty", f"{_text(row.get('Category'), limit)} ({_text(row.get('Level'), limit)})", f"{_text(row.get('Reason'), limit)} Required action: {_text(row.get('Required_Action'), limit)}", _text(row.get("Evidence_Fields"), limit), True)
        if isinstance(endpoint_evidence_ledger, pd.DataFrame) and not endpoint_evidence_ledger.empty:
            member_keys = set()
            for column in ("DTXSID", "EPA_CTX_DTXSID", "InChIKey", "CAS", "CASRN"):
                series = _column_series(frame, column)
                if series is not None:
                    member_keys.update(_text(value, limit) for value in series if _text(value, limit))
            identity_key_series = _column_series(endpoint_evidence_ledger, "identity_key")
            rows = endpoint_evidence_ledger[identity_key_series.astype(str).isin(member_keys)] if identity_key_series is not None else pd.DataFrame()
            for _, row in rows.sort_values(["endpoint", "source_record"], kind="stable").iterrows():
                add("Biological", "Retrieved endpoint evidence", f"Endpoint={_text(row.get('endpoint'), limit)}; Result={_text(row.get('result'), limit)}; Source={_text(row.get('source'), limit)}; Record={_text(row.get('source_record'), limit)}; Review={_text(row.get('reviewer_status'), limit)}", "Endpoint Evidence Ledger", True)
    matrix = pd.DataFrame(records, columns=MATRIX_COLUMNS)
    return matrix.sort_values(["Cluster_ID", "Evidence_Domain", "Evidence_Type", "Evidence_ID"], kind="stable").reset_index(drop=True) if not matrix.empty else matrix


def matrix_hash(matrix: pd.DataFrame, cluster_id: str) -> str:
    rows = matrix[matrix["Cluster_ID"].astype(str) == str(cluster_id)].to_dict(orient="records")
    return hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def validate_grounded_output(content: str, matrix: pd.DataFrame, cluster_id: str, rules_path: str | Path | None = None) -> dict[str, Any]:
    rules = load_stage6_rules(rules_path)
    try:
        payload = json.loads(content)
    except (TypeError, json.JSONDecodeError) as exc:
        raise GroundedOutputError("LLM response is not valid JSON.") from exc
    if not isinstance(payload, dict) or any(field not in payload for field in rules["required_output_fields"]):
        raise GroundedOutputError("LLM response does not contain the required grounded-output fields.")
    scoped = matrix[matrix["Cluster_ID"].astype(str) == str(cluster_id)]
    known_evidence = set(scoped["Evidence_ID"].dropna())
    known_uncertainty = set(scoped["Uncertainty_ID"].dropna()) - {""}
    evidence_ids = payload["evidence_ids"]
    uncertainty_ids = payload["uncertainty_ids"]
    if not isinstance(evidence_ids, list) or not evidence_ids or not isinstance(uncertainty_ids, list) or not set(evidence_ids).issubset(known_evidence) or not set(uncertainty_ids).issubset(known_uncertainty):
        raise GroundedOutputError("LLM response references evidence not present in the supplied matrix.")
    narrative_fields = ("cluster_grouping_rationale", "assessment_boundary", "uncertainty")
    claim_evidence = payload["claim_evidence"]
    if not isinstance(claim_evidence, dict) or set(claim_evidence) != set(narrative_fields):
        raise GroundedOutputError("LLM response must link every narrative field to supplied evidence IDs.")
    linked_ids: set[str] = set()
    for field in narrative_fields:
        if not isinstance(payload[field], str) or not payload[field].strip() or not isinstance(claim_evidence[field], list) or not claim_evidence[field] or not set(claim_evidence[field]).issubset(known_evidence):
            raise GroundedOutputError("Each narrative field requires non-empty, valid supporting evidence IDs.")
        linked_ids.update(claim_evidence[field])
    if set(evidence_ids) != linked_ids:
        raise GroundedOutputError("Top-level evidence_ids must exactly match the claim-level evidence links.")
    rendered = "\n\n".join(f"{label}: {payload[field].strip()}" for label, field in (("Cluster Grouping Rationale", "cluster_grouping_rationale"), ("Assessment Boundary", "assessment_boundary"), ("Uncertainty", "uncertainty")))
    if len(rendered) > int(rules["max_rendered_characters"]):
        raise GroundedOutputError("LLM response exceeds the approved rendered length.")
    if any(re.search(pattern, rendered, flags=re.IGNORECASE) for pattern in rules["banned_claim_patterns"]):
        raise GroundedOutputError("LLM response contains a prohibited safety or regulatory claim.")
    for term in rules["endpoint_claim_terms"]:
        if term in rendered.casefold():
            raise GroundedOutputError(f"LLM response contains a prohibited endpoint conclusion term: {term}.")
    payload["rendered"] = rendered
    return payload
