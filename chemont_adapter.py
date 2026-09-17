"""Stage 7 optional ClassyFire/ChemOnt structural-taxonomy enrichment.

The adapter is deliberately cache-first and review-only.  A ClassyFire result
is a structural taxonomy label, not evidence of toxicological equivalence,
endpoint outcome, safety, or regulatory acceptability.
"""

from __future__ import annotations

import hashlib
import json
import os
import threading
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import requests
import os
import time
from dotenv import load_dotenv

from models import ChemontClassificationCacheEntry


DEFAULT_RULES_PATH = Path(__file__).with_name("config") / "stage7_chemont_rules.json"
_REQUEST_LOCK = threading.Lock()
_LAST_REQUEST_MONOTONIC = 0.0
OUTPUT_COLUMNS = (
    "ChemOnt_Kingdom", "ChemOnt_Superclass", "ChemOnt_Class", "ChemOnt_Subclass",
    "ChemOnt_Direct_Parent", "ChemOnt_Molecular_Framework", "ChemOnt_Substituents",
    "ChemOnt_Classification_Version", "ChemOnt_Source", "ChemOnt_Retrieval_Status",
)


class ChemontError(RuntimeError):
    pass


def load_chemont_rules(path: str | Path | None = None) -> dict[str, Any]:
    with open(path or DEFAULT_RULES_PATH, encoding="utf-8") as handle:
        rules = json.load(handle)
    required = {"rule_version", "source", "base_url", "entity_path_template", "timeout_seconds", "cache_ttl_days", "min_request_interval_seconds", "allowed_identifier_columns", "required_output_fields", "status_version"}
    missing = required.difference(rules)
    if missing:
        raise ChemontError("ChemOnt rules missing: " + ", ".join(sorted(missing)))
    if rules["allowed_identifier_columns"] != ["InChIKey"]:
        raise ChemontError("ChemOnt requests may use only the public InChIKey identifier.")
    return rules


def _load_environment() -> None:
    load_dotenv(Path(__file__).with_name(".env"))


def chemont_enabled() -> bool:
    _load_environment()
    return os.getenv("CHEMONT_CLASSYFIRE_ENABLED", "false").strip().casefold() == "true"


def _canonical_hash(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def _node_name(value: Any) -> str:
    if isinstance(value, dict):
        return str(value.get("name", "")).strip()
    return str(value).strip() if isinstance(value, str) else ""


def normalise_classification(payload: Any, rules: dict[str, Any] | None = None) -> dict[str, str]:
    """Reduce a raw ClassyFire response to the approved report fields."""
    rules = rules or load_chemont_rules()
    if not isinstance(payload, dict):
        raise ChemontError("ClassyFire returned a non-object JSON response.")
    values = {
        "ChemOnt_Kingdom": _node_name(payload.get("kingdom")),
        "ChemOnt_Superclass": _node_name(payload.get("superclass")),
        "ChemOnt_Class": _node_name(payload.get("class")),
        "ChemOnt_Subclass": _node_name(payload.get("subclass")),
        "ChemOnt_Direct_Parent": _node_name(payload.get("direct_parent")),
        "ChemOnt_Molecular_Framework": _node_name(payload.get("molecular_framework")),
        "ChemOnt_Substituents": "; ".join(sorted({_node_name(item) for item in payload.get("substituents", []) if _node_name(item)})),
        "ChemOnt_Classification_Version": str(payload.get("classification_version") or payload.get("version") or "Not reported"),
    }
    if not any(values[key] for key in ("ChemOnt_Kingdom", "ChemOnt_Superclass", "ChemOnt_Class", "ChemOnt_Direct_Parent")):
        raise ChemontError("ClassyFire response contains no classification hierarchy.")
    return values


def _respect_request_interval(seconds: float) -> None:
    global _LAST_REQUEST_MONOTONIC
    with _REQUEST_LOCK:
        delay = max(0.0, _LAST_REQUEST_MONOTONIC + max(0.0, seconds) - time.monotonic())
        if delay:
            time.sleep(delay)
        _LAST_REQUEST_MONOTONIC = time.monotonic()


def _request_classification(inchikey: str, rules: dict[str, Any], get: Callable[..., Any] | None = None) -> Any:
    _respect_request_interval(float(rules["min_request_interval_seconds"]))
    requester = get or requests.get
    path = rules["entity_path_template"].format(inchikey=inchikey)
    response = requester(rules["base_url"].rstrip("/") + path, headers={"Accept": "application/json"}, timeout=(5, int(rules["timeout_seconds"])))
    if response.status_code == 404:
        return None
    if response.status_code >= 400:
        raise ChemontError(f"ClassyFire request failed with HTTP {response.status_code}.")
    try:
        return response.json()
    except ValueError as exc:
        raise ChemontError("ClassyFire returned a non-JSON response.") from exc


def fetch_chemont_by_inchikey(inchikey: str, db_session, rules: dict[str, Any] | None = None, get: Callable[..., Any] | None = None, now: datetime | None = None) -> tuple[Any | None, str]:
    """Return a cached response or retrieve one public InChIKey classification."""
    rules = rules or load_chemont_rules()
    inchikey = str(inchikey or "").strip().upper()
    if not inchikey:
        return None, "Not requested: public InChIKey unavailable"
    now = now or datetime.now(timezone.utc)
    existing = db_session.query(ChemontClassificationCacheEntry).filter_by(inchikey=inchikey, source=rules["source"]).one_or_none()
    if existing:
        retrieved = datetime.fromisoformat(existing.retrieved_at.replace("Z", "+00:00"))
        if retrieved + timedelta(days=int(rules["cache_ttl_days"])) >= now:
            return existing.response_json, "Cached"
    payload = _request_classification(inchikey, rules, get)
    if payload is None:
        return None, "No classification returned"
    record = {
        "inchikey": inchikey, "source": rules["source"], "standardized_smiles": str(payload.get("smiles", "")) if isinstance(payload, dict) else "",
        "classification_version": str(payload.get("classification_version") or payload.get("version") or "Not reported") if isinstance(payload, dict) else "Not reported",
        "retrieved_at": now.isoformat(), "http_status": 200, "response_hash": _canonical_hash(payload), "response_json": payload,
    }
    if existing:
        for key, value in record.items():
            setattr(existing, key, value)
    else:
        db_session.add(ChemontClassificationCacheEntry(**record))
    db_session.commit()
    return payload, "Retrieved"


def enrich_with_chemont(df: pd.DataFrame, db_session, rules: dict[str, Any] | None = None, get: Callable[..., Any] | None = None) -> pd.DataFrame:
    """Append hierarchy/status columns without changing clustering or decisions."""
    result = df.copy()
    rules = rules or load_chemont_rules()
    for column in OUTPUT_COLUMNS:
        result[column] = ""
    result["ChemOnt_Source"] = rules["source"]
    if not chemont_enabled():
        result["ChemOnt_Retrieval_Status"] = "Unavailable: CHEMONT_CLASSYFIRE_ENABLED is false"
        result.attrs["chemont_classifications"] = result[list(OUTPUT_COLUMNS)].copy()
        return result
    deadline = time.monotonic() + float(os.environ.get("PHASE_B_STAGE_DEADLINE_SECONDS", "300"))
    for index, row in result.iterrows():
        if time.monotonic() >= deadline:
            result.loc[index:, "ChemOnt_Retrieval_Status"] = "Unavailable: StageDeadlineExceeded"
            break
        inchikey = str(row.get("InChIKey", "") or "").strip()
        try:
            payload, status = fetch_chemont_by_inchikey(inchikey, db_session, rules, get)
            if payload is not None:
                result_values = normalise_classification(payload, rules)
                for key, value in result_values.items():
                    result.at[index, key] = value
            result.at[index, "ChemOnt_Retrieval_Status"] = status
        except Exception as exc:
            result.at[index, "ChemOnt_Retrieval_Status"] = f"Unavailable: {type(exc).__name__}"
    result.attrs["chemont_classifications"] = result[[column for column in result.columns if column in OUTPUT_COLUMNS or column in ("InChIKey", "Standardized SMILES")]].copy()
    return result
