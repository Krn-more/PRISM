"""Stage 5 read-only EPA CTX adapter with DTXSID cache and conservative ledger mapping."""

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
from dotenv import load_dotenv
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from evidence_ledger import ingest_evidence_records
from models import EndpointEvidenceRecord, EpaCtxCacheEntry


DEFAULT_RULES_PATH = Path(__file__).with_name("config") / "stage5_epa_ctx_rules.json"
_REQUEST_LOCK = threading.Lock()
_LAST_REQUEST_MONOTONIC = 0.0


class EpaCtxError(RuntimeError):
    pass


def load_epa_ctx_rules(path: str | Path | None = None) -> dict[str, Any]:
    with open(path or DEFAULT_RULES_PATH, encoding="utf-8") as handle:
        rules = json.load(handle)
    required = {"rule_version", "base_url", "chemical_equal_path", "hazard_by_dtxsid_path", "source_version", "timeout_seconds", "max_retries", "min_request_interval_seconds", "cache_ttl_days", "allowed_identifier_columns"}
    missing = required.difference(rules)
    if missing:
        raise EpaCtxError("EPA CTX rules missing: " + ", ".join(sorted(missing)))
    return rules


def _load_environment() -> None:
    load_dotenv(Path(__file__).with_name(".env"))


def epa_ctx_enabled() -> bool:
    _load_environment()
    return os.getenv("EPA_CTX_ENABLED", "false").strip().casefold() == "true" and bool(os.getenv("EPA_CTX_API_KEY", "").strip())


def _api_key() -> str:
    _load_environment()
    key = os.getenv("EPA_CTX_API_KEY", "").strip()
    if not key:
        raise EpaCtxError("EPA CTX API key is not configured.")
    return key


def _session(rules: dict[str, Any]) -> requests.Session:
    session = requests.Session()
    session.trust_env = False
    retry = Retry(total=int(rules["max_retries"]), backoff_factor=0.5, status_forcelist=(429, 500, 502, 503, 504), allowed_methods=("POST",))
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def _canonical_hash(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def _records(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        for value in payload.values():
            if isinstance(value, list) and all(isinstance(item, dict) for item in value):
                return value
        return [payload]
    return []


def _extract_dtxsid(payload: Any) -> str | None:
    for record in _records(payload):
        for key, value in record.items():
            if str(key).casefold() == "dtxsid" and str(value).upper().startswith("DTXSID"):
                return str(value).upper()
    return None


def _respect_request_interval(seconds: float) -> None:
    """Serialise outgoing CTX calls at a conservative configurable pace."""
    global _LAST_REQUEST_MONOTONIC
    with _REQUEST_LOCK:
        wait_seconds = max(0.0, _LAST_REQUEST_MONOTONIC + max(0.0, seconds) - time.monotonic())
        if wait_seconds:
            time.sleep(wait_seconds)
        _LAST_REQUEST_MONOTONIC = time.monotonic()


def _request_json(path: str, payload: list[str], rules: dict[str, Any], post: Callable[..., Any] | None = None) -> Any:
    requester = post or _session(rules).post
    _respect_request_interval(float(rules["min_request_interval_seconds"]))
    response = requester(
        rules["base_url"].rstrip("/") + path,
        headers={"x-api-key": _api_key(), "Accept": "application/json"},
        json=payload,
        timeout=(10, int(rules["timeout_seconds"])),
    )
    if response.status_code >= 400:
        raise EpaCtxError(f"EPA CTX request failed with HTTP {response.status_code}.")
    try:
        return response.json()
    except ValueError as exc:
        raise EpaCtxError("EPA CTX returned a non-JSON response.") from exc


def _request_json_get(path: str, rules: dict[str, Any], get: Callable[..., Any] | None = None) -> Any:
    requester = get or _session(rules).get
    _respect_request_interval(float(rules["min_request_interval_seconds"]))
    url = rules["base_url"].rstrip("/") + path
    response = requester(
        url,
        headers={"x-api-key": _api_key(), "Accept": "application/json"},
        timeout=(10, int(rules["timeout_seconds"])),
    )
    if response.status_code >= 400:
        raise EpaCtxError(f"EPA CTX request failed with HTTP {response.status_code}.")
    try:
        return response.json()
    except ValueError as exc:
        raise EpaCtxError("EPA CTX returned a non-JSON response.") from exc


def resolve_dtxsid(identifier: str, rules: dict[str, Any] | None = None, post: Callable[..., Any] | None = None) -> str | None:
    """Resolve an approved public identifier through EPA Chemical equality search."""
    identifier = str(identifier or "").strip()
    if not identifier:
        return None
    if identifier.upper().startswith("DTXSID"):
        return identifier.upper()
    rules = rules or load_epa_ctx_rules()
    encoded_identifier = requests.utils.quote(identifier, safe="")
    return _extract_dtxsid(_request_json_get(rules["chemical_equal_path"] + encoded_identifier, rules, post))


def fetch_hazard_by_dtxsid(dtxsid: str, db_session, rules: dict[str, Any] | None = None, post: Callable[..., Any] | None = None, now: datetime | None = None) -> tuple[Any, str]:
    """Return cached or live EPA hazard payload. Cache entries contain public data only."""
    rules = rules or load_epa_ctx_rules()
    dtxsid = str(dtxsid).upper()
    now = now or datetime.now(timezone.utc)
    existing = db_session.query(EpaCtxCacheEntry).filter_by(dtxsid=dtxsid, source_version=rules["source_version"], resource_path=rules["hazard_by_dtxsid_path"]).one_or_none()
    if existing:
        retrieved = datetime.fromisoformat(existing.retrieved_at.replace("Z", "+00:00"))
        if retrieved + timedelta(days=int(rules["cache_ttl_days"])) >= now:
            return existing.response_json, "Cached"
    payload = _request_json_get(rules["hazard_by_dtxsid_path"] + requests.utils.quote(dtxsid, safe=""), rules, post)
    record = {
        "dtxsid": dtxsid, "source_version": rules["source_version"], "resource_path": rules["hazard_by_dtxsid_path"],
        "retrieved_at": now.isoformat(), "http_status": 200, "response_hash": _canonical_hash(payload), "response_json": payload,
    }
    if existing:
        for key, value in record.items():
            setattr(existing, key, value)
    else:
        db_session.add(EpaCtxCacheEntry(**record))
    db_session.commit()
    return payload, "Retrieved"


def _value(record: dict[str, Any], *names: str) -> str:
    lookup = {str(key).casefold(): value for key, value in record.items()}
    for name in names:
        value = lookup.get(name.casefold())
        if value not in (None, ""):
            return str(value).strip()
    return ""


def _endpoint_family(record: dict[str, Any]) -> str | None:
    evidence_text = " ".join(_value(record, "effect", "effectName", "endpoint", "toxvalType", "studyType", "title", "description").casefold().split())
    mappings = (("Genotoxicity", ("genotox", "mutagen", "genetic")), ("Sensitisation", ("sensiti", "sensitiz")), ("Reproductive/developmental toxicity", ("repro", "development", "fertility", "terato")), ("Carcinogenicity", ("carcin", "cancer")))
    for family, terms in mappings:
        if any(term in evidence_text for term in terms):
            return family
    return None


def map_hazard_payload_to_ledger(dtxsid: str, payload: Any, rules: dict[str, Any] | None = None) -> list[dict[str, str]]:
    """Map only unambiguous endpoint families; unclassified EPA records remain cached raw."""
    rules = rules or load_epa_ctx_rules()
    ledger_records: list[dict[str, str]] = []
    for record in _records(payload):
        endpoint = _endpoint_family(record)
        if not endpoint:
            continue
        numeric = _value(record, "toxvalNumeric", "value")
        qualifier = _value(record, "toxvalNumericQualifier", "qualifier")
        dose = " ".join(part for part in (qualifier, numeric) if part)
        source_record = _value(record, "id", "recordId", "sourceRecord") or "EPA-CTX-" + _canonical_hash(record)[:16]
        species = _value(record, "species") or "Not reported"
        route = _value(record, "exposureRoute", "route") or "Not reported"
        study_type = _value(record, "studyType") or "Other"
        if study_type not in {"In vitro", "In vivo", "Human", "In silico", "Read-across", "Other", "Not reported"}:
            study_type = "Other"
        if species not in {"Human", "Mouse", "Rat", "Rabbit", "Guinea pig", "Other", "Not applicable", "Not reported"}:
            species = "Other"
        if route not in {"Oral", "Dermal", "Inhalation", "Parenteral", "In vitro", "Other", "Not applicable", "Not reported"}:
            route = "Other"
        ledger_records.append({
            "identity_key": dtxsid, "compound_mapping_method": "DTXSID", "endpoint": endpoint, "result": "Not reported",
            "study_type": study_type, "species": species, "route": route, "dose": dose, "dose_unit": _value(record, "toxvalUnits", "units"),
            "duration": _value(record, "studyDuration", "duration"), "source": "EPA CTX Hazard API", "source_record": source_record,
            "source_url": rules["base_url"], "source_database_version": rules["source_version"], "reliability": "Not assignable",
            "evidence_status": "Retrieved", "retrieved_at": "", "reviewer_status": "Not reviewed",
            "notes": "Mapped conservatively from EPA CTX raw record; reviewer interpretation required.",
        })
    return ledger_records


def _link_supersessions(db_session, records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Preserve prior EPA interpretations when a source record changes.

    Stage 4 is immutable: changed content for the same EPA source record must
    create a new row that names the prior row it supersedes. Exact duplicates
    remain idempotent in ``ingest_evidence_records``.
    """
    for record in records:
        existing = (
            db_session.query(EndpointEvidenceRecord)
            .filter_by(
                identity_key=record["identity_key"],
                source=record["source"],
                source_record=record["source_record"],
            )
            .order_by(EndpointEvidenceRecord.imported_at.desc(), EndpointEvidenceRecord.id.desc())
            .first()
        )
        if existing:
            record["supersedes_evidence_id"] = existing.evidence_id
    return records


def enrich_with_epa_ctx(df: pd.DataFrame, db_session, rules: dict[str, Any] | None = None, post: Callable[..., Any] | None = None) -> pd.DataFrame:
    """Append status columns and, when enabled, retain EPA evidence without changing decisions."""
    result = df.copy()
    rules = rules or load_epa_ctx_rules()
    result["EPA_CTX_Status"] = "Unavailable"
    result["EPA_CTX_DTXSID"] = ""
    result["EPA_CTX_Source_Version"] = rules["source_version"]
    if not epa_ctx_enabled():
        result["EPA_CTX_Status"] = "Unavailable: EPA_CTX_ENABLED is false or EPA_CTX_API_KEY is not configured"
        return result
    deadline = time.monotonic() + float(os.environ.get("PHASE_B_STAGE_DEADLINE_SECONDS", "300"))
    for index, row in result.iterrows():
        if time.monotonic() >= deadline:
            result.loc[index:, "EPA_CTX_Status"] = "Unavailable: StageDeadlineExceeded"
            break
        identifier = next((str(row[column]).strip() for column in rules["allowed_identifier_columns"] if column in result.columns and pd.notna(row[column]) and str(row[column]).strip()), "")
        if not identifier:
            result.at[index, "EPA_CTX_Status"] = "Not requested: approved public identifier unavailable"
            continue
        try:
            dtxsid = resolve_dtxsid(identifier, rules, post)
            if not dtxsid:
                result.at[index, "EPA_CTX_Status"] = "No DTXSID resolved"
                continue
            payload, cache_status = fetch_hazard_by_dtxsid(dtxsid, db_session, rules, post)
            mapped = map_hazard_payload_to_ledger(dtxsid, payload, rules)
            if mapped:
                ingest_evidence_records(db_session, _link_supersessions(db_session, mapped), imported_by="EPA_CTX")
            result.at[index, "EPA_CTX_DTXSID"] = dtxsid
            result.at[index, "EPA_CTX_Status"] = f"{cache_status}: {len(mapped)} Stage 4-mappable record(s)"
        except Exception as exc:
            result.at[index, "EPA_CTX_Status"] = f"Unavailable: {type(exc).__name__}"
    return result
