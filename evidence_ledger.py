"""Stage 4 local, provenance-preserving endpoint-evidence ledger."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from sqlalchemy.orm import Session

from models import EndpointEvidenceRecord


DEFAULT_VOCABULARY_PATH = Path(__file__).with_name("config") / "stage4_evidence_vocabularies.json"
REQUIRED_FIELDS = (
    "identity_key", "compound_mapping_method", "endpoint", "result", "study_type",
    "species", "route", "source", "source_record", "reliability", "evidence_status",
)
EXPORT_COLUMNS = (
    "evidence_id", "identity_key", "compound_mapping_method", "endpoint", "result",
    "study_type", "species", "route", "dose", "dose_unit", "duration", "source",
    "source_record", "source_url", "source_database_version", "reliability",
    "evidence_status", "raw_record_hash", "vocabulary_version", "retrieved_at",
    "imported_at", "imported_by", "reviewer_status", "supersedes_evidence_id", "notes",
)


class EvidenceValidationError(ValueError):
    """Raised when an evidence record cannot be safely retained in the ledger."""


def load_vocabularies(path: str | Path | None = None) -> dict[str, Any]:
    with open(path or DEFAULT_VOCABULARY_PATH, encoding="utf-8") as handle:
        vocabularies = json.load(handle)
    required = {"vocabulary_version", "endpoint_families", "results", "study_types", "species", "routes", "reliability", "evidence_statuses", "reviewer_statuses", "mapping_methods"}
    missing = required.difference(vocabularies)
    if missing:
        raise EvidenceValidationError("Controlled vocabulary is missing: " + ", ".join(sorted(missing)))
    return vocabularies


def _text(value: object) -> str:
    if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
        return ""
    return str(value).strip()


def _canonical_payload(record: dict[str, Any]) -> dict[str, str]:
    return {field: _text(record.get(field)) for field in EXPORT_COLUMNS if field not in {"evidence_id", "raw_record_hash", "imported_at", "imported_by", "vocabulary_version"}}


def _hash_record(record: dict[str, Any]) -> str:
    serialized = json.dumps(_canonical_payload(record), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _validate_dose(dose: str, dose_unit: str) -> None:
    if not dose:
        return
    plain_numeric = re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)", dose)
    if plain_numeric:
        numeric = float(dose)
        if not math.isfinite(numeric) or numeric < 0:
            raise EvidenceValidationError("Dose must be a finite non-negative value.")
        if not dose_unit:
            raise EvidenceValidationError("A plain numeric dose requires dose_unit; retain qualitative/range/concentration text in dose.")
    if re.search(r"(?:^|\s)-\d", dose):
        raise EvidenceValidationError("Dose range contains a negative value.")


def validate_evidence_record(record: dict[str, Any], vocabularies: dict[str, Any] | None = None) -> dict[str, str]:
    """Validate and normalise one record without inferring a hazard conclusion."""
    vocabularies = vocabularies or load_vocabularies()
    normalized = {field: _text(record.get(field)) for field in EXPORT_COLUMNS}
    for field in REQUIRED_FIELDS:
        if not normalized[field]:
            raise EvidenceValidationError(f"Missing required evidence field: {field}")
    controlled = {
        "compound_mapping_method": "mapping_methods", "endpoint": "endpoint_families", "result": "results",
        "study_type": "study_types", "species": "species", "route": "routes",
        "reliability": "reliability", "evidence_status": "evidence_statuses",
    }
    for field, vocabulary in controlled.items():
        if normalized[field] not in vocabularies[vocabulary]:
            raise EvidenceValidationError(f"Invalid {field}: {normalized[field]!r}. Allowed values are controlled by {vocabulary}.")
    reviewer_status = normalized["reviewer_status"] or "Not reviewed"
    if reviewer_status not in vocabularies["reviewer_statuses"]:
        raise EvidenceValidationError(f"Invalid reviewer_status: {reviewer_status!r}.")
    normalized["reviewer_status"] = reviewer_status
    _validate_dose(normalized["dose"], normalized["dose_unit"])
    if normalized["evidence_status"] == "No record retrieved" and normalized["result"] == "Negative":
        raise EvidenceValidationError("No record retrieved is not a negative study; use Not applicable or Not reported as result.")
    if not normalized["source_url"] and normalized["evidence_status"] == "Retrieved":
        normalized["evidence_status"] = "Incomplete"
        suffix = "Source URL not supplied; provenance is incomplete."
        normalized["notes"] = (normalized["notes"] + " " + suffix).strip()
    normalized["vocabulary_version"] = vocabularies["vocabulary_version"]
    normalized["raw_record_hash"] = _hash_record(normalized)
    normalized["evidence_id"] = normalized["evidence_id"] or "EV-" + normalized["raw_record_hash"][:16].upper()
    return normalized


def ingest_evidence_records(session: Session, records: Iterable[dict[str, Any]], imported_by: str, imported_at: str | None = None) -> list[dict[str, str]]:
    """Idempotently retain validated records; corrections must supersede, never overwrite."""
    if not _text(imported_by):
        raise EvidenceValidationError("imported_by is required for audit provenance.")
    timestamp = imported_at or datetime.now(timezone.utc).isoformat()
    vocabularies = load_vocabularies()
    outcomes: list[dict[str, str]] = []
    try:
        for source_record in records:
            normalized = validate_evidence_record(source_record, vocabularies)
            existing = session.query(EndpointEvidenceRecord).filter_by(raw_record_hash=normalized["raw_record_hash"]).one_or_none()
            if existing:
                outcomes.append({"evidence_id": existing.evidence_id, "status": "Duplicate", "raw_record_hash": existing.raw_record_hash})
                continue
            if session.query(EndpointEvidenceRecord).filter_by(evidence_id=normalized["evidence_id"]).one_or_none():
                raise EvidenceValidationError(f"evidence_id {normalized['evidence_id']} already exists with different content; create a new record that supersedes it.")
            prior_records = session.query(EndpointEvidenceRecord).filter_by(
                identity_key=normalized["identity_key"], source=normalized["source"], source_record=normalized["source_record"]
            ).all()
            supersedes = normalized["supersedes_evidence_id"]
            if prior_records and not supersedes:
                raise EvidenceValidationError("Changed content for an existing source record requires supersedes_evidence_id; original evidence is retained unchanged.")
            if supersedes:
                superseded = session.query(EndpointEvidenceRecord).filter_by(evidence_id=supersedes).one_or_none()
                if not superseded:
                    raise EvidenceValidationError(f"supersedes_evidence_id {supersedes!r} does not exist.")
                if (superseded.identity_key, superseded.source, superseded.source_record) != (normalized["identity_key"], normalized["source"], normalized["source_record"]):
                    raise EvidenceValidationError("supersedes_evidence_id must refer to the same identity_key, source, and source_record.")
            normalized["imported_at"] = timestamp
            normalized["imported_by"] = _text(imported_by)
            session.add(EndpointEvidenceRecord(**normalized))
            outcomes.append({"evidence_id": normalized["evidence_id"], "status": "Imported", "raw_record_hash": normalized["raw_record_hash"]})
        session.commit()
        return outcomes
    except Exception:
        session.rollback()
        raise


def import_evidence_csv(session: Session, csv_path: str | Path, imported_by: str) -> list[dict[str, str]]:
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:
        return ingest_evidence_records(session, csv.DictReader(handle), imported_by)


def export_evidence_dataframe(session: Session, identity_keys: Iterable[str] | None = None) -> pd.DataFrame:
    query = session.query(EndpointEvidenceRecord)
    keys = sorted({_text(key) for key in identity_keys or [] if _text(key)})
    if keys:
        query = query.filter(EndpointEvidenceRecord.identity_key.in_(keys))
    records = query.order_by(EndpointEvidenceRecord.imported_at, EndpointEvidenceRecord.evidence_id).all()
    return pd.DataFrame([{field: getattr(record, field) for field in EXPORT_COLUMNS} for record in records], columns=EXPORT_COLUMNS)


def export_evidence_csv(session: Session, identity_keys: Iterable[str] | None = None) -> str:
    """Return a complete, reproducible CSV audit export including provenance fields."""
    return export_evidence_dataframe(session, identity_keys).to_csv(index=False)


def attach_evidence_ledger(df: pd.DataFrame, session_factory=None) -> pd.DataFrame:
    """Attach report-relevant local evidence as an Excel-only dataframe attribute."""
    from database import SessionLocal
    keys = []
    # Stage 5 can resolve a supplied CAS/InChIKey to DTXSID. Include that
    # resolved public identifier so the same report receives its EPA-backed
    # ledger rows even when its original input did not contain DTXSID.
    for column in ("InChIKey", "CAS", "CASRN", "DTXSID", "EPA_CTX_DTXSID"):
        if column in df.columns:
            keys.extend(df[column].tolist())
    session = (session_factory or SessionLocal)()
    try:
        evidence = export_evidence_dataframe(session, keys)
    finally:
        session.close()
    result = df.copy()
    result.attrs["endpoint_evidence_ledger"] = evidence
    return result
