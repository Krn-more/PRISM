"""Offline exact-match lookup for the Chemical Taxonomies reference source.

The database is deliberately an evidence layer: an exact, structure-specific
ClassyFire taxonomy can enrich PRISM's ontology context, while PRISM's local
SMARTS rules remain the authority for functional-group interpretation and
review decisions.
"""

from __future__ import annotations

import os
import sqlite3
import sys
import tempfile
from pathlib import Path
from typing import Any


DATABASE_VERSION = "ChemicalTaxonomies v1.1"
SIDECAR_FILENAME = "chemical_taxonomies_v1_1.sqlite"
# The explicit PRISM_TAXONOMY_DB setting is preferred for a deployed,
# persistent sidecar.  The default must nevertheless be writable for the
# portable beta application, including corporate installations where the app
# directory and LocalAppData root can be read-only to the launched process.
DEFAULT_DATABASE_PATH = Path(tempfile.gettempdir()) / "PRISM" / "taxonomy" / "chemical_taxonomies_v1_1.sqlite"
HIERARCHY_COLUMNS = ("Kingdom", "Superclass", "Class", "Subclass", "Level_5", "Level_6", "Level_7", "Level_8", "Level_9", "Level_10")


def configured_database_path() -> Path:
    """Return the deployment override, portable sidecar, or local cache path."""
    override = os.environ.get("PRISM_TAXONOMY_DB", "").strip()
    if override:
        return Path(override).expanduser()

    # Beta distribution is deliberately portable: when the frozen executable
    # is unzipped beside the prepared SQLite reference file, no environment
    # variable or installer step is required.  Do not return a missing path;
    # the writable cache remains the safe fallback for existing deployments.
    if getattr(sys, "frozen", False):
        portable_sidecar = Path(sys.executable).resolve().parent / SIDECAR_FILENAME
        if portable_sidecar.is_file():
            return portable_sidecar

    return DEFAULT_DATABASE_PATH


def normalise_inchikey(value: object) -> str:
    return str(value or "").strip().upper()


def database_available(db_path: str | Path | None = None) -> bool:
    path = Path(db_path) if db_path else configured_database_path()
    return path.is_file()


def lookup_by_inchikey(inchikey: object, db_path: str | Path | None = None) -> dict[str, str] | None:
    """Look up one exact InChIKey without loading the reference database."""
    key = normalise_inchikey(inchikey)
    path = Path(db_path) if db_path else configured_database_path()
    if not key or not path.is_file():
        return None
    try:
        with sqlite3.connect(f"file:{path.resolve().as_posix()}?mode=ro", uri=True) as connection:
            connection.row_factory = sqlite3.Row
            row = connection.execute(
                """
                SELECT record_id, chemical_name, casrn, dtxsid, inchikey, smiles,
                       kingdom, superclass, class, subclass, level_5, level_6,
                       level_7, level_8, level_9, level_10, data_version,
                       date_extracted, reference, url
                FROM taxonomy_records
                WHERE inchikey = ?
                ORDER BY record_id
                LIMIT 1
                """,
                (key,),
            ).fetchone()
    except (sqlite3.Error, OSError):
        # A missing/corrupt optional evidence store must never block PRISM's
        # deterministic classifier.
        return None
    return dict(row) if row is not None else None


def taxonomy_path(record: dict[str, Any]) -> list[str]:
    """Return the populated reference hierarchy from broad to specific."""
    return [str(record.get(column.lower(), "")).strip() for column in HIERARCHY_COLUMNS if str(record.get(column.lower(), "")).strip() and str(record.get(column.lower(), "")).strip().upper() != "NA"]


def ontology_evidence(inchikey: object, db_path: str | Path | None = None) -> dict[str, str] | None:
    """Adapt an exact local taxonomy match to PRISM's ontology evidence shape."""
    record = lookup_by_inchikey(inchikey, db_path=db_path)
    if record is None:
        return None
    path = taxonomy_path(record)
    direct_parent = path[-1] if path else ""
    version = str(record.get("data_version") or DATABASE_VERSION)
    reference = str(record.get("reference") or "Chemical Taxonomies")
    source_url = str(record.get("url") or "")
    evidence = (
        f"Exact local Chemical Taxonomies match by InChIKey {record['inchikey']}; "
        f"record {record.get('record_id', '')}; reference {reference}."
    )
    if source_url:
        evidence += f" Source: {source_url}."
    return {
        "mapping_status": "Exact local taxonomy match",
        "ontology_mapping_version": version,
        "kingdom": str(record.get("kingdom") or ""),
        "superclass": str(record.get("superclass") or ""),
        "class": str(record.get("class") or ""),
        "subclass": str(record.get("subclass") or ""),
        "direct_parent": direct_parent,
        "mapping_evidence": evidence,
        "source": "Chemical Taxonomies local SQLite",
        "record_id": str(record.get("record_id") or ""),
        "taxonomy_path": " → ".join(path),
        "reference": reference,
        "url": source_url,
    }
