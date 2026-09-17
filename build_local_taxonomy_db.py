"""Build PRISM's indexed local Chemical Taxonomies SQLite sidecar.

Usage:
  python build_local_taxonomy_db.py --source ChemicalTaxonomies_v1.1.txt
"""

from __future__ import annotations

import argparse
import csv
import os
import sqlite3
from pathlib import Path

from local_taxonomy_db import DEFAULT_DATABASE_PATH, normalise_inchikey


EXPECTED_COLUMNS = [
    "Data_Name", "Data_Type", "Data_Version", "Record_ID", "Chemical_Name", "CASRN", "DTXSID", "InChiKey", "SMILES",
    "Kingdom", "Superclass", "Class", "Subclass", "Level_5", "Level_6", "Level_7", "Level_8", "Level_9", "Level_10",
    "Structure_Used", "Tool_Used", "Date_Extracted", "Reference", "URL",
]


def _value(row: dict[str, str], column: str) -> str:
    return str(row.get(column, "") or "").strip()


def build_database(source: str | Path, output: str | Path = DEFAULT_DATABASE_PATH, batch_size: int = 10_000) -> dict[str, int | str]:
    source_path = Path(source)
    output_path = Path(output)
    if not source_path.is_file():
        raise FileNotFoundError(f"Chemical Taxonomies source not found: {source_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    staging_path = output_path.with_suffix(output_path.suffix + ".building")
    if staging_path.exists():
        staging_path.unlink()

    connection = sqlite3.connect(staging_path)
    try:
        connection.executescript(
            """
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
            PRAGMA temp_store=MEMORY;
            CREATE TABLE taxonomy_records (
                record_id TEXT NOT NULL, chemical_name TEXT, casrn TEXT, dtxsid TEXT,
                inchikey TEXT NOT NULL, smiles TEXT, kingdom TEXT, superclass TEXT,
                class TEXT, subclass TEXT, level_5 TEXT, level_6 TEXT, level_7 TEXT,
                level_8 TEXT, level_9 TEXT, level_10 TEXT, data_version TEXT,
                date_extracted TEXT, reference TEXT, url TEXT
            );
            CREATE TABLE taxonomy_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            """
        )
        inserted = 0
        skipped = 0
        batch: list[tuple[str, ...]] = []
        with source_path.open("r", encoding="utf-8-sig", newline="") as source_file:
            reader = csv.DictReader(source_file, delimiter="\t")
            missing = [column for column in EXPECTED_COLUMNS if column not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"Source is missing expected columns: {', '.join(missing)}")
            for row in reader:
                inchikey = normalise_inchikey(_value(row, "InChiKey"))
                if not inchikey:
                    skipped += 1
                    continue
                batch.append((
                    _value(row, "Record_ID"), _value(row, "Chemical_Name"), _value(row, "CASRN"), _value(row, "DTXSID"),
                    inchikey, _value(row, "SMILES"), _value(row, "Kingdom"), _value(row, "Superclass"),
                    _value(row, "Class"), _value(row, "Subclass"), _value(row, "Level_5"), _value(row, "Level_6"),
                    _value(row, "Level_7"), _value(row, "Level_8"), _value(row, "Level_9"), _value(row, "Level_10"),
                    _value(row, "Data_Version"), _value(row, "Date_Extracted"), _value(row, "Reference"), _value(row, "URL"),
                ))
                if len(batch) >= batch_size:
                    connection.executemany("INSERT INTO taxonomy_records VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", batch)
                    inserted += len(batch)
                    batch.clear()
        if batch:
            connection.executemany("INSERT INTO taxonomy_records VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", batch)
            inserted += len(batch)
        connection.execute("CREATE INDEX idx_taxonomy_inchikey ON taxonomy_records(inchikey)")
        connection.execute("CREATE INDEX idx_taxonomy_casrn ON taxonomy_records(casrn)")
        connection.execute("CREATE INDEX idx_taxonomy_dtxsid ON taxonomy_records(dtxsid)")
        connection.executemany(
            "INSERT INTO taxonomy_metadata VALUES (?, ?)",
            [
                ("source_filename", source_path.name),
                ("source_bytes", str(source_path.stat().st_size)),
                ("record_count", str(inserted)),
                ("skipped_without_inchikey", str(skipped)),
                ("build_version", "1.0.0"),
            ],
        )
        connection.commit()
    finally:
        connection.close()

    if output_path.exists():
        output_path.unlink()
    os.replace(staging_path, output_path)
    return {"database": str(output_path), "records": inserted, "skipped": skipped}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the PRISM local Chemical Taxonomies SQLite database.")
    parser.add_argument("--source", required=True, help="Path to ChemicalTaxonomies_v1.1.txt")
    parser.add_argument("--output", default=str(DEFAULT_DATABASE_PATH), help="SQLite sidecar output path")
    arguments = parser.parse_args()
    print(build_database(arguments.source, arguments.output))
