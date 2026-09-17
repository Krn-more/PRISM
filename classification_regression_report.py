"""Compare the approved classification baseline with the active local rules.

This tool is deliberately offline.  It treats the assessed workbook as a
baseline, re-runs the current deterministic classifier for each source row,
and writes only changed governance fields for review before a rule-set release.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from build_classification_assessment_data import build_rows
from structure_classification import classification_export_fields, classify_structure


TRACKED_FIELDS = (
    "Corrected Chemical Class",
    "Classification Rule ID",
    "Classification Rule Version",
    "Classification Status",
    "Classification Scope",
    "Taxonomy Path",
    "Direct Parent",
    "Manual Review Flag",
    "Manual Review Reason",
    "Classification Suggested Action",
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def current_classification_rows(source: str | Path) -> pd.DataFrame:
    """Materialize current governance fields using the shared classifier."""
    source_rows = build_rows(source)["rows"]
    records: list[dict[str, object]] = []
    for source_row in source_rows:
        result = classify_structure(
            source_row.get("Standardized SMILES", ""),
            source_row.get("Original SMILES", ""),
        )
        fields = classification_export_fields(result)
        records.append({
            "Source Row": source_row["Source Row"],
            "Compound Name": source_row.get("Compound Name", ""),
            "CASRN": source_row.get("CASRN", ""),
            **fields,
        })
    return pd.DataFrame(records)


def build_regression_report(source: str | Path, baseline: str | Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """Return changed classifier fields, keyed by the stable source-row ID."""
    current = current_classification_rows(source)
    baseline_df = pd.read_excel(baseline, sheet_name="PRISM Classification")
    baseline_df = baseline_df.drop_duplicates(subset=["Source Row"], keep="last")
    merged = current.merge(baseline_df, on="Source Row", how="left", suffixes=(" Current", " Baseline"))

    changes: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        for field in TRACKED_FIELDS:
            current_value = _text(row.get(f"{field} Current", row.get(field, "")))
            baseline_value = _text(row.get(f"{field} Baseline", ""))
            # Newly added governance fields are recorded as additions, not
            # silently treated as baseline mismatches.
            if current_value != baseline_value:
                changes.append({
                    "Source Row": row["Source Row"],
                    "Compound Name": row.get("Compound Name Current", ""),
                    "CASRN": row.get("CASRN Current", ""),
                    "Field": field,
                    "Baseline Value": baseline_value,
                    "Current Value": current_value,
                    "Change Type": "Added field" if not baseline_value else "Changed",
                })

    report = pd.DataFrame(changes, columns=[
        "Source Row", "Compound Name", "CASRN", "Field", "Baseline Value", "Current Value", "Change Type",
    ])
    summary = {
        "source_rows": len(current),
        "changed_rows": int(report["Source Row"].nunique()) if not report.empty else 0,
        "changed_fields": len(report),
    }
    return report, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a PRISM classification regression report.")
    parser.add_argument("--source", default="List of Inchikey and SMILES.xlsx")
    parser.add_argument("--baseline", default="outputs/classification_corpus_assessed.xlsx")
    parser.add_argument("--output", default="outputs/classification_regression_report.csv")
    args = parser.parse_args()

    report, summary = build_regression_report(args.source, args.baseline)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(output, index=False, encoding="utf-8-sig")
    print({**summary, "output": str(output)})


if __name__ == "__main__":
    main()
