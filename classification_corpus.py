"""Reproducible local baseline and curation outputs for a classification corpus.

This module never modifies the source workbook.  It records the current local
classification result, identifies records that require identity curation, and
derives candidate InChIKeys only from already-valid local SMILES.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from structure_classification import classify_structure


SOURCE_SHEET = "Compounds"
REQUIRED_COLUMNS = {
    "Compound Name", "CASRN", "Source", "Original SMILES",
    "Standardized SMILES", "InChIKey",
}
AMBIGUOUS_NAME_PATTERN = re.compile(r"\b(related compound|loss of|c=\d+)\b", re.IGNORECASE)


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def _first_valid_structure(standardized_smiles: object, original_smiles: object) -> tuple[Chem.Mol | None, str, str]:
    for source, value in (("Standardized SMILES", standardized_smiles), ("Original SMILES", original_smiles)):
        text = _clean_text(value)
        if not text:
            continue
        mol = Chem.MolFromSmiles(text)
        if mol is not None:
            return mol, source, text
    return None, "", ""


def _structure_status(row: pd.Series, result: dict[str, Any]) -> tuple[str, str]:
    standardized = _clean_text(row.get("Standardized SMILES"))
    original = _clean_text(row.get("Original SMILES"))
    name = _clean_text(row.get("Compound Name"))
    if not standardized and not original:
        reason = "No standardized or original SMILES supplied"
    elif result.get("manual_review_flag"):
        reason = "Supplied structure could not be parsed"
    else:
        return "Structured", "Valid local structure available"
    if AMBIGUOUS_NAME_PATTERN.search(name):
        return "Ambiguous identity", f"{reason}; name contains an ambiguity annotation"
    return "Needs structure" if not standardized and not original else "Invalid structure", reason


def _resolution_route(curation_status: str, supplied_cas: str) -> tuple[str, str, str]:
    """Choose a reviewer route without treating a name as a chemical identity."""
    if curation_status == "Ambiguous identity":
        return (
            "Insufficient identity evidence",
            "Review the source table and analytical evidence; do not create a structure from the annotation text",
            "Ineligible until identity is approved",
        )
    if supplied_cas and supplied_cas.lower() not in {"not given", "nan", "none"}:
        return (
            "CASRN supplied — verify scope",
            "Resolve through an approved reference source by CASRN, then verify structure against the reported analytical context",
            "Ineligible until structure is reviewer-approved",
        )
    return (
        "Name-only evidence",
        "Obtain an authoritative structure from source documentation or an approved reference source; retain name-only status if no unique identity exists",
        "Ineligible until structure is reviewer-approved",
    )


def _source_sha256(source_path: Path) -> str:
    digest = hashlib.sha256()
    with source_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_corpus_outputs(source_path: str | Path, output_dir: str | Path) -> dict[str, Any]:
    """Generate deterministic CSV/JSON outputs from the supplied corpus workbook."""
    source = Path(source_path)
    target = Path(output_dir)
    frame = pd.read_excel(source, sheet_name=SOURCE_SHEET)
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"Corpus is missing required columns: {', '.join(sorted(missing))}")

    baseline_rows: list[dict[str, Any]] = []
    curation_rows: list[dict[str, Any]] = []
    derived_rows: list[dict[str, Any]] = []
    rule_review_rows: list[dict[str, Any]] = []

    for index, row in frame.iterrows():
        source_row = int(index) + 2
        result = classify_structure(row.get("Standardized SMILES"), row.get("Original SMILES"))
        structure_status, structure_reason = _structure_status(row, result)
        supplied_inchikey = _clean_text(row.get("InChIKey"))
        mol, selected_structure_source, selected_smiles = _first_valid_structure(
            row.get("Standardized SMILES"), row.get("Original SMILES")
        )
        derived_inchikey = ""
        derived_status = "Not needed" if supplied_inchikey else "Unavailable"
        if not supplied_inchikey and mol is not None:
            try:
                derived_inchikey = Chem.MolToInchiKey(mol)
                derived_status = "Derived locally" if derived_inchikey else "Unavailable"
            except Exception:
                derived_status = "Unavailable"
        ontology = result.get("ontology_compatibility", {})
        detected_features = result.get("detected_feature_profile", [])
        baseline_rows.append({
            "Source Row": source_row,
            "Compound Name": _clean_text(row.get("Compound Name")),
            "CASRN": _clean_text(row.get("CASRN")),
            "Supplied InChIKey": supplied_inchikey,
            "Derived InChIKey": derived_inchikey,
            "Derived InChIKey Status": derived_status,
            "Structure Status": structure_status,
            "Structure Status Reason": structure_reason,
            "Classification Input Source": result.get("classification_input_source", ""),
            "Corrected Chemical Class": result.get("corrected_chemical_class", ""),
            "Classification Rule ID": result.get("classification_rule_id", ""),
            "Classification Rule Version": result.get("classification_rule_version", ""),
            "Manual Review Flag": bool(result.get("manual_review_flag", False)),
            "Manual Review Reason": result.get("manual_review_reason", ""),
            "Taxonomy Path": result.get("taxonomy_path", ""),
            "Ontology Mapping Status": ontology.get("mapping_status", ""),
            "Ontology Mapping Version": ontology.get("ontology_mapping_version", ""),
        })
        if not supplied_inchikey and derived_inchikey:
            derived_rows.append({
                "Source Row": source_row,
                "Compound Name": _clean_text(row.get("Compound Name")),
                "CASRN": _clean_text(row.get("CASRN")),
                "Structure Source": selected_structure_source,
                "SMILES Used": selected_smiles,
                "Derived InChIKey": derived_inchikey,
                "Derivation Method": "RDKit MolToInchiKey from locally validated SMILES",
                "Review Status": "Candidate — retain supplied field unchanged until approved",
            })
        if structure_status != "Structured":
            evidence_tier, resolution_path, eligibility = _resolution_route(structure_status, _clean_text(row.get("CASRN")))
            curation_rows.append({
                "Source Row": source_row,
                "Compound Name": _clean_text(row.get("Compound Name")),
                "CASRN": _clean_text(row.get("CASRN")),
                "Source": _clean_text(row.get("Source")),
                "Standardized SMILES": _clean_text(row.get("Standardized SMILES")),
                "Original SMILES": _clean_text(row.get("Original SMILES")),
                "Supplied InChIKey": supplied_inchikey,
                "Curation Status": structure_status,
                "Reason": structure_reason,
                "Evidence Tier": evidence_tier,
                "Resolution Path": resolution_path,
                "Automation Eligibility": eligibility,
                "Suggested Next Action": "Review source evidence and obtain an approved structure; do not infer from name text",
            })
        generic_fallback = str(result.get("classification_rule_id", "")).endswith("-UNK-001")
        if structure_status == "Structured" and generic_fallback:
            rule_review_rows.append({
                "Source Row": source_row,
                "Compound Name": _clean_text(row.get("Compound Name")),
                "CASRN": _clean_text(row.get("CASRN")),
                "Standardized SMILES": _clean_text(row.get("Standardized SMILES")),
                "Corrected Chemical Class": result.get("corrected_chemical_class", ""),
                "Classification Rule ID": result.get("classification_rule_id", ""),
                "Detected Feature Profile": "; ".join(detected_features),
                "Review Reason": "Broad elemental fallback selected despite a valid structure",
                "Recommended Review": "Approve a narrower structural rule only with positive, negative, and multifunctional controls",
            })

    target.mkdir(parents=True, exist_ok=True)
    baseline = pd.DataFrame(baseline_rows)
    curation = pd.DataFrame(curation_rows)
    derived = pd.DataFrame(derived_rows)
    rule_review = pd.DataFrame(rule_review_rows)
    baseline_path = target / "classification_corpus_baseline.csv"
    curation_path = target / "classification_identity_curation_queue.csv"
    derived_path = target / "classification_derived_inchikey_candidates.csv"
    rule_review_path = target / "classification_rule_review_queue.csv"
    baseline.to_csv(baseline_path, index=False, encoding="utf-8")
    curation.to_csv(curation_path, index=False, encoding="utf-8")
    derived.to_csv(derived_path, index=False, encoding="utf-8")
    rule_review.to_csv(rule_review_path, index=False, encoding="utf-8")

    summary = {
        "source_workbook": source.name,
        "source_sheet": SOURCE_SHEET,
        "source_sha256": _source_sha256(source),
        "row_count": len(frame),
        "classified_rows": int((~baseline["Manual Review Flag"]).sum()),
        "manual_review_rows": int(baseline["Manual Review Flag"].sum()),
        "missing_supplied_inchikey_rows": int((baseline["Supplied InChIKey"] == "").sum()),
        "derived_inchikey_candidate_rows": len(derived),
        "structure_status_counts": dict(Counter(baseline["Structure Status"])),
        "curation_evidence_tier_counts": dict(Counter(curation["Evidence Tier"])) if not curation.empty else {},
        "rule_review_queue_rows": len(rule_review),
        "class_counts": dict(Counter(baseline["Corrected Chemical Class"])),
        "rule_counts": dict(Counter(baseline["Classification Rule ID"])),
        "ontology_mapping_status_counts": dict(Counter(baseline["Ontology Mapping Status"])),
        "outputs": {
            "baseline": baseline_path.name,
            "curation_queue": curation_path.name,
            "derived_inchikey_candidates": derived_path.name,
            "rule_review_queue": rule_review_path.name,
        },
    }
    summary_path = target / "classification_corpus_baseline.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build local classification-corpus baseline outputs.")
    parser.add_argument("--source", default="List of Inchikey and SMILES.xlsx")
    parser.add_argument("--output-dir", default="fixtures")
    args = parser.parse_args()
    print(json.dumps(build_corpus_outputs(args.source, args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
