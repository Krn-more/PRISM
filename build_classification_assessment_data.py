"""Build a JSON payload for the assessed classification workbook.

The source workbook remains read-only.  This materializes every source record
and its local deterministic classification so the workbook renderer can add
auditable result sheets without using Excel as a calculation engine.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
from rdkit import Chem

from classification_corpus import (
    _clean_text,
    _first_valid_structure,
    _resolution_route,
    _structure_status,
)
from structure_classification import classify_structure
from structure_classification import classification_export_fields
from structure_classification import classification_secondary_functional_groups


SOURCE_SHEET = "Compounds"


def _json_value(value: object) -> object:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    return value


def build_rows(source: str | Path) -> dict[str, object]:
    frame = pd.read_excel(source, sheet_name=SOURCE_SHEET)
    result_rows: list[dict[str, object]] = []
    curation_rows: list[dict[str, object]] = []

    for index, row in frame.iterrows():
        source_values = {column: _json_value(row[column]) for column in frame.columns}
        result = classify_structure(row.get("Standardized SMILES"), row.get("Original SMILES"))
        structure_status, structure_reason = _structure_status(row, result)
        hierarchy = result.get("taxonomy_hierarchy", {})
        topology = result.get("topology_profile", {})
        ontology = result.get("ontology_compatibility", {})
        features = result.get("feature_profile", {})
        final_details = result.get("final_classification_details", {})
        export_fields = classification_export_fields(result)
        secondary_groups, secondary_basis = classification_secondary_functional_groups(result)
        supplied_inchikey = _clean_text(row.get("InChIKey"))
        mol, structure_source, selected_smiles = _first_valid_structure(
            row.get("Standardized SMILES"), row.get("Original SMILES")
        )
        derived_inchikey = ""
        derived_status = "Not needed" if supplied_inchikey else "Unavailable"
        derivation_method = ""
        if not supplied_inchikey and mol is not None:
            try:
                derived_inchikey = Chem.MolToInchiKey(mol) or ""
                derived_status = "Derived locally" if derived_inchikey else "Unavailable"
                if derived_inchikey:
                    derivation_method = f"RDKit MolToInchiKey from {structure_source}"
            except Exception:
                derived_status = "Unavailable"
        output = {
            "Source Row": int(index) + 2,
            **source_values,
            "Structure Status": structure_status,
            "Structure Status Reason": structure_reason,
            "Classification Input Source": export_fields.get("Classification Input Source", ""),
            "Classification Structure SMILES": export_fields.get("Classification Structure SMILES", ""),
            "Classification Standardization Version": export_fields.get("Classification Standardization Version", ""),
            "Classification Stereochemistry Status": export_fields.get("Classification Stereochemistry Status", ""),
            "Classification Status": export_fields.get("Classification Status", ""),
            "Classification Scope": export_fields.get("Classification Scope", ""),
            "Classification Review Recommended": export_fields.get("Classification Review Recommended", False),
            "Classification Review Recommendation": export_fields.get("Classification Review Recommendation", ""),
            "Classification Suggested Action": export_fields.get("Classification Suggested Action", ""),
            "Corrected Chemical Class": export_fields.get("Corrected Chemical Class", ""),
            "Classification Rule ID": export_fields.get("Classification Rule ID", ""),
            "Classification Rule Version": export_fields.get("Classification Rule Version", ""),
            "Detected Feature Profile": export_fields.get("Detected Feature Profile", ""),
            "Taxonomy Path": export_fields.get("Taxonomy Path", ""),
            "Taxonomy Path Steps": export_fields.get("Taxonomy Path Steps", ""),
            "Matched Categories": export_fields.get("Matched Categories", ""),
            "Direct Parent": export_fields.get("Direct Parent", ""),
            "Taxonomy Dictionary Version": export_fields.get("Taxonomy Dictionary Version", ""),
            "Parent Class": hierarchy.get("Parent Class", ""),
            "Primary Functional Group": hierarchy.get("Functional Group", ""),
            "Secondary Functional Groups": secondary_groups,
            "Secondary Functional Groups Basis": secondary_basis,
            "Subclass": hierarchy.get("Subclass", ""),
            "Structural Type": hierarchy.get("Structural Type", ""),
            "Topology": final_details.get("Topology Profile", ""),
            "Ring Profile": topology.get("Ring Profile", ""),
            "Molecular Weight": features.get("molecular_weight", ""),
            "H-Bond Donors": features.get("hbd", ""),
            "H-Bond Acceptors": features.get("hba", ""),
            "Rotatable Bonds": features.get("rotatable_bonds", ""),
            "Manual Review Flag": bool(result.get("manual_review_flag", False)),
            "Manual Review Reason": result.get("manual_review_reason", ""),
            "Supplied InChIKey": supplied_inchikey,
            "Derived InChIKey Candidate": derived_inchikey,
            "Derived InChIKey Status": derived_status,
            "InChIKey Derivation Method": derivation_method,
            "Ontology Mapping Status": ontology.get("mapping_status", ""),
            "Ontology Mapping Version": ontology.get("ontology_mapping_version", ""),
            "Ontology Kingdom": ontology.get("kingdom", ""),
            "Ontology Superclass": ontology.get("superclass", ""),
            "Ontology Class": ontology.get("class", ""),
            "Ontology Subclass": ontology.get("subclass", ""),
            "Ontology Direct Parent": ontology.get("direct_parent", ""),
            "Ontology Mapping Evidence": ontology.get("mapping_evidence", ""),
            "Final Classification Record": result.get("final_classification_record", ""),
        }
        result_rows.append(output)
        if structure_status != "Structured":
            evidence_tier, resolution_path, eligibility = _resolution_route(
                structure_status, _clean_text(row.get("CASRN"))
            )
            output["Curation Evidence Tier"] = evidence_tier
            output["Curation Resolution Path"] = resolution_path
            output["Automation Eligibility"] = eligibility
            output["Suggested Next Action"] = (
                "Review source evidence and obtain an approved structure; do not infer from name text"
            )
            curation_rows.append(output)

    return {
        "source_columns": list(frame.columns),
        "result_columns": list(result_rows[0]) if result_rows else [],
        "rows": result_rows,
        "curation_rows": curation_rows,
    }


if __name__ == "__main__":
    source_path = Path("List of Inchikey and SMILES.xlsx")
    target_path = Path("fixtures/classification_workbook_rows.json")
    target_path.write_text(json.dumps(build_rows(source_path), ensure_ascii=False), encoding="utf-8")
    print(target_path)
