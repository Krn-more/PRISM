from __future__ import annotations

from typing import Any, Mapping


DESCRIPTOR_EXPORT_FIELDS = [
    "Structural_Evidence_Status",
    "Structural_Evidence_Version",
    "MW",
    "LogP",
    "TPSA",
    "Min_EState",
    "Max_EState",
    "3D_Asphericity",
    "3D_PMI1",
    "3D_PMI2",
    "3D_PMI3",
    "3D_RadiusOfGyration",
    "HBD",
    "HBA",
    "Rotatable_Bonds",
    "Formal_Charge",
    "Ring_Count",
    "Aromatic_Ring_Count",
    "Fraction_CSP3",
    "Heavy_Atom_Count",
    "Molecular_Refractivity",
    "Ionisation_Indicator",
    "Toxicophore_Profile",
    "All_Structural_Alerts",
    "Alerts",
]


def build_descriptor_export_fields(
    structural_evidence: Mapping[str, Any] | None,
    identity: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the Step 5/6 workbook mirror fields from a single evidence source."""
    evidence = dict(structural_evidence or {})
    identity = dict(identity or {})
    return {
        "Structural Evidence Status": evidence.get("Structural_Evidence_Status", ""),
        "Structural Evidence Version": evidence.get("Structural_Evidence_Version", ""),
        "MW": evidence.get("MW", identity.get("molecular_weight")),
        "LogP": evidence.get("LogP"),
        "TPSA": evidence.get("TPSA"),
        "Min_EState": evidence.get("Min_EState"),
        "Max_EState": evidence.get("Max_EState"),
        "3D_Asphericity": evidence.get("3D_Asphericity"),
        "3D_PMI1": evidence.get("3D_PMI1"),
        "3D_PMI2": evidence.get("3D_PMI2"),
        "3D_PMI3": evidence.get("3D_PMI3"),
        "3D_RadiusOfGyration": evidence.get("3D_RadiusOfGyration"),
        "HBD": evidence.get("HBD"),
        "HBA": evidence.get("HBA"),
        "Rotatable Bonds": evidence.get("Rotatable_Bonds"),
        "Formal Charge": evidence.get("Formal_Charge"),
        "Ring Count": evidence.get("Ring_Count"),
        "Aromatic Ring Count": evidence.get("Aromatic_Ring_Count"),
        "Fraction CSP3": evidence.get("Fraction_CSP3"),
        "Heavy Atom Count": evidence.get("Heavy_Atom_Count"),
        "Molecular Refractivity": evidence.get("Molecular_Refractivity"),
        "Ionisation Indicator": evidence.get("Ionisation_Indicator"),
        "Toxicophore Profile": evidence.get("Toxicophore_Profile"),
        "All Structural Alerts": evidence.get("All_Structural_Alerts"),
        "Alerts": evidence.get("Alerts"),
    }
