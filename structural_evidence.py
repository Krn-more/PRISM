"""Deterministic Stage 1 structural-evidence calculations.

All values in this module are calculated from standardized SMILES using RDKit.
They are supporting grouping evidence only and do not alter clustering decisions.
"""

from __future__ import annotations

import os
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, Lipinski, rdMolDescriptors, Descriptors3D, AllChem


STRUCTURAL_EVIDENCE_VERSION = "1.2.0"

# RDKit's legacy BRENK catalogue labels organosilicon structures as
# ``heavy_metal``. Silicon is a metalloid, so that screen is not valid for a
# siloxane or other organosilicon compound unless a genuine heavy-metal atom
# is also present. Keep genuine metal alerts while suppressing this known
# catalogue false positive.
HEAVY_METAL_ATOMIC_NUMBERS = {
    22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
}

# These are structural indicators, not pKa predictions or measured ionisation.
IONISATION_PATTERNS = {
    "Likely acidic": [
        "[CX3](=O)[OX2H1,-]",  # carboxylic acid/carboxylate
        "[SX4](=O)(=O)[OX2H1,-]",  # sulfonic acid/sulfonate
        "[PX4](=O)([OX2H1,-])[OX2H1,-]",  # phosphate/phosphonate motif
    ],
    "Likely basic": [
        "[NX3;!$(N[C,S]=O);!$(N=*)]",  # non-amide amines
        "[nH0;+0]",  # neutral aromatic nitrogens, e.g. pyridine-like
        "[NX2]=[CX2]",  # amidine-like motif
    ],
    "Permanently charged": ["[N+;!$(N[O-])](*) (*) (*)", "[P+](*) (*) (*) (*)", "[S+](*) (*) (*)"],
}

TOXICOPHORE_PATTERNS = {
    "Aldehyde": "[CX3H1](=O)[#6]",
    "Alkyl halide": "[CX4][Cl,Br,I]",
    "Aromatic amine": "[NX3;H0,H1,H2][a]",
    "Azide": "[$([N-]=[N+]=N),$([N]=[N+]=[N-])]",
    "Epoxide": "[OX2r3]1[CX4r3][CX4r3]1",
    "Hydrazine": "[NX3,NX2][NX3,NX2]",
    "Michael acceptor": "[CX3]=[CX3][CX3](=O)[#6,O,N,S]",
    "Nitrosamine": "[NX3][NX2]=O",
    "Peroxide": "[OX2][OX2]",
    "Sulfonyl halide": "[SX4](=O)(=O)[Cl,Br]",
}


def _compiled_patterns(patterns: dict[str, Any]) -> dict[str, list[Chem.Mol]]:
    compiled: dict[str, list[Chem.Mol]] = {}
    for name, value in patterns.items():
        values = value if isinstance(value, list) else [value]
        compiled[name] = [pattern for smarts in values if (pattern := Chem.MolFromSmarts(smarts))]
    return compiled


IONISATION_SMARTS = _compiled_patterns(IONISATION_PATTERNS)
TOXICOPHORE_SMARTS = _compiled_patterns(TOXICOPHORE_PATTERNS)


def _matches_any(mol: Chem.Mol, patterns: list[Chem.Mol]) -> bool:
    return any(mol.HasSubstructMatch(pattern) for pattern in patterns)


def _all_alerts(mol: Chem.Mol) -> list[str]:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)
    return sorted({entry.GetDescription() for entry in catalog.GetMatches(mol) if _retain_alert(entry.GetDescription(), mol)})


def _first_alert(mol: Chem.Mol) -> str:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)
    for match in catalog.GetMatches(mol):
        if _retain_alert(match.GetDescription(), mol):
            return match.GetDescription()
    return "None"


def _retain_alert(description: str, mol: Chem.Mol) -> bool:
    """Filter only the documented organosilicon heavy-metal false positive."""
    if str(description).strip().casefold() != "heavy_metal":
        return True
    return any(atom.GetAtomicNum() in HEAVY_METAL_ATOMIC_NUMBERS for atom in mol.GetAtoms())


def _ionisation_indicator(mol: Chem.Mol) -> str:
    if _matches_any(mol, IONISATION_SMARTS["Permanently charged"]):
        return "Permanently charged"
    acidic = _matches_any(mol, IONISATION_SMARTS["Likely acidic"])
    basic = _matches_any(mol, IONISATION_SMARTS["Likely basic"])
    if acidic and basic:
        return "Amphoteric"
    if acidic:
        return "Likely acidic"
    if basic:
        return "Likely basic"
    return "Neutral"


def calculate_structural_evidence(smiles: object) -> dict[str, Any]:
    """Return Stage 1 evidence for one standardized SMILES string."""
    if smiles is None or pd.isna(smiles) or not str(smiles).strip():
        return _empty_evidence("Missing structure")
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return _empty_evidence("Invalid structure")

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    aromatic_ring_count = sum(all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring) for ring in atom_rings)
    toxicophores = sorted(name for name, patterns in TOXICOPHORE_SMARTS.items() if _matches_any(mol, patterns))
    alerts = _all_alerts(mol)

    # Match the workbook workflow by exposing the same descriptor family used
    # for Step 5 in the Excel report.
    mw = round(Descriptors.ExactMolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    tpsa = round(Descriptors.TPSA(mol), 2)
    min_estate = round(Descriptors.MinEStateIndex(mol), 2)
    max_estate = round(Descriptors.MaxEStateIndex(mol), 2)
    enable_3d = os.environ.get("PRISM_ENABLE_3D_EVIDENCE", "false").strip().casefold() in {"1", "true", "yes", "on"}
    asphericity = pmi1 = pmi2 = pmi3 = radius_gyration = None
    evidence_3d_status = "Skipped in fast deterministic mode"
    if enable_3d:
        try:
            mol3d = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol3d, maxAttempts=10, randomSeed=42) != -1:
                asphericity = round(Descriptors3D.Asphericity(mol3d), 3)
                pmi1 = round(Descriptors3D.PMI1(mol3d), 3)
                pmi2 = round(Descriptors3D.PMI2(mol3d), 3)
                pmi3 = round(Descriptors3D.PMI3(mol3d), 3)
                radius_gyration = round(Descriptors3D.RadiusOfGyration(mol3d), 3)
                evidence_3d_status = "Calculated"
            else:
                asphericity = pmi1 = pmi2 = pmi3 = radius_gyration = 0.0
                evidence_3d_status = "Embedding unavailable"
        except Exception:
            asphericity = pmi1 = pmi2 = pmi3 = radius_gyration = 0.0
            evidence_3d_status = "Embedding unavailable"

    return {
        "Structural_Evidence_Status": "Calculated",
        "Structural_Evidence_Version": STRUCTURAL_EVIDENCE_VERSION,
        "3D_Evidence_Status": evidence_3d_status,
        "MW": mw,
        "LogP": logp,
        "TPSA": tpsa,
        "Min_EState": min_estate,
        "Max_EState": max_estate,
        "3D_Asphericity": asphericity,
        "3D_PMI1": pmi1,
        "3D_PMI2": pmi2,
        "3D_PMI3": pmi3,
        "3D_RadiusOfGyration": radius_gyration,
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Rotatable_Bonds": Lipinski.NumRotatableBonds(mol),
        "Formal_Charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "Ring_Count": len(atom_rings),
        "Aromatic_Ring_Count": aromatic_ring_count,
        "Fraction_CSP3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "Heavy_Atom_Count": mol.GetNumHeavyAtoms(),
        "Molecular_Refractivity": round(Crippen.MolMR(mol), 3),
        "Ionisation_Indicator": _ionisation_indicator(mol),
        "Toxicophore_Profile": "; ".join(toxicophores) if toxicophores else "None",
        "All_Structural_Alerts": "; ".join(alerts) if alerts else "None",
        "Alerts": _first_alert(mol),
    }


def _empty_evidence(status: str) -> dict[str, Any]:
    values = {
        "Structural_Evidence_Status": status,
        "Structural_Evidence_Version": STRUCTURAL_EVIDENCE_VERSION,
        "3D_Evidence_Status": "Not assessable",
        "MW": None, "LogP": None, "TPSA": None, "Min_EState": None, "Max_EState": None,
        "3D_Asphericity": None, "3D_PMI1": None, "3D_PMI2": None, "3D_PMI3": None,
        "3D_RadiusOfGyration": None,
        "HBD": None, "HBA": None, "Rotatable_Bonds": None, "Formal_Charge": None,
        "Ring_Count": None, "Aromatic_Ring_Count": None, "Fraction_CSP3": None,
        "Heavy_Atom_Count": None, "Molecular_Refractivity": None,
        "Ionisation_Indicator": "Unknown", "Toxicophore_Profile": "Not assessable",
        "All_Structural_Alerts": "Not assessable", "Alerts": "Not assessable",
    }
    return values


def add_structural_evidence(df: pd.DataFrame, smiles_col: str = "Standardized SMILES") -> pd.DataFrame:
    """Refresh Stage 1 evidence columns without creating duplicate headers."""
    result = df.copy()
    if smiles_col not in result.columns:
        for column, value in _empty_evidence("SMILES column missing").items():
            result[column] = value
        return result
    evidence = pd.DataFrame([calculate_structural_evidence(smiles) for smiles in result[smiles_col]])
    evidence.index = result.index
    # A Phase A workbook already carries these calculated fields.  Assigning
    # by column refreshes stale values while preserving one unambiguous header
    # for the B1a consensus copy-back.
    for column in evidence.columns:
        result[column] = evidence[column]
    return result
