"""Safe local helpers for using OPSIN as an unverified name-parser fallback."""

from __future__ import annotations

import re
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


OPSIN_WEB_SERVICE = "https://www.ebi.ac.uk/opsin/ws"

# These inputs identify a family, degradation product, mixture, or otherwise
# incomplete analytical annotation. A parser must never invent a discrete
# identity for them.
_INELIGIBLE_PATTERNS = (
    r"\brelated\s+compound\b",
    r"\bloss\s+of\b",
    r"\bdegradant\b",
    r"\badduct\b",
    r"\b(?:co)?polymer\b",
    r"\bpoly[a-z]*\b",
    r"\b(?:di|tri|tetra)mer\b",
    r"\bunknown\b",
    r"\bunidentified\b",
    r"\bnot\s+given\b",
    r"\bmixture\b",
    r"\bcompound\s*\d+\b",
    r"\b(?:cyanox|irganox|tinuvin)\b",
    r"\(\s*[nc]\s*=\s*\d+\s*\)",
    r"\bc\s*=\s*\d+\b",
)


def opsin_name_is_eligible(name: object) -> tuple[bool, str]:
    """Return whether a supplied name is safe to send to OPSIN."""
    text = str(name or "").strip()
    if not text:
        return False, "No chemical name supplied"
    if len(text) > 180:
        return False, "Name is too long for a safe deterministic parse"
    normalized = text.casefold()
    for pattern in _INELIGIBLE_PATTERNS:
        if re.search(pattern, normalized, re.IGNORECASE):
            return False, "Generic, polymeric, or incomplete analytical name"
    return True, "Eligible systematic-name candidate"


def opsin_identity_from_payload(name: str, payload: dict[str, Any]) -> dict[str, Any] | None:
    """Validate an OPSIN JSON payload and return the standard identity shape."""
    if str(payload.get("status", "")).upper() != "SUCCESS":
        return None
    smiles = str(payload.get("smiles") or "").strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchi = str(payload.get("inchi") or "").strip()
    inchikey = str(payload.get("stdInChIKey") or payload.get("inchiKey") or "").strip()
    if not inchi:
        try:
            inchi = Chem.MolToInchi(mol) or ""
        except Exception:
            inchi = ""
    if not inchikey:
        try:
            inchikey = Chem.MolToInchiKey(mol) or ""
        except Exception:
            inchikey = ""
    return {
        "status": "OPSIN parsed — identity review recommended",
        "source": "OPSIN systematic-name parser (unverified identity)",
        "name": name,
        "iupac_name": name,
        "cid": "",
        "synonyms": [],
        "cas": "",
        "smiles": canonical_smiles,
        "inchi": inchi,
        "inchikey": inchikey,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": float(Descriptors.MolWt(mol)),
    }
