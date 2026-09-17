"""Deterministic local structure classification.

This module is intentionally offline, rule-driven, and audit-friendly. It
returns one primary corrected class plus a concise final record, while
preserving all detected supporting features for traceability.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from local_taxonomy_db import ontology_evidence


DEFAULT_RULES_PATH = Path(__file__).resolve().parent / "config" / "structure_classification_rules.json"
DEFAULT_ONTOLOGY_MAPPING_PATH = Path(__file__).resolve().parent / "config" / "ontology_compatibility_rules.json"
STANDARDIZATION_VERSION = "RDKit canonical-isomeric-v1"

# These rules prove that a valid structure belongs to a broad elemental family,
# but they do not claim a fully approved specific class.  Keep them separate
# from an identity/structure failure so downstream clustering can use the
# structure while reviewers can still prioritize taxonomy expansion.
BROAD_PROVISIONAL_RULE_IDS = {
    "STC-P-UNK-001",
    "STC-S-UNK-001",
    "STC-SI-UNK-001",
    "STC-N-UNK-001",
    "STC-O-UNK-001",
}


def _config_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "structure_classification_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _config_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                return json.load(handle)
    raise FileNotFoundError("structure_classification_rules.json not found")


@lru_cache(maxsize=2)
def load_ontology_compatibility_mappings(mapping_path: str | Path | None = None) -> dict[str, Any]:
    """Load approved local ontology mappings; never call a remote service."""
    candidate = Path(mapping_path) if mapping_path else DEFAULT_ONTOLOGY_MAPPING_PATH
    if not candidate.exists():
        return {"ontology_mapping_version": "unavailable", "mappings": {}}
    with open(candidate, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload.get("mappings"), dict):
        raise ValueError("ontology_compatibility_rules.json must contain a mappings object")
    return payload


def _ontology_compatibility(rule_id: str, corrected_class: str) -> dict[str, str]:
    """Return a supplementary, approved ontology mapping for one local rule."""
    payload = load_ontology_compatibility_mappings()
    mapping = payload.get("mappings", {}).get(rule_id)
    version = str(payload.get("ontology_mapping_version", "unavailable"))
    if not mapping:
        return {
            "mapping_status": "Not mapped",
            "ontology_mapping_version": version,
            "local_rule_id": rule_id,
            "local_class": corrected_class,
            "mapping_evidence": "No approved local ontology mapping for this structural rule",
        }
    if mapping.get("local_class") != corrected_class:
        raise ValueError(f"Ontology mapping class mismatch for {rule_id}")
    return {
        "mapping_status": str(mapping.get("mapping_status", "Mapped")),
        "ontology_mapping_version": version,
        "local_rule_id": rule_id,
        "local_class": corrected_class,
        "kingdom": str(mapping.get("kingdom", "")),
        "superclass": str(mapping.get("superclass", "")),
        "class": str(mapping.get("class", "")),
        "subclass": str(mapping.get("subclass", "")),
        "direct_parent": str(mapping.get("direct_parent", "")),
        "mapping_evidence": str(mapping.get("mapping_evidence", "")),
    }


def _mol_from_smiles(smiles: object) -> Chem.Mol | None:
    if smiles is None or pd.isna(smiles):
        return None
    text = str(smiles).strip()
    if not text:
        return None
    return Chem.MolFromSmiles(text)


def _ring_summary(mol: Chem.Mol) -> tuple[int, int, list[tuple[int, ...]]]:
    ring_info = mol.GetRingInfo()
    atom_rings = list(ring_info.AtomRings())
    aromatic_ring_count = sum(
        all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring)
        for ring in atom_rings
    )
    return len(atom_rings), aromatic_ring_count, atom_rings


def _topology_profile_for_mol(mol: Chem.Mol, feature_profile: dict[str, Any]) -> dict[str, Any]:
    ring_count, aromatic_ring_count, atom_rings = _ring_summary(mol)
    ring_sizes = [len(ring) for ring in atom_rings]
    ring_atoms = {atom_idx for ring in atom_rings for atom_idx in ring}
    hetero_ring_atom_count = sum(
        1
        for atom_idx in ring_atoms
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() not in {1, 6}
    )

    shared_atom_sizes: list[int] = []
    spiro_pair_count = 0
    fused_pair_count = 0
    for i, ring_a in enumerate(atom_rings):
        set_a = set(ring_a)
        for ring_b in atom_rings[i + 1:]:
            shared = set_a.intersection(ring_b)
            if not shared:
                continue
            shared_size = len(shared)
            shared_atom_sizes.append(shared_size)
            if shared_size == 1:
                spiro_pair_count += 1
            elif shared_size >= 2:
                fused_pair_count += 1

    max_ring_size = max(ring_sizes, default=0)
    topology_class = "Acyclic"
    topology_modifiers: list[str] = []
    alcohol_substitution_labels: list[str] = []
    alcohol_context_labels: list[str] = []

    alcohol_matches = _matches(mol, "[OX2H][C;!$(C=O)]")
    for match in alcohol_matches:
        if len(match) < 2:
            continue
        oxygen_idx, carbon_idx = match[0], match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [
            neighbor for neighbor in carbon_atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != oxygen_idx
        ]
        if len(carbon_neighbors) >= 3:
            alcohol_substitution_labels.append("Tertiary Alcohol")
        elif len(carbon_neighbors) == 2:
            alcohol_substitution_labels.append("Secondary Alcohol")
        else:
            alcohol_substitution_labels.append("Primary Alcohol")

        if carbon_atom.IsInRing():
            alcohol_context_labels.append("Cyclic Alcohol")
        elif any(neighbor.IsInRing() for neighbor in carbon_neighbors):
            alcohol_context_labels.append("Ring-substituted Alcohol")
        else:
            alcohol_context_labels.append("Acyclic Alcohol")

    if ring_count:
        if max_ring_size >= 12:
            topology_class = "Macrocyclic"
        elif ring_count == 1:
            topology_class = "Monocyclic"
        elif fused_pair_count:
            topology_class = "Fused Polycyclic"
        elif spiro_pair_count:
            topology_class = "Spirocyclic"
        elif ring_count == 2:
            topology_class = "Bicyclic"
        else:
            topology_class = "Polycyclic"

        if aromatic_ring_count:
            topology_modifiers.append("Aromatic")
        if hetero_ring_atom_count:
            topology_modifiers.append("Heterocyclic")
        if fused_pair_count:
            topology_modifiers.append("Fused")
        if spiro_pair_count:
            topology_modifiers.append("Spiro")
        if max_ring_size >= 12:
            topology_modifiers.append("Macrocyclic")

    ring_profile_text = f"{ring_count} ring(s), {aromatic_ring_count} aromatic ring(s)"
    ring_system_text = "Open-chain" if ring_count == 0 else topology_class

    return {
        "Topology Class": topology_class,
        "Topology Modifiers": "; ".join(topology_modifiers) if topology_modifiers else "None",
        "Ring System": ring_system_text,
        "Ring Profile": ring_profile_text,
        "Ring Count": ring_count,
        "Aromatic Ring Count": aromatic_ring_count,
        "Ring Sizes": ring_sizes,
        "Hetero Ring Atom Count": hetero_ring_atom_count,
        "Shared Ring Atom Sizes": shared_atom_sizes,
        "Alcohol Substitution": (
            "Tertiary Alcohol"
            if "Tertiary Alcohol" in alcohol_substitution_labels
            else "Secondary Alcohol"
            if "Secondary Alcohol" in alcohol_substitution_labels
            else "Primary Alcohol"
            if "Primary Alcohol" in alcohol_substitution_labels
            else "Not applicable"
        ),
        "Alcohol Context": (
            "Cyclic Alcohol"
            if "Cyclic Alcohol" in alcohol_context_labels
            else "Ring-substituted Alcohol"
            if "Ring-substituted Alcohol" in alcohol_context_labels
            else "Acyclic Alcohol"
            if "Acyclic Alcohol" in alcohol_context_labels
            else "Not applicable"
        ),
        "Alcohol Contexts": sorted(set(alcohol_context_labels)),
    }


def _compiled(smarts: str) -> Chem.Mol | None:
    return Chem.MolFromSmarts(smarts)


def _has_match(mol: Chem.Mol, smarts: str) -> bool:
    pattern = _compiled(smarts)
    return bool(pattern and mol.HasSubstructMatch(pattern))


def _matches(mol: Chem.Mol, smarts: str) -> list[tuple[int, ...]]:
    pattern = _compiled(smarts)
    return list(mol.GetSubstructMatches(pattern)) if pattern else []


def _feature_profile(mol: Chem.Mol, rules: dict[str, Any]) -> dict[str, Any]:
    thresholds = rules.get("thresholds", {})
    ring_count, aromatic_ring_count, atom_rings = _ring_summary(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
    oxygen_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8}

    features: list[str] = []
    detected_atoms: dict[str, list[int]] = {}

    def add_feature(label: str, atom_indices: list[int] | None = None) -> None:
        if label not in features:
            features.append(label)
        if atom_indices is not None:
            detected_atoms[label] = sorted(set(atom_indices))

    # Core functional motifs
    aldehyde_matches = _matches(mol, "[CX3H1](=O)[#6]")
    ketone_matches = _matches(mol, "[#6][CX3](=O)[#6]")
    ester_matches = _matches(mol, "[#6][CX3](=O)[OX2H0][#6]")
    acid_matches = _matches(mol, "[CX3](=O)[OX2H1,OX1-]")
    amide_matches = _matches(mol, "[NX3][CX3](=[OX1])[#6]")
    nitrile_matches = _matches(mol, "[CX2]#[NX1]")
    ether_matches = _matches(mol, "[OD2]([#6])[#6]")
    alcohol_matches = _matches(mol, "[OX2H][C;!$(C=O)]")
    phenol_matches = _matches(mol, "[OX2H][c]")
    # Exclude nitro nitrogen.  Its neutral-valence representation can satisfy
    # a broad NX3 amine pattern, but it is not an amine functional group.
    amine_matches = _matches(mol, "[NX3;H2,H1,H0;!$(NC=O);!$([N+](=O)[O-])]")
    siloxane_matches = _matches(mol, "[Si][OX2][Si]")
    silanol_matches = _matches(mol, "[Si][OX2H]")
    silyl_alkyl_peroxide_matches = _matches(mol, "[#6][OX2][OX2][Si]")
    phosphate_matches = _matches(mol, "[PX4](=O)([OX2H,OX2H0])([OX2H,OX2H0])[OX2H,OX2H0]")
    # P(=O) with three direct carbon substituents. This excludes phosphate
    # esters, which carry oxygen substituents on phosphorus.
    trialkylphosphine_oxide_matches = _matches(mol, "[PX4](=[OX1])([#6])([#6])[#6]")
    organic_hydroperoxide_matches = _matches(mol, "[#6][OX2][OX2H]")
    dialkyl_peroxide_matches = _matches(mol, "[#6][OX2][OX2][#6]")
    sulfonamide_matches = _matches(mol, "[SX4](=O)(=O)[NX3]")
    sulfonic_acid_matches = _matches(mol, "[SX4](=O)(=O)[OX2H,OX1-]")
    halogen_matches = _matches(mol, "[#9,#17,#35,#53]")

    if aldehyde_matches:
        add_feature("Aldehyde", [idx for match in aldehyde_matches for idx in match])
    if ketone_matches:
        add_feature("Ketone", [idx for match in ketone_matches for idx in match])
    if ester_matches:
        add_feature("Ester", [idx for match in ester_matches for idx in match])
    if acid_matches:
        add_feature("Carboxylic acid", [idx for match in acid_matches for idx in match])
    if amide_matches:
        add_feature("Amide", [idx for match in amide_matches for idx in match])
    if nitrile_matches:
        add_feature("Nitrile", [idx for match in nitrile_matches for idx in match])
    if ether_matches:
        add_feature("Ether", [idx for match in ether_matches for idx in match])
    if alcohol_matches:
        add_feature("Alcohol", [idx for match in alcohol_matches for idx in match])
    if phenol_matches:
        add_feature("Phenol", [idx for match in phenol_matches for idx in match])
    if amine_matches:
        add_feature("Amine", [idx for match in amine_matches for idx in match])
    if siloxane_matches:
        add_feature("Siloxane", [idx for match in siloxane_matches for idx in match])
    if silanol_matches:
        add_feature("Silanol", [idx for match in silanol_matches for idx in match])
    if silyl_alkyl_peroxide_matches:
        add_feature("Silyl alkyl peroxide", [idx for match in silyl_alkyl_peroxide_matches for idx in match])

    alcohol_cyclic = False
    alcohol_ring_substituted = False
    for match in alcohol_matches:
        if len(match) < 2:
            continue
        alcohol_carbon_idx = match[1]
        alcohol_carbon = mol.GetAtomWithIdx(alcohol_carbon_idx)
        if alcohol_carbon.IsInRing():
            alcohol_cyclic = True
            continue
        if any(neighbor.IsInRing() for neighbor in alcohol_carbon.GetNeighbors()):
            alcohol_ring_substituted = True
    if alcohol_cyclic:
        add_feature("Cyclic alcohol")
    if alcohol_ring_substituted:
        add_feature("Ring-substituted alcohol")
    if phosphate_matches:
        add_feature("Phosphate ester", [idx for match in phosphate_matches for idx in match])
    if trialkylphosphine_oxide_matches:
        add_feature("Trialkylphosphine oxide", [idx for match in trialkylphosphine_oxide_matches for idx in match])
    if organic_hydroperoxide_matches:
        add_feature("Organic hydroperoxide", [idx for match in organic_hydroperoxide_matches for idx in match])
    if dialkyl_peroxide_matches:
        add_feature("Dialkyl peroxide", [idx for match in dialkyl_peroxide_matches for idx in match])
    if sulfonamide_matches:
        add_feature("Sulfonamide", [idx for match in sulfonamide_matches for idx in match])
    if sulfonic_acid_matches:
        add_feature("Sulfonic acid", [idx for match in sulfonic_acid_matches for idx in match])
    if halogen_matches:
        add_feature("Halogenated compound", [idx for match in halogen_matches for idx in match])

    # Topology-sensitive refinements
    if any(mol.GetAtomWithIdx(match[1]).IsInRing() for match in ketone_matches if len(match) > 1):
        add_feature("Cyclic ketone")
    if any(mol.GetAtomWithIdx(match[0]).IsInRing() for match in aldehyde_matches if len(match) > 0):
        add_feature("Cyclic aldehyde")
    if any(mol.GetAtomWithIdx(match[1]).IsInRing() or mol.GetAtomWithIdx(match[2]).IsInRing() for match in ether_matches if len(match) > 2):
        add_feature("Cyclic ether")
    if ether_matches and len(ether_matches) >= int(thresholds.get("polyether_min_ether_bridges", 2)):
        add_feature("Polyether")
        if ring_count and max(len(ring) for ring in atom_rings) >= int(thresholds.get("macrocyclic_min_ring_size", 12)) and len(ether_matches) >= int(thresholds.get("macrocyclic_min_ether_bridges", 3)):
            add_feature("Macrocyclic polyether")
    # A cyclic siloxane must have silicon and oxygen in the same perceived
    # ring.  A ring elsewhere in a molecule must not refine a linear Si-O-Si
    # chain into a cyclic siloxane.
    siloxane_rings = [
        ring for ring in atom_rings
        if any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 14 for atom_idx in ring)
        and any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8 for atom_idx in ring)
    ]
    if siloxane_matches and siloxane_rings:
        siloxane_ring_atoms = sorted({atom_idx for ring in siloxane_rings for atom_idx in ring})
        add_feature("Cyclic siloxane", siloxane_ring_atoms)
        if max(len(ring) for ring in siloxane_rings) >= int(thresholds.get("macrocyclic_min_ring_size", 12)):
            add_feature("Macrocyclic siloxane", siloxane_ring_atoms)
    # A 1,2-oxasilinane is a six-membered saturated ring containing adjacent
    # Si and O atoms, with the other ring members carbon. It is not a
    # siloxane because it lacks a Si-O-Si linkage.
    oxasilolane_rings = []
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(atom_idx) for atom_idx in ring]
        if sum(atom.GetAtomicNum() == 14 for atom in atoms) != 1 or sum(atom.GetAtomicNum() == 8 for atom in atoms) != 1:
            continue
        if any(atom.GetAtomicNum() not in {6, 8, 14} for atom in atoms):
            continue
        if not any(
            {atoms[position].GetAtomicNum(), atoms[(position + 1) % len(atoms)].GetAtomicNum()} == {8, 14}
            for position in range(len(atoms))
        ):
            continue
        oxasilolane_rings.append(ring)
    if oxasilolane_rings:
        add_feature("1,2-Oxasilinane", sorted({atom_idx for ring in oxasilolane_rings for atom_idx in ring}))
    if acid_matches and any(mol.GetAtomWithIdx(match[0]).IsInRing() for match in acid_matches if len(match) > 0):
        add_feature("Lactone")
    if ester_matches and any(mol.GetAtomWithIdx(match[1]).IsInRing() for match in ester_matches if len(match) > 1):
        add_feature("Lactone")
    if amide_matches and any(mol.GetAtomWithIdx(match[1]).IsInRing() for match in amide_matches if len(match) > 1):
        add_feature("Lactam")
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        carbonyl_neighbors = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            double_bond_oxygen = any(
                bond.GetBondType() == Chem.BondType.DOUBLE and other.GetAtomicNum() == 8
                for bond in neighbor.GetBonds()
                for other in (bond.GetBeginAtom(), bond.GetEndAtom())
                if other.GetIdx() != neighbor.GetIdx()
            )
            if double_bond_oxygen:
                carbonyl_neighbors.append(neighbor.GetIdx())
        if len(carbonyl_neighbors) >= 2:
            add_feature("Imide")
            break

    acetal_candidates: list[int] = []
    ketal_candidates: list[int] = []
    cyclic_acetal = False
    cyclic_ketal = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0:
            continue
        o_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
        if len(o_neighbors) < 2:
            continue
        # Acetals and ketals require two ether-like O substituents.  A carbon
        # bearing two hydroxyl O atoms is a geminal diol, not an acetal/ketal.
        if any(nbr.GetTotalNumHs() > 0 for nbr in o_neighbors):
            continue
        if atom.GetTotalNumHs() > 0:
            acetal_candidates.append(atom.GetIdx())
            if atom.IsInRing() or any(nbr.IsInRing() for nbr in o_neighbors):
                cyclic_acetal = True
        else:
            ketal_candidates.append(atom.GetIdx())
            if atom.IsInRing() or any(nbr.IsInRing() for nbr in o_neighbors):
                cyclic_ketal = True
    if acetal_candidates:
        add_feature("Acetal", acetal_candidates)
    if ketal_candidates:
        add_feature("Ketal", ketal_candidates)
    if cyclic_acetal:
        add_feature("Cyclic acetal")
    if cyclic_ketal:
        add_feature("Cyclic ketal")

    if phosphate_matches:
        phosphate_alkoxy = 0
        phosphate_atoms = set()
        for match in phosphate_matches:
            phosphate_atoms.update(match)
            p_idx = match[0]
            p_atom = mol.GetAtomWithIdx(p_idx)
            oxygens = [nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            phosphate_alkoxy = max(
                phosphate_alkoxy,
                sum(
                    1 for oxygen in oxygens
                    if any(
                        nei.GetAtomicNum() == 6 and not nei.GetIsAromatic()
                        and nei.GetHybridization() == Chem.HybridizationType.SP3
                        for nei in oxygen.GetNeighbors() if nei.GetIdx() != p_idx
                    )
                ),
            )
        if phosphate_alkoxy >= 3:
            add_feature("Trialkyl phosphate", sorted(phosphate_atoms))

    if any(atom.GetAtomicNum() == 14 and atom.GetDegree() == 2 for atom in mol.GetAtoms()):
        if siloxane_matches:
            add_feature("Siloxane")
    if silanol_matches:
        add_feature("Silanol")

    # A compound is a hydrocarbon only when every heavy atom is carbon. A
    # missing SMARTS match must never turn an O/N/P/S/Si compound into a
    # hydrocarbon. Broad elemental fallbacks retain valid chemistry without
    # pretending to know an unsupported functional group.
    specific_features = {
        "Aldehyde", "Ketone", "Ester", "Carboxylic acid", "Amide", "Ether", "Nitrile",
        "Sulfonamide", "Sulfonic acid", "Siloxane", "Cyclic siloxane", "Macrocyclic siloxane", "Silanol", "1,2-Oxasilinane", "Silyl alkyl peroxide", "Phosphate ester",
        "Trialkylphosphine oxide", "Organic hydroperoxide", "Dialkyl peroxide", "Lactone", "Lactam", "Imide",
        "Acetal", "Ketal", "Cyclic acetal", "Cyclic ketal", "Phenol", "Alcohol", "Amine",
        "Trialkyl phosphate", "Polyether", "Macrocyclic polyether",
    }
    heavy_atomic_numbers = {atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1}
    if not any(feature in features for feature in specific_features):
        if heavy_atomic_numbers and heavy_atomic_numbers == {6}:
            add_feature("Hydrocarbon")
        elif 15 in heavy_atomic_numbers:
            add_feature("Organophosphorus compound — not further classified")
        elif 16 in heavy_atomic_numbers:
            add_feature("Sulfur-containing compound — not further classified")
        elif 14 in heavy_atomic_numbers:
            add_feature("Silicon-containing compound — not further classified")
        elif 7 in heavy_atomic_numbers:
            add_feature("Nitrogen-containing compound — not further classified")
        elif 8 in heavy_atomic_numbers:
            add_feature("Oxygen-containing compound — not further classified")

    if halogen_matches and "Hydrocarbon" not in features:
        add_feature("Halogenated compound")

    # Support metadata
    return {
        "features": features,
        "feature_atom_indices": detected_atoms,
        "ring_count": ring_count,
        "aromatic_ring_count": aromatic_ring_count,
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "molecular_weight": round(Descriptors.MolWt(mol), 3),
        "molecular_refractivity": round(Crippen.MolMR(mol), 3),
    }


def _select_primary_class(feature_profile: dict[str, Any]) -> tuple[str, str, str]:
    features = set(feature_profile["features"])
    checks: list[tuple[str, str, str]] = [
        ("STC-SIL-004", "Macrocyclic siloxane", "Macrocyclic silicon-oxygen ring detected"),
        ("STC-SIL-003", "Cyclic siloxane", "Cyclic silicon-oxygen ring detected"),
        ("STC-SIL-006", "Silyl alkyl peroxide", "Carbon-oxygen-oxygen-silicon motif detected"),
        ("STC-SIL-005", "1,2-Oxasilinane", "Six-membered adjacent silicon-oxygen ring detected"),
        ("STC-SIL-001", "Siloxane", "Siloxane linkage detected"),
        ("STC-SIL-002", "Silanol", "Silanol motif detected"),
        ("STC-PHO-001", "Trialkyl phosphate", "Trialkyl phosphate motif detected"),
        ("STC-PHO-002", "Phosphate ester", "Phosphate ester motif detected"),
        ("STC-PHO-003", "Trialkylphosphine oxide", "Trialkylphosphine oxide motif detected"),
        ("STC-HYP-001", "Organic hydroperoxide", "Organic hydroperoxide motif detected"),
        ("STC-PER-001", "Dialkyl peroxide", "Dialkyl peroxide motif detected"),
        ("STC-SUL-001", "Sulfonic acid", "Sulfonic acid motif detected"),
        ("STC-SUL-002", "Sulfonamide", "Sulfonamide motif detected"),
        ("STC-IMD-001", "Imide", "Imide motif detected"),
        ("STC-LAC-001", "Lactam", "Lactam motif detected"),
        ("STC-AMD-001", "Amide", "Amide motif detected"),
        ("STC-LAC-002", "Lactone", "Lactone motif detected"),
        ("STC-ACD-001", "Carboxylic acid", "Carboxylic acid motif detected"),
        ("STC-ACE-001", "Cyclic acetal", "Cyclic acetal motif detected"),
        ("STC-ACE-002", "Cyclic ketal", "Cyclic ketal motif detected"),
        ("STC-ACE-003", "Acetal", "Acetal motif detected"),
        ("STC-ACE-004", "Ketal", "Ketal motif detected"),
        ("STC-EST-001", "Ester", "Ester motif detected"),
        ("STC-KET-001", "Aromatic ketone", "Aromatic ketone motif detected"),
        ("STC-KET-002", "Cyclic ketone", "Cyclic ketone motif detected"),
        ("STC-KET-003", "Ketone", "Ketone motif detected"),
        ("STC-ALD-001", "Aldehyde", "Aldehyde motif detected"),
        ("STC-ETH-001", "Macrocyclic polyether", "Macrocyclic polyether motif detected"),
        ("STC-ETH-002", "Polyether", "Polyether motif detected"),
        ("STC-ETH-003", "Cyclic ether", "Cyclic ether motif detected"),
        ("STC-ETH-004", "Ether", "Ether motif detected"),
        ("STC-NIT-001", "Nitrile", "Nitrile motif detected"),
        ("STC-PHE-001", "Phenol", "Phenol motif detected"),
        ("STC-ALC-001", "Alcohol", "Alcohol motif detected"),
        ("STC-AMI-001", "Amine", "Amine motif detected"),
        ("STC-HAL-001", "Halogenated compound", "Halogenated compound fallback"),
        ("STC-P-UNK-001", "Organophosphorus compound — not further classified", "Phosphorus present but no approved specific motif matched"),
        ("STC-S-UNK-001", "Sulfur-containing compound — not further classified", "Sulfur present but no approved specific motif matched"),
        ("STC-SI-UNK-001", "Silicon-containing compound — not further classified", "Silicon present but no approved specific motif matched"),
        ("STC-N-UNK-001", "Nitrogen-containing compound — not further classified", "Nitrogen present but no approved specific motif matched"),
        ("STC-O-UNK-001", "Oxygen-containing compound — not further classified", "Oxygen present but no approved specific motif matched"),
        ("STC-HYD-001", "Hydrocarbon", "Hydrocarbon fallback"),
    ]
    for rule_id, label, reason in checks:
        if label in features:
            return rule_id, label, reason
    return "STC-UNK-001", "Unclassified", "No priority rule matched"


def classification_scope(rule_id: str) -> tuple[str, bool, str]:
    """Describe whether a result is specific, provisional, or unavailable."""
    if rule_id == "STC-UNK-000":
        return "Structure unavailable", False, "Manual review required before classification"
    if rule_id == "STC-UNK-001":
        return "Unclassified", True, "No approved structural rule matched; review recommended"
    if rule_id in BROAD_PROVISIONAL_RULE_IDS:
        return "Broad provisional class", True, "Valid structure has no approved specific rule; review recommended"
    return "Specific rule", False, ""


def classification_suggested_action(rule_id: str) -> str:
    """Return the next safe curation action without changing the decision."""
    if rule_id == "STC-UNK-000":
        return "Resolve identity or provide a valid structure before classification"
    if rule_id == "STC-UNK-001":
        return "Chemistry review required; consider an approved new structural rule"
    if rule_id in BROAD_PROVISIONAL_RULE_IDS:
        return "Review the feature evidence and curate a specific structural rule if warranted"
    return "No classification action required"


def classification_secondary_functional_groups(classification: dict[str, Any]) -> tuple[str, str]:
    """Return independent secondary features without double-counting SMARTS overlap."""
    feature_names = list(classification.get("detected_feature_profile") or [])
    feature_profile = classification.get("feature_profile") or {}
    feature_atoms = feature_profile.get("feature_atom_indices") or {}
    hierarchy = classification.get("taxonomy_hierarchy") or {}
    primary_names = {
        str(classification.get("corrected_chemical_class", "")),
        str(hierarchy.get("Functional Group", "")),
    }
    primary_atoms: set[int] = set()
    for name in primary_names:
        primary_atoms.update(feature_atoms.get(name, []))
    secondary_features: list[str] = []
    for name in feature_names:
        if name in primary_names:
            continue
        atom_indices = set(feature_atoms.get(name, []))
        if not atom_indices or atom_indices.intersection(primary_atoms):
            continue
        secondary_features.append(name)
    if not secondary_features:
        return "Not available", "No non-overlapping secondary feature detected"
    return "; ".join(secondary_features), "Non-overlapping detected feature atoms"


def classification_export_fields(classification: dict[str, Any]) -> dict[str, Any]:
    """Flatten the deterministic classification contract for dataframe exports.

    This is deliberately shared by Phase A paths so corrected class, taxonomy,
    topology, and ontology values cannot be assembled differently by parallel
    and sequential workflows.
    """
    hierarchy = classification.get("taxonomy_hierarchy") or {}
    topology = classification.get("topology_profile") or {}
    ontology = classification.get("ontology_compatibility") or {}
    secondary_groups, secondary_basis = classification_secondary_functional_groups(classification)
    scope, review_recommended, recommendation = classification_scope(
        str(classification.get("classification_rule_id", ""))
    )
    fields = {
        "Corrected Chemical Class": classification.get("corrected_chemical_class", "Unclassified"),
        # Keep PRISM's deterministic structural interpretation visibly
        # separate from an exact external/reference taxonomy match.  Neither
        # is silently substituted for the other in the exported workbook.
        "PRISM Structural Class": classification.get("corrected_chemical_class", "Unclassified"),
        "Final Classification Record": classification.get("final_classification_record", ""),
        "Classification Input Source": classification.get("classification_input_source", ""),
        "Classification Structure SMILES": classification.get("classification_structure_smiles", ""),
        "Classification Standardization Version": classification.get("classification_standardization_version", ""),
        "Classification Stereochemistry Status": classification.get("classification_stereochemistry_status", ""),
        "Classification Status": classification.get("classification_status", ""),
        "Classification Scope": scope,
        "Classification Review Recommended": review_recommended,
        "Classification Review Recommendation": recommendation,
        "Classification Suggested Action": classification.get("classification_suggested_action", ""),
        "Classification Rule ID": classification.get("classification_rule_id", ""),
        "Classification Rule Version": classification.get("classification_rule_version", ""),
        "Manual Review Flag": bool(classification.get("manual_review_flag", False)),
        "Manual Review Reason": classification.get("manual_review_reason", ""),
        "Detected Feature Profile": "; ".join(classification.get("detected_feature_profile") or []),
        "Primary Functional Group": hierarchy.get("Functional Group", ""),
        "Secondary Functional Groups": secondary_groups,
        "Secondary Functional Groups Basis": secondary_basis,
        "Taxonomy Path": classification.get("taxonomy_path", ""),
        "Taxonomy Path Steps": " → ".join(classification.get("taxonomy_path_steps") or []),
        "Matched Categories": " → ".join(classification.get("matched_categories") or []),
        "Direct Parent": classification.get("direct_parent", ""),
        "Taxonomy Dictionary Version": classification.get("taxonomy_dictionary_version", ""),
        "Taxonomy Hierarchy": " | ".join(f"{key}: {value}" for key, value in hierarchy.items()),
        "Topology Profile": " | ".join(f"{key}: {value}" for key, value in topology.items()),
        "Ontology Mapping Status": ontology.get("mapping_status", ""),
        "Ontology Mapping Version": ontology.get("ontology_mapping_version", ""),
        "Ontology Kingdom": ontology.get("kingdom", ""),
        "Ontology Superclass": ontology.get("superclass", ""),
        "Ontology Class": ontology.get("class", ""),
        "Ontology Subclass": ontology.get("subclass", ""),
        "Ontology Direct Parent": ontology.get("direct_parent", ""),
        "Ontology Mapping Evidence": ontology.get("mapping_evidence", ""),
        "Local Taxonomy Match Status": ontology.get("mapping_status", "") if ontology.get("source") == "Chemical Taxonomies local SQLite" else "Not matched",
        "Local Taxonomy Source": ontology.get("source", ""),
        "Local Taxonomy Record ID": ontology.get("record_id", ""),
        "Local Taxonomy Path": ontology.get("taxonomy_path", ""),
        "Reference Taxonomy Path": ontology.get("taxonomy_path", ""),
        "Local Taxonomy Reference": ontology.get("reference", ""),
        "Local Taxonomy URL": ontology.get("url", ""),
    }
    return fields


def apply_local_taxonomy_evidence(classification: dict[str, Any], inchikey: object) -> dict[str, Any]:
    """Attach an exact offline taxonomy match without replacing PRISM's rule result.

    A reference match is molecular-identity evidence. PRISM's SMARTS-derived
    corrected class, functional group, and review flag remain authoritative
    because the source taxonomy does not encode this application's
    endpoint-specific interpretation policy.
    """
    evidence = ontology_evidence(inchikey)
    if evidence is None:
        return classification
    enriched = dict(classification)
    existing = dict(classification.get("ontology_compatibility") or {})
    existing.update(evidence)
    # Preserve the local-rule context while exposing the exact reference
    # hierarchy in the same audited ontology block.
    existing["local_rule_id"] = str(classification.get("classification_rule_id") or "")
    existing["local_class"] = str(classification.get("corrected_chemical_class") or "")
    enriched["ontology_compatibility"] = existing
    return enriched


def _row_inchikey(row: pd.Series, classification: dict[str, Any]) -> str:
    """Prefer supplied identity evidence, then derive a key from valid structure."""
    for column in ("InChIKey", "InChiKey", "inchikey"):
        value = str(row.get(column, "") or "").strip()
        if value:
            return value
    smiles = str(classification.get("classification_structure_smiles", "") or "").strip()
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return ""
    try:
        return Chem.MolToInchiKey(molecule) or ""
    except Exception:
        return ""


def apply_classification_contract(
    dataframe: pd.DataFrame,
    standardized_col: str = "Standardized SMILES",
    original_col: str = "Original SMILES",
) -> pd.DataFrame:
    """Refresh shared classification columns from the current local rules.

    This is safe for Phase B resume/retry because it preserves all non-
    classification columns and deterministically replaces only the contract
    fields derived from the available structure.
    """
    result = dataframe.copy()
    records = []
    for _, row in result.iterrows():
        classification = classify_structure(
            row.get(standardized_col, ""), row.get(original_col, "")
        )
        classification = apply_local_taxonomy_evidence(
            classification, _row_inchikey(row, classification)
        )
        records.append(classification_export_fields(classification))
    if not records:
        return result
    contract = pd.DataFrame(records, index=result.index)
    for column in contract.columns:
        result[column] = contract[column]
    return result


def _taxonomy_family(corrected_class: str) -> str:
    family_map = {
        "Siloxane": "Silicon-containing Compound",
        "Cyclic siloxane": "Silicon-containing Compound",
        "Macrocyclic siloxane": "Silicon-containing Compound",
        "1,2-Oxasilinane": "Silicon-containing Compound",
        "Silyl alkyl peroxide": "Silicon-containing Compound",
        "Silanol": "Silicon-containing Compound",
        "Trialkyl phosphate": "Phosphorus-containing Compound",
        "Phosphate ester": "Phosphorus-containing Compound",
        "Trialkylphosphine oxide": "Phosphorus-containing Compound",
        "Organophosphorus compound — not further classified": "Phosphorus-containing Compound",
        "Sulfonic acid": "Sulfur-containing Compound",
        "Sulfonamide": "Sulfur-containing Compound",
        "Imide": "Nitrogen-containing Compound",
        "Lactam": "Nitrogen-containing Compound",
        "Amide": "Nitrogen-containing Compound",
        "Amine": "Nitrogen-containing Compound",
        "Nitrile": "Nitrogen-containing Compound",
        "Aldehyde": "Oxygen-containing Compound",
        "Ketone": "Oxygen-containing Compound",
        "Cyclic ketone": "Oxygen-containing Compound",
        "Aromatic ketone": "Oxygen-containing Compound",
        "Ester": "Oxygen-containing Compound",
        "Lactone": "Oxygen-containing Compound",
        "Carboxylic acid": "Oxygen-containing Compound",
        "Acetal": "Oxygen-containing Compound",
        "Cyclic acetal": "Oxygen-containing Compound",
        "Ketal": "Oxygen-containing Compound",
        "Cyclic ketal": "Oxygen-containing Compound",
        "Ether": "Oxygen-containing Compound",
        "Cyclic ether": "Oxygen-containing Compound",
        "Polyether": "Oxygen-containing Compound",
        "Macrocyclic polyether": "Oxygen-containing Compound",
        "Alcohol": "Oxygen-containing Compound",
        "Phenol": "Oxygen-containing Compound",
        "Organic hydroperoxide": "Oxygen-containing Compound",
        "Dialkyl peroxide": "Oxygen-containing Compound",
        "Oxygen-containing compound — not further classified": "Oxygen-containing Compound",
        "Nitrogen-containing compound — not further classified": "Nitrogen-containing Compound",
        "Sulfur-containing compound — not further classified": "Sulfur-containing Compound",
        "Silicon-containing compound — not further classified": "Silicon-containing Compound",
        "Halogenated compound": "Halogenated Compound",
        "Hydrocarbon": "Hydrocarbon",
        "Unclassified": "Unclassified",
    }
    return family_map.get(corrected_class, "Unclassified")


def _alcohol_topology_label(topology_profile: dict[str, Any]) -> str:
    topology_class = topology_profile.get("Topology Class", "Acyclic")
    if topology_class == "Macrocyclic":
        return "Macrocyclic Alcohol"
    if topology_class == "Fused Polycyclic":
        return "Fused Polycyclic Alcohol"
    if topology_class == "Spirocyclic":
        return "Spirocyclic Alcohol"
    if topology_class == "Polycyclic":
        return "Polycyclic Alcohol"
    if topology_class == "Bicyclic":
        return "Bicyclic Alcohol"
    if topology_class == "Monocyclic":
        return "Monocyclic Alcohol"
    return "Acyclic Alcohol"


def _taxonomy_path_for_class(corrected_class: str, feature_profile: dict[str, Any], topology_profile: dict[str, Any] | None = None) -> list[str]:
    features = set(feature_profile.get("features", []))
    family = _taxonomy_family(corrected_class)
    topology_profile = topology_profile or {}

    if corrected_class == "Unclassified":
        return ["Unclassified"]

    generic_fallbacks = {
        "Organophosphorus compound — not further classified",
        "Sulfur-containing compound — not further classified",
        "Silicon-containing compound — not further classified",
        "Nitrogen-containing compound — not further classified",
        "Oxygen-containing compound — not further classified",
    }
    if corrected_class in generic_fallbacks:
        return [family, "Not further classified"]

    if family == "Silicon-containing Compound":
        if corrected_class == "1,2-Oxasilinane":
            return [family, "Cyclic Organosilicon Ether", "1,2-Oxasilinane"]
        if corrected_class == "Silyl alkyl peroxide":
            return [family, "Organosilicon Peroxide", "Silyl Alkyl Peroxide"]
        if corrected_class == "Macrocyclic siloxane":
            return [family, "Siloxane", "Cyclic Siloxane", "Macrocyclic Siloxane"]
        if corrected_class == "Cyclic siloxane":
            return [family, "Siloxane", "Cyclic Siloxane"]
        if corrected_class == "Siloxane":
            return [family, "Siloxane"]
        if corrected_class == "Silanol":
            return [family, "Silanol"]
        return [family, corrected_class]

    if family == "Phosphorus-containing Compound":
        if corrected_class == "Trialkyl phosphate":
            return [family, "Phosphate Ester", "Trialkyl Phosphate"]
        if corrected_class == "Phosphate ester":
            return [family, "Phosphate Ester"]
        if corrected_class == "Trialkylphosphine oxide":
            return [family, "Phosphine Oxide", "Trialkylphosphine Oxide"]
        return [family, corrected_class]

    if family == "Sulfur-containing Compound":
        if corrected_class == "Sulfonamide":
            return [family, "Sulfonamide"]
        if corrected_class == "Sulfonic acid":
            return [family, "Sulfonic Acid"]
        return [family, corrected_class]

    if family == "Nitrogen-containing Compound":
        if corrected_class == "Lactam":
            return [family, "Amide", "Lactam", "Cyclic Lactam"]
        if corrected_class == "Imide":
            return [family, "Amide", "Imide"]
        if corrected_class == "Amide":
            return [family, "Amide"]
        if corrected_class == "Amine":
            return [family, "Amine"]
        if corrected_class == "Nitrile":
            return [family, "Nitrile"]
        return [family, corrected_class]

    if family == "Halogenated Compound":
        return [family]

    if family == "Hydrocarbon":
        if feature_profile.get("aromatic_ring_count", 0):
            return [family, "Aromatic Hydrocarbon"]
        topology_class = topology_profile.get("Topology Class", "Acyclic")
        if topology_class == "Monocyclic":
            return [family, "Monocyclic Hydrocarbon"]
        if topology_class in {"Bicyclic", "Polycyclic", "Fused Polycyclic", "Spirocyclic"}:
            return [family, f"{topology_class} Hydrocarbon"]
        return [family, "Acyclic Hydrocarbon"]

    if family == "Oxygen-containing Compound":
        if corrected_class in {"Cyclic ketal", "Ketal", "Cyclic acetal", "Acetal"}:
            path = [family, "Acetal"]
            if corrected_class in {"Ketal", "Cyclic ketal"} or "Ketal" in features:
                path.append("Ketal")
            if corrected_class == "Cyclic ketal" or "Cyclic ketal" in features:
                path.append("Cyclic Ketal")
            elif corrected_class == "Cyclic acetal" or "Cyclic acetal" in features:
                path.append("Cyclic Acetal")
            return path
        if corrected_class in {"Cyclic ether", "Ether", "Polyether", "Macrocyclic polyether"}:
            path = [family, "Ether"]
            if corrected_class == "Polyether" or "Polyether" in features:
                path.append("Polyether")
            if corrected_class == "Macrocyclic polyether" or "Macrocyclic polyether" in features:
                path.append("Macrocyclic Polyether")
            if corrected_class == "Cyclic ether" or "Cyclic ether" in features:
                path.append("Cyclic Ether")
            return path
        if corrected_class in {"Cyclic ketone", "Ketone", "Aromatic ketone"}:
            path = [family, "Ketone"]
            if corrected_class == "Aromatic ketone" or "Aromatic ketone" in features:
                path.append("Aromatic Ketone")
            if corrected_class == "Cyclic ketone" or "Cyclic ketone" in features:
                path.append("Cyclic Ketone")
            return path
        if corrected_class == "Lactone":
            return [family, "Ester", "Lactone", "Cyclic Lactone"]
        if corrected_class == "Ester":
            return [family, "Ester"]
        if corrected_class == "Carboxylic acid":
            if "Aromatic Ring" in features or feature_profile.get("aromatic_ring_count", 0):
                path = [family, "Carboxylic Acid", "Aromatic Carboxylic Acid"]
                if feature_profile.get("ring_count", 0) == 1:
                    path.append("Monocyclic Aromatic Carboxylic Acid")
                else:
                    path.append("Aromatic Carboxylic Acid")
                return path
            return [family, "Carboxylic Acid"]
        if corrected_class == "Aldehyde":
            return [family, "Aldehyde"]
        if corrected_class == "Alcohol":
            alcohol_topology = _alcohol_topology_label(topology_profile)
            alcohol_substitution = topology_profile.get("Alcohol Substitution", "Not applicable")
            path = [family, "Alcohol", alcohol_topology]
            if alcohol_substitution != "Not applicable":
                path.append(alcohol_substitution)
            return path
        if corrected_class == "Phenol":
            return [family, "Phenol"]
        if corrected_class == "Organic hydroperoxide":
            return [family, "Hydroperoxide", "Organic Hydroperoxide"]
        if corrected_class == "Dialkyl peroxide":
            return [family, "Peroxide", "Organic Peroxide", "Dialkyl Peroxide"]
        return [family, corrected_class]

    return [family, corrected_class] if family != corrected_class else [family]


def _taxonomy_hierarchy_for_class(corrected_class: str, feature_profile: dict[str, Any], topology_profile: dict[str, Any] | None = None) -> dict[str, str]:
    features = set(feature_profile.get("features", []))
    family = _taxonomy_family(corrected_class)
    topology_profile = topology_profile or {}

    hierarchy = {
        "Parent Class": family,
        "Functional Group": corrected_class,
        "Subclass": corrected_class,
        "Structural Type": corrected_class,
    }

    if corrected_class == "Unclassified":
        return {
            "Parent Class": "Unclassified",
            "Functional Group": "Unclassified",
            "Subclass": "Unclassified",
            "Structural Type": "Unclassified",
        }

    generic_fallbacks = {
        "Organophosphorus compound — not further classified",
        "Sulfur-containing compound — not further classified",
        "Silicon-containing compound — not further classified",
        "Nitrogen-containing compound — not further classified",
        "Oxygen-containing compound — not further classified",
    }
    if corrected_class in generic_fallbacks:
        hierarchy["Functional Group"] = "Not further classified"
        hierarchy["Subclass"] = "Not further classified"
        hierarchy["Structural Type"] = "Not further subclassified"
        return hierarchy

    if family == "Phosphorus-containing Compound":
        hierarchy["Functional Group"] = "Phosphate Ester"
        if corrected_class == "Trialkyl phosphate":
            hierarchy["Subclass"] = "Phosphate Ester"
            hierarchy["Structural Type"] = "Trialkyl Phosphate"
        elif corrected_class == "Phosphate ester":
            hierarchy["Subclass"] = "Phosphate Ester"
        elif corrected_class == "Trialkylphosphine oxide":
            hierarchy["Functional Group"] = "Phosphine Oxide"
            hierarchy["Subclass"] = "Trialkylphosphine Oxide"
            hierarchy["Structural Type"] = "Acyclic Phosphine Oxide"
        return hierarchy

    if family == "Nitrogen-containing Compound":
        if corrected_class == "Lactam":
            hierarchy["Functional Group"] = "Amide"
            hierarchy["Subclass"] = "Lactam"
            hierarchy["Structural Type"] = "Cyclic Lactam"
        elif corrected_class == "Imide":
            hierarchy["Functional Group"] = "Amide"
            hierarchy["Subclass"] = "Imide"
        elif corrected_class == "Amide":
            hierarchy["Functional Group"] = "Amide"
        elif corrected_class == "Amine":
            hierarchy["Functional Group"] = "Amine"
        elif corrected_class == "Nitrile":
            hierarchy["Functional Group"] = "Nitrile"
        return hierarchy

    if family == "Oxygen-containing Compound":
        if corrected_class in {"Cyclic ketal", "Ketal", "Cyclic acetal", "Acetal"}:
            hierarchy["Functional Group"] = "Acetal"
            hierarchy["Subclass"] = "Ketal"
            hierarchy["Structural Type"] = "Cyclic Ketal" if corrected_class == "Cyclic ketal" or "Cyclic ketal" in features else ("Cyclic Acetal" if corrected_class == "Cyclic acetal" or "Cyclic acetal" in features else corrected_class)
        elif corrected_class in {"Cyclic ether", "Ether", "Polyether", "Macrocyclic polyether"}:
            hierarchy["Functional Group"] = "Ether"
            hierarchy["Subclass"] = "Polyether" if corrected_class in {"Polyether", "Macrocyclic polyether"} or "Polyether" in features else "Ether"
            hierarchy["Structural Type"] = "Macrocyclic Polyether" if corrected_class == "Macrocyclic polyether" or "Macrocyclic polyether" in features else ("Cyclic Ether" if corrected_class == "Cyclic ether" or "Cyclic ether" in features else corrected_class)
        elif corrected_class in {"Cyclic ketone", "Ketone", "Aromatic ketone"}:
            hierarchy["Functional Group"] = "Ketone"
            hierarchy["Subclass"] = "Cyclic Ketone" if corrected_class == "Cyclic ketone" or "Cyclic ketone" in features else "Ketone"
            hierarchy["Structural Type"] = "Aromatic Ketone" if corrected_class == "Aromatic ketone" or "Aromatic ketone" in features else hierarchy["Subclass"]
        elif corrected_class == "Lactone":
            hierarchy["Functional Group"] = "Ester"
            hierarchy["Subclass"] = "Lactone"
            hierarchy["Structural Type"] = "Cyclic Lactone"
        elif corrected_class == "Ester":
            hierarchy["Functional Group"] = "Ester"
            hierarchy["Subclass"] = "Ester"
        elif corrected_class == "Carboxylic acid":
            hierarchy["Functional Group"] = "Carboxylic Acid"
            if "Aromatic Ring" in features or feature_profile.get("aromatic_ring_count", 0):
                hierarchy["Subclass"] = "Aromatic Carboxylic Acid"
                hierarchy["Structural Type"] = "Monocyclic Aromatic Carboxylic Acid" if feature_profile.get("ring_count", 0) == 1 else "Aromatic Carboxylic Acid"
        elif corrected_class == "Aldehyde":
            hierarchy["Functional Group"] = "Aldehyde"
        elif corrected_class == "Alcohol":
            hierarchy["Functional Group"] = "Alcohol"
            hierarchy["Subclass"] = _alcohol_topology_label(topology_profile)
            hierarchy["Structural Type"] = topology_profile.get("Alcohol Substitution", "Not applicable")
        elif corrected_class == "Phenol":
            hierarchy["Functional Group"] = "Phenol"
        elif corrected_class == "Organic hydroperoxide":
            hierarchy["Functional Group"] = "Hydroperoxide"
            hierarchy["Subclass"] = "Organic Hydroperoxide"
            hierarchy["Structural Type"] = "Organic Hydroperoxide"
        elif corrected_class == "Dialkyl peroxide":
            hierarchy["Functional Group"] = "Peroxide"
            hierarchy["Subclass"] = "Organic Peroxide"
            hierarchy["Structural Type"] = "Dialkyl Peroxide"
        return hierarchy

    if family == "Hydrocarbon":
        hierarchy["Functional Group"] = "Hydrocarbon"
        hierarchy["Subclass"] = _taxonomy_path_for_class(corrected_class, feature_profile, topology_profile)[-1]
        hierarchy["Structural Type"] = "Not further subclassified"
        return hierarchy

    if family == "Silicon-containing Compound":
        if corrected_class == "1,2-Oxasilinane":
            hierarchy["Functional Group"] = "Cyclic Organosilicon Ether"
            hierarchy["Subclass"] = "1,2-Oxasilinane"
            hierarchy["Structural Type"] = "1,2-Oxasilinane"
        elif corrected_class == "Silyl alkyl peroxide":
            hierarchy["Functional Group"] = "Organosilicon Peroxide"
            hierarchy["Subclass"] = "Silyl Alkyl Peroxide"
            hierarchy["Structural Type"] = "Silyl Alkyl Peroxide"
        elif corrected_class == "Macrocyclic siloxane":
            hierarchy["Functional Group"] = "Siloxane"
            hierarchy["Subclass"] = "Cyclic Siloxane"
            hierarchy["Structural Type"] = "Macrocyclic Siloxane"
        elif corrected_class == "Cyclic siloxane":
            hierarchy["Functional Group"] = "Siloxane"
            hierarchy["Subclass"] = "Cyclic Siloxane"
            hierarchy["Structural Type"] = "Cyclic Siloxane"
        return hierarchy

    if family in {"Sulfur-containing Compound", "Halogenated Compound"}:
        return hierarchy

    if hierarchy["Structural Type"] == hierarchy["Subclass"]:
        hierarchy["Structural Type"] = "Not further subclassified"

    return hierarchy


def classify_structure(standardized_smiles: object, original_smiles: object | None = None, rules_path: str | Path | None = None) -> dict[str, Any]:
    """Return a deterministic local structure classification."""
    rules = load_rules(rules_path)

    input_source = "Standardized SMILES"
    candidate_smiles = standardized_smiles
    mol = _mol_from_smiles(candidate_smiles)
    if mol is None:
        candidate_smiles = original_smiles
        input_source = "Original SMILES" if _mol_from_smiles(candidate_smiles) is not None else "Unavailable"
        mol = _mol_from_smiles(candidate_smiles)

    if mol is None:
        scope, review_recommended, recommendation = classification_scope("STC-UNK-000")
        return {
            "corrected_chemical_class": "Unclassified",
            "final_classification_record": "Taxonomy path: Unclassified; Structure unavailable -> Manual Review",
            "final_classification_details": {
                "Taxonomy Path": "Unclassified",
                "Primary Class": "Unclassified",
                "Rule ID": "STC-UNK-000",
                "Rule Source": "Unavailable",
                "Feature Evidence": "None",
                "Ring Profile": "Not assessable",
                "Topology Profile": "Not assessable",
                "State": "Not assessable",
                "Review": "Manual review required",
            },
            "manual_review_flag": True,
            "manual_review_reason": "No valid standardized or original SMILES could be parsed.",
            "classification_input_source": "Unavailable",
            "classification_structure_smiles": "",
            "classification_standardization_version": STANDARDIZATION_VERSION,
            "classification_stereochemistry_status": "Unavailable",
            "classification_status": "Manual review",
            "classification_rule_id": "STC-UNK-000",
            "classification_rule_version": rules.get("rule_version", "unknown"),
            "classification_scope": scope,
            "classification_review_recommended": review_recommended,
            "classification_review_recommendation": recommendation,
            "classification_suggested_action": classification_suggested_action("STC-UNK-000"),
            "detected_feature_profile": [],
            "taxonomy_path": "Unclassified",
            "taxonomy_path_steps": ["Unclassified"],
            "matched_categories": ["Unclassified"],
            "direct_parent": "",
            "taxonomy_dictionary_version": rules.get("taxonomy_dictionary_version", "unknown"),
            "taxonomy_hierarchy": {
                "Parent Class": "Unclassified",
                "Functional Group": "Unclassified",
                "Subclass": "Unclassified",
                "Structural Type": "Not assessable",
            },
            "topology_profile": {
                "Topology Class": "Unclassified",
                "Topology Modifiers": "Unclassified",
                "Ring System": "Unclassified",
                "Ring Profile": "Not assessable",
                "Ring Count": 0,
                "Aromatic Ring Count": 0,
                "Ring Sizes": [],
                "Hetero Ring Atom Count": 0,
                "Shared Ring Atom Sizes": [],
            },
            "feature_profile_text": "Not assessable",
            "feature_profile": {},
            "ontology_compatibility": {
                "mapping_status": "Not applicable",
                "ontology_mapping_version": load_ontology_compatibility_mappings().get("ontology_mapping_version", "unavailable"),
                "local_rule_id": "STC-UNK-000",
                "local_class": "Unclassified",
                "mapping_evidence": "Structure unavailable",
            },
        }

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    has_defined_stereochemistry = any(
        atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()
    ) or any(bond.GetStereo() != Chem.rdchem.BondStereo.STEREONONE for bond in mol.GetBonds())
    feature_profile = _feature_profile(mol, rules)
    rule_id, corrected_class, rule_reason = _select_primary_class(feature_profile)
    topology_profile = _topology_profile_for_mol(mol, feature_profile)
    taxonomy_path_steps = _taxonomy_path_for_class(corrected_class, feature_profile, topology_profile)
    taxonomy_path = " → ".join(taxonomy_path_steps)
    # The direct parent is the immediately less-specific approved category.
    # It is derived from the same de-duplicated local path that drives the
    # hierarchy, rather than from an optional external ontology crosswalk.
    direct_parent = taxonomy_path_steps[-2] if len(taxonomy_path_steps) > 1 else ""
    taxonomy_hierarchy = _taxonomy_hierarchy_for_class(corrected_class, feature_profile, topology_profile)

    selected_source = input_source if input_source != "Unavailable" else "Original SMILES"
    detected_features = feature_profile["features"]
    aromatic_note = "aromatic" if feature_profile["aromatic_ring_count"] else "non-aromatic"
    scope, review_recommended, recommendation = classification_scope(rule_id)
    review_text = "Manual review required" if corrected_class == "Unclassified" else (
        recommendation if review_recommended else "No manual review required"
    )
    final_record = (
        f"Taxonomy path: {taxonomy_path}; primary class: {corrected_class}; rule: {rule_id}; source: {selected_source}; "
        f"features: {', '.join(detected_features) if detected_features else 'none'}; "
        f"ring profile: {feature_profile['ring_count']} ring(s), {feature_profile['aromatic_ring_count']} aromatic ring(s); "
        f"topology: {topology_profile['Topology Class']}; modifiers: {topology_profile['Topology Modifiers']}; "
        f"state: {aromatic_note}; review: {review_text}"
    )

    if corrected_class == "Unclassified":
        manual_review_flag = True
        manual_review_reason = "No priority class matched the detected feature profile."
        classification_status = "Manual review"
    elif rule_id in BROAD_PROVISIONAL_RULE_IDS:
        manual_review_flag = False
        manual_review_reason = ""
        classification_status = "Calculated - broad provisional class; review recommended"
    else:
        manual_review_flag = False
        manual_review_reason = ""
        classification_status = "Calculated"

    final_classification_details = {
        "Taxonomy Path": taxonomy_path,
        "Primary Class": corrected_class,
        "Rule ID": rule_id,
        "Rule Source": selected_source,
        "Feature Evidence": ", ".join(detected_features) if detected_features else "None",
        "Ring Profile": f"{feature_profile['ring_count']} ring(s), {feature_profile['aromatic_ring_count']} aromatic ring(s)",
        "Topology Profile": f"{topology_profile['Topology Class']}; {topology_profile['Topology Modifiers']}",
        "Alcohol Substitution": topology_profile.get("Alcohol Substitution", "Not applicable"),
        "State": aromatic_note,
        "Review": "Manual review required" if manual_review_flag else "No manual review required",
    }
    ontology_compatibility = _ontology_compatibility(rule_id, corrected_class)

    return {
        "corrected_chemical_class": corrected_class,
        "final_classification_record": final_record,
        "final_classification_details": final_classification_details,
        "manual_review_flag": manual_review_flag,
        "manual_review_reason": manual_review_reason,
        "classification_input_source": selected_source,
        "classification_structure_smiles": canonical_smiles,
        "classification_standardization_version": STANDARDIZATION_VERSION,
        "classification_stereochemistry_status": "Defined stereochemistry retained" if has_defined_stereochemistry else "No defined stereochemistry",
        "classification_status": classification_status,
        "classification_rule_id": rule_id,
        "classification_rule_version": rules.get("rule_version", "unknown"),
        "classification_scope": scope,
        "classification_review_recommended": review_recommended,
        "classification_review_recommendation": recommendation,
        "classification_suggested_action": classification_suggested_action(rule_id),
        "detected_feature_profile": detected_features,
        "taxonomy_path": taxonomy_path,
        "taxonomy_path_steps": taxonomy_path_steps,
        "matched_categories": taxonomy_path_steps,
        "direct_parent": direct_parent,
        "taxonomy_dictionary_version": rules.get("taxonomy_dictionary_version", "unknown"),
        "taxonomy_hierarchy": taxonomy_hierarchy,
        "topology_profile": topology_profile,
        "feature_profile_text": "; ".join(detected_features) if detected_features else "None",
        "feature_profile": feature_profile,
        "ontology_compatibility": ontology_compatibility,
    }


def add_structure_classification(df: pd.DataFrame, standardized_col: str = "Standardized SMILES", original_col: str = "Original SMILES", rules_path: str | Path | None = None) -> pd.DataFrame:
    """Append classification columns to a dataframe without changing legacy fields."""
    result = df.copy()
    records = []
    for _, row in result.iterrows():
        classification = classify_structure(
            row.get(standardized_col, ""), row.get(original_col, ""), rules_path=rules_path
        )
        records.append(apply_local_taxonomy_evidence(classification, _row_inchikey(row, classification)))
    classified = pd.DataFrame(records, index=result.index)
    classified = classified.rename(columns={
        "corrected_chemical_class": "Corrected Chemical Class",
        "final_classification_record": "Final Classification Record",
        "manual_review_flag": "Manual Review Flag",
        "manual_review_reason": "Manual Review Reason",
        "classification_input_source": "Classification Input Source",
        "classification_status": "Classification Status",
        "classification_rule_id": "Classification Rule ID",
        "classification_rule_version": "Classification Rule Version",
        "detected_feature_profile": "Detected Feature Profile",
        "taxonomy_path": "Taxonomy Path",
        "taxonomy_path_steps": "Taxonomy Path Steps",
        "taxonomy_hierarchy": "Taxonomy Hierarchy",
        "topology_profile": "Topology Profile",
        "ontology_compatibility": "Ontology Compatibility",
    })
    if "Detected Feature Profile" in classified.columns:
        classified["Detected Feature Profile"] = classified["Detected Feature Profile"].apply(
            lambda value: "; ".join(value) if isinstance(value, list) else value
        )
    if "Taxonomy Path Steps" in classified.columns:
        classified["Taxonomy Path Steps"] = classified["Taxonomy Path Steps"].apply(
            lambda value: " → ".join(value) if isinstance(value, list) else value
        )
    if "Ontology Compatibility" in classified.columns:
        ontology_fields = {
            "mapping_status": "Ontology Mapping Status",
            "ontology_mapping_version": "Ontology Mapping Version",
            "kingdom": "Ontology Kingdom",
            "superclass": "Ontology Superclass",
            "class": "Ontology Class",
            "subclass": "Ontology Subclass",
            "direct_parent": "Ontology Direct Parent",
            "mapping_evidence": "Ontology Mapping Evidence",
        }
        ontology_values = classified["Ontology Compatibility"].apply(
            lambda value: value if isinstance(value, dict) else {}
        )
        for source_field, output_column in ontology_fields.items():
            classified[output_column] = ontology_values.apply(
                lambda value, field=source_field: value.get(field, "")
            )
        classified = classified.drop(columns=["Ontology Compatibility"])
    if "Topology Profile" in classified.columns:
        classified["Topology Profile"] = classified["Topology Profile"].apply(
            lambda value: " | ".join(f"{k}: {v}" for k, v in value.items()) if isinstance(value, dict) else value
        )
    return pd.concat([result, classified], axis=1)
