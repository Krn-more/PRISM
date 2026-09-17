"""Deterministic, compound-specific evidence for explaining cluster membership.

This module never changes cluster membership, scores, or decisions.  It makes
the already-calculated grouping evidence explicit so a reviewer (and an
optional AI narrative) can explain why a particular row is in a cluster.
"""

from __future__ import annotations

import json
from typing import Any

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


RULE_VERSION = "1.0.0-deterministic-membership"
_FINGERPRINT = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _text(value: Any, default: str = "Not supplied") -> str:
    if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
        return default
    value = str(value).strip()
    return value or default


def _identifier(row: pd.Series) -> str:
    name = _text(row.get("Compound Name"), "Unnamed compound")
    key = _text(row.get("InChIKey"), "")
    return f"{name} [{key}]" if key else name


def _fingerprint(smiles: Any):
    try:
        molecule = Chem.MolFromSmiles(_text(smiles, ""))
        return _FINGERPRINT.GetFingerprint(molecule) if molecule is not None else None
    except Exception:
        return None


def _sme_boundary(row: pd.Series) -> str:
    triggers: list[str] = []
    domain_status = _text(row.get("Domain_Status"), "Not assessable")
    if domain_status.casefold() != "inside":
        triggers.append(f"domain status is {domain_status}")
    cluster_status = _text(row.get("Cluster_Status"), "")
    if "unstable" in cluster_status.casefold():
        triggers.append("cluster stability requires review")
    property_flags = _text(row.get("Property_Outlier_Flags"), "None")
    if property_flags.casefold().startswith("none; not assessed"):
        triggers.append("some cluster property measures are not assessed")
    elif property_flags.casefold() not in {"none", "not supplied", "not assessable"}:
        triggers.append("property-outlier flags are present")
    ionisation = _text(row.get("Ionisation_Consistency"), "Consistent").casefold()
    if ionisation not in {"consistent", "not assessable", "not supplied"}:
        triggers.append("ionisation differs from the cluster majority")
    toxicophore = _text(row.get("Toxicophore_Consistency"), "Consistent").casefold()
    if toxicophore not in {"consistent", "not assessable", "not supplied"}:
        triggers.append("toxicophore profile differs from the cluster majority")
    if _text(row.get("Identity Status"), "Resolved").casefold() not in {"resolved", "exact match"}:
        triggers.append("identity status requires confirmation")
    reason = "; ".join(triggers) if triggers else "the grouping is a technical similarity result, not an endpoint conclusion"
    return (
        "SME review is required before this grouping is used for analogue selection, "
        "read-across, or any endpoint-specific assessment because " + reason + "."
    )


def add_cluster_membership_evidence(df: pd.DataFrame, smiles_col: str = "Standardized SMILES") -> pd.DataFrame:
    """Append deterministic per-row evidence without recalculating clusters."""
    result = df.copy()
    fields = [
        "Cluster_Reasoning_Status",
        "Cluster_Nearest_Members",
        "Cluster_Nearest_Tanimoto",
        "Cluster_Consensus_Method",
        "Cluster_Membership_Parameters",
        "Cluster_Membership_Rationale",
        "SME_Review_Boundary",
        "Cluster_Reasoning_Rule_Version",
    ]
    for field in fields:
        result[field] = ""
    if result.empty or "Cluster ID" not in result.columns or smiles_col not in result.columns:
        return result

    fingerprints = {index: _fingerprint(result.at[index, smiles_col]) for index in result.index}
    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        assessable = str(cluster_id).startswith("Cluster_")
        for index in indices:
            row = result.loc[index]
            fingerprint = fingerprints[index]
            if not assessable or fingerprint is None:
                result.at[index, "Cluster_Reasoning_Status"] = "Not assessable"
                result.at[index, "Cluster_Membership_Rationale"] = (
                    "No cluster-membership rationale is available because a valid standardized structure and "
                    "an assessable cluster assignment are required."
                )
                result.at[index, "SME_Review_Boundary"] = _sme_boundary(row)
                result.at[index, "Cluster_Reasoning_Rule_Version"] = RULE_VERSION
                continue

            neighbours: list[tuple[float, Any]] = []
            for other_index in indices:
                if other_index == index or fingerprints[other_index] is None:
                    continue
                neighbours.append((DataStructs.TanimotoSimilarity(fingerprint, fingerprints[other_index]), other_index))
            neighbours.sort(key=lambda item: (-item[0], str(item[1])))
            top_neighbours = neighbours[:3]
            nearest_text = "; ".join(
                f"{_identifier(result.loc[other_index])} (Tanimoto {score:.3f})"
                for score, other_index in top_neighbours
            ) or "No in-cluster structural neighbour (singleton cluster)"
            nearest_score = f"{top_neighbours[0][0]:.3f}" if top_neighbours else "Not available"
            scaffold = _text(row.get("Scaffold"), "Not supplied")
            representative = _text(row.get("Cluster_Representative_Scaffold"), "Not supplied")
            scaffold_match = "yes" if scaffold != "Not supplied" and scaffold == representative else "no"
            consensus_method = _text(
                row.get("Cluster_Consensus_Method"),
                "Compatibility-gated Morgan/ECFP4 Tanimoto Butina clustering with complete-linkage refinement",
            )
            membership_parameters = _text(
                row.get("Cluster_Membership_Parameters"),
                "Morgan radius=2; fpSize=1024; Tanimoto >= 0.60; complete-linkage refinement",
            )
            parameters = {
                "cluster_id": str(cluster_id),
                "cluster_size": _text(row.get("Cluster Size")),
                "consensus_method": consensus_method,
                "membership_parameters": membership_parameters,
                "nearest_members": nearest_text,
                "nearest_member_tanimoto": nearest_score,
                "representative_scaffold": representative,
                "representative_scaffold_match": scaffold_match,
                "scaffold_coverage": _text(row.get("Cluster_Scaffold_Coverage")),
                "domain_status": _text(row.get("Domain_Status")),
                "ionisation_consistency": _text(row.get("Ionisation_Consistency")),
                "toxicophore_consistency": _text(row.get("Toxicophore_Consistency")),
                "property_outlier_flags": _text(row.get("Property_Outlier_Flags")),
                "classification": _text(row.get("Corrected Chemical Class")),
            }
            rationale = (
                f"Assigned to {cluster_id} by {consensus_method}. "
                f"Cluster size: {parameters['cluster_size']}. Nearest in-cluster member(s): {nearest_text}. "
                f"Representative scaffold match: {scaffold_match}; shared-scaffold coverage: {parameters['scaffold_coverage']}. "
                f"Domain status: {parameters['domain_status']}; ionisation: {parameters['ionisation_consistency']}; "
                f"toxicophore profile: {parameters['toxicophore_consistency']}; "
                f"property-outlier flags: {parameters['property_outlier_flags']}."
            )
            result.at[index, "Cluster_Reasoning_Status"] = "Calculated deterministic evidence"
            result.at[index, "Cluster_Nearest_Members"] = nearest_text
            result.at[index, "Cluster_Nearest_Tanimoto"] = nearest_score
            result.at[index, "Cluster_Consensus_Method"] = parameters["consensus_method"]
            result.at[index, "Cluster_Membership_Parameters"] = json.dumps(parameters, sort_keys=True)
            result.at[index, "Cluster_Membership_Rationale"] = rationale
            result.at[index, "SME_Review_Boundary"] = _sme_boundary(row)
            result.at[index, "Cluster_Reasoning_Rule_Version"] = RULE_VERSION
    return result


def compound_reasoning_packet(row: dict[str, Any]) -> dict[str, Any]:
    """Return an allow-listed, non-authoritative packet for an AI narrative."""
    fields = (
        "Compound Name", "CASRN", "InChIKey", "Standardized SMILES", "Corrected Chemical Class",
        "Primary Functional Group", "Taxonomy Path", "Cluster ID", "Cluster Size", "Cluster_Status",
        "Cluster_Reasoning_Status",
        "Cluster_Membership_Rationale", "Cluster_Membership_Parameters", "Cluster_Nearest_Members",
        "Cluster_Nearest_Tanimoto", "Cluster_Consensus_Method", "Domain_Status", "Domain_Reasons",
        "Ionisation_Consistency", "Toxicophore_Consistency", "Property_Outlier_Flags",
        "Uncertainty_Summary", "SME_Review_Boundary",
    )
    packet = {field: _text(row.get(field), "Not supplied") for field in fields}
    return {"evidence_fields": list(packet), "evidence": packet, "rule_version": RULE_VERSION}
