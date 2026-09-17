"""Auditable, reviewer-facing summaries for revalidated structural clusters."""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


CLUSTER_SUMMARY_COLUMNS = [
    "Compatibility Group",
    "Cluster ID",
    "Compatibility Gate",
    "Cluster Status",
    "Source Row Count",
    "Unique Structure Count",
    "Duplicate or Alias Row Count",
    "Representative Compound",
    "Representative Standardized SMILES",
    "Member Chemical Classes",
    "Primary Functional Group Profile",
    "Minimum Pairwise Tanimoto",
    "Median Pairwise Tanimoto",
    "Nearest Neighbour Tanimoto Range",
    "Scaffold Profile",
    "Ionisation Profile",
    "Structural Alert Profile",
    "Property Compatibility",
    "Property Compatibility Reason",
    "Identity Structure Consistency",
    "Compatibility Gate Reason",
    "Membership Rationale",
    "Stability Assessment",
    "SME Review Status",
    "SME Decision",
    "SME Review Comments",
    "Assessment Boundary",
]


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "not available", "not assessable"} else text


def _unique_values(frame: pd.DataFrame, column: str, limit: int = 6) -> str:
    if column not in frame.columns:
        return "Not assessed"
    values = []
    for value in frame[column].tolist():
        text = _text(value)
        if text and text not in values:
            values.append(text)
    if not values:
        return "Not assessed"
    suffix = "; additional values present" if len(values) > limit else ""
    return "; ".join(values[:limit]) + suffix


def _numeric_range(frame: pd.DataFrame, column: str) -> str:
    if column not in frame.columns:
        return "Not assessable"
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if values.empty:
        return "Not assessable"
    return f"{values.min():.3f} to {values.max():.3f}"


def _representative(frame: pd.DataFrame, smiles_col: str = "Standardized SMILES") -> tuple[str, str]:
    unique = frame.drop_duplicates(smiles_col, keep="first") if smiles_col in frame.columns else frame.iloc[0:0]
    if unique.empty:
        return "Not assessable", "Not assessable"
    valid = []
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    for index, row in unique.iterrows():
        smiles = _text(row.get(smiles_col))
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is not None:
            valid.append((index, generator.GetFingerprint(molecule)))
    if not valid:
        row = unique.iloc[0]
        return _text(row.get("Compound Name")) or "Unnamed compound", _text(row.get(smiles_col)) or "Not assessable"
    if len(valid) == 1:
        index = valid[0][0]
    else:
        scores = []
        for current, fingerprint in valid:
            similarities = [DataStructs.TanimotoSimilarity(fingerprint, other) for _, other in valid]
            scores.append((float(np.mean(similarities)), current))
        index = max(scores)[1]
    row = unique.loc[index]
    return _text(row.get("Compound Name")) or "Unnamed compound", _text(row.get(smiles_col)) or "Not assessable"


def build_cluster_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Create one explicit SME-review record per calculated cluster."""
    if df.empty or "Cluster ID" not in df.columns:
        return pd.DataFrame(columns=CLUSTER_SUMMARY_COLUMNS)
    source = df.copy()
    if "Compatibility Group" not in source.columns and "Cluster Family" in source.columns:
        source = source.rename(columns={"Cluster Family": "Compatibility Group"})
    clustered = source[
        source["Cluster ID"].astype(str).str.startswith("Cluster_", na=False)
    ].copy()
    if clustered.empty:
        return pd.DataFrame(columns=CLUSTER_SUMMARY_COLUMNS)

    records = []
    for cluster_id, group in clustered.groupby("Cluster ID", sort=True):
        unique_structures = group.drop_duplicates("Standardized SMILES", keep="first") if "Standardized SMILES" in group.columns else group
        representative_name, representative_smiles = _representative(group)
        # ``DataFrame.get`` returns a scalar when the column is absent, which
        # cannot be treated as a Series.  Sparse/legacy exports legitimately
        # omit these optional diagnostics, so retain an empty Series instead.
        median_values = (
            group["Cluster_Median_Pairwise_Tanimoto"]
            if "Cluster_Median_Pairwise_Tanimoto" in group.columns
            else pd.Series(dtype=float)
        )
        min_similarity_values = (
            group["Cluster_Min_Pairwise_Tanimoto"]
            if "Cluster_Min_Pairwise_Tanimoto" in group.columns
            else pd.Series(dtype=float)
        )
        median = pd.to_numeric(median_values, errors="coerce").dropna()
        min_similarity = pd.to_numeric(min_similarity_values, errors="coerce").dropna()
        source_count = len(group)
        unique_count = len(unique_structures)
        review_reason = _unique_values(group, "Compatibility Gate Reason", limit=1)
        rationale = _unique_values(group, "Cluster_Membership_Rationale", limit=1)
        status = _unique_values(group, "Cluster_Status", limit=1)
        stability = "Not assessable"
        if unique_count == 1:
            stability = "Not assessable: one unique structure"
        elif not median.empty and median.iloc[0] >= 0.60:
            stability = "Structurally coherent under primary threshold; sensitivity review pending"
        else:
            stability = "Borderline: SME review required before use"
        records.append({
            "Compatibility Group": _unique_values(group, "Compatibility Group", limit=1),
            "Cluster ID": cluster_id,
            "Compatibility Gate": _unique_values(group, "Compatibility Gate", limit=1),
            "Cluster Status": status,
            "Source Row Count": source_count,
            "Unique Structure Count": unique_count,
            "Duplicate or Alias Row Count": source_count - unique_count,
            "Representative Compound": representative_name,
            "Representative Standardized SMILES": representative_smiles,
            "Member Chemical Classes": _unique_values(group, "Corrected Chemical Class"),
            "Primary Functional Group Profile": _unique_values(group, "Primary Functional Group"),
            "Minimum Pairwise Tanimoto": "Not assessable" if min_similarity.empty else round(float(min_similarity.iloc[0]), 3),
            "Median Pairwise Tanimoto": "Not assessable" if median.empty else round(float(median.iloc[0]), 3),
            "Nearest Neighbour Tanimoto Range": _numeric_range(unique_structures, "Nearest_Neighbour_Tanimoto"),
            "Scaffold Profile": _unique_values(unique_structures, "Scaffold"),
            "Ionisation Profile": _unique_values(unique_structures, "Ionisation_Indicator"),
            "Structural Alert Profile": _unique_values(unique_structures, "All_Structural_Alerts"),
            "Property Compatibility": _unique_values(group, "Cluster Property Compatibility", limit=1),
            "Property Compatibility Reason": _unique_values(group, "Cluster Property Compatibility Reason", limit=1),
            "Identity Structure Consistency": _unique_values(group, "Identity Structure Consistency"),
            "Compatibility Gate Reason": review_reason,
            "Membership Rationale": rationale,
            "Stability Assessment": stability,
            "SME Review Status": "Pending",
            "SME Decision": "Pending",
            "SME Review Comments": "",
            "Assessment Boundary": "Structural grouping only. SME review is required before analogue selection, read-across, or endpoint-specific assessment.",
        })
    return pd.DataFrame(records, columns=CLUSTER_SUMMARY_COLUMNS)
