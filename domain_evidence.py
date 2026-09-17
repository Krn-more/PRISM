"""Stage 2 deterministic, review-only grouping-domain evidence.

This module must run after consensus clustering. It appends evidence columns and
does not alter existing cluster membership, `Domain`, or `Decision` values.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


DEFAULT_RULES_PATH = Path(__file__).with_name("config") / "stage2_domain_rules.json"
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_DOMAIN_STATUSES = {"Inside", "Borderline", "Outside", "Not assessable"}


def load_domain_rules(path: str | Path | None = None) -> dict[str, Any]:
    """Load the versioned review-only domain rules."""
    with open(path or DEFAULT_RULES_PATH, encoding="utf-8") as handle:
        rules = json.load(handle)
    if rules.get("mode") != "review_only":
        raise ValueError("Stage 2 rules must remain review_only until separately approved.")
    if rules.get("descriptor_outlier_method") != "iqr_1.5":
        raise ValueError("Stage 2 currently supports only the iqr_1.5 descriptor outlier method.")
    for key in ("minimum_nearest_neighbour_similarity", "scaffold_coverage_threshold"):
        if not isinstance(rules.get(key), (int, float)) or not 0 <= rules[key] <= 1:
            raise ValueError(f"Stage 2 rule {key} must be a number between 0 and 1.")
    if rules.get("singleton_status") not in _DOMAIN_STATUSES:
        raise ValueError("Stage 2 singleton_status must be a recognised domain status.")
    if rules.get("missing_structure_status") not in _DOMAIN_STATUSES:
        raise ValueError("Stage 2 missing_structure_status must be a recognised domain status.")
    return rules


def _fingerprint(smiles: object):
    if smiles is None or pd.isna(smiles) or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    return _MORGAN.GetFingerprint(mol) if mol is not None else None


def _is_missing_scalar(value: object) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _column_series(frame: pd.DataFrame, column: str) -> pd.Series | None:
    """Return a single Series even when the source frame contains duplicate labels."""
    if column not in frame.columns:
        return None
    subset = frame.loc[:, frame.columns == column]
    if isinstance(subset, pd.Series):
        return subset
    if subset.shape[1] == 1:
        return subset.iloc[:, 0]

    def _pick_first_non_missing(row: pd.Series) -> object:
        for value in row.tolist():
            if not _is_missing_scalar(value) and str(value).strip():
                return value
        for value in row.tolist():
            if not _is_missing_scalar(value):
                return value
        return np.nan

    return subset.apply(_pick_first_non_missing, axis=1)


def _column_value(frame: pd.DataFrame, index: int, column: str, default: object = None) -> object:
    series = _column_series(frame, column)
    if series is None or index not in series.index:
        return default
    value = series.at[index]
    return default if _is_missing_scalar(value) else value


def _profile(value: object) -> set[str]:
    if value is None or pd.isna(value):
        return set()
    text = str(value).strip()
    if not text or text in {"None", "Not assessable", "N/A"}:
        return set()
    return {item.strip() for item in text.split(";") if item.strip()}


def _unique_mode(values: list[Any]) -> Any | None:
    """Return the sole most-common value, or None when no unique majority exists."""
    counts = Counter(values)
    if not counts:
        return None
    highest_count = max(counts.values())
    modes = [value for value, count in counts.items() if count == highest_count]
    return modes[0] if len(modes) == 1 else None


def _identity_evidence(value: object, rules: dict[str, Any]) -> tuple[str | None, str | None]:
    """Return an overriding status/reason only when an identity status is supplied."""
    if value is None or pd.isna(value) or not str(value).strip():
        return None, None
    status = str(value).strip()
    normalised = status.casefold()
    if any(term.casefold() in normalised for term in rules["identity_not_assessable_terms"]):
        return "Not assessable", f"Identity status prevents assessment: {status}"
    if any(term.casefold() in normalised for term in rules["identity_review_terms"]):
        return "Borderline", f"Identity status requires review: {status}"
    return None, None


def _descriptor_outlier_values(series: pd.Series, multiplier: float) -> set[Any]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if len(values) < 4:
        return set()
    q1, q3 = values.quantile(0.25), values.quantile(0.75)
    iqr = q3 - q1
    if iqr == 0:
        return set()
    return set(values[(values < q1 - multiplier * iqr) | (values > q3 + multiplier * iqr)].index)


def _new_evidence_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    result["Domain_Status"] = "Not assessable"
    result["Domain_Reasons"] = ""
    result["Domain_Rule_Version"] = ""
    result["Nearest_Neighbour_Tanimoto"] = np.nan
    result["Cluster_Min_Pairwise_Tanimoto"] = np.nan
    result["Cluster_Median_Pairwise_Tanimoto"] = np.nan
    result["Cluster_Scaffold_Coverage"] = np.nan
    result["Cluster_Representative_Scaffold"] = ""
    result["Property_Outlier_Flags"] = ""
    result["Ionisation_Consistency"] = "Not assessable"
    result["Toxicophore_Consistency"] = "Not assessable"
    return result


def add_domain_evidence(df: pd.DataFrame, smiles_col: str = "Standardized SMILES", rules_path: str | Path | None = None) -> pd.DataFrame:
    """Append per-member Stage 2 evidence grouped by existing `Cluster ID`."""
    rules = load_domain_rules(rules_path)
    result = _new_evidence_columns(df)
    if result.empty or "Cluster ID" not in result.columns or smiles_col not in result.columns:
        result["Domain_Reasons"] = "Cluster ID or standardized structure unavailable"
        result["Domain_Rule_Version"] = rules["rule_version"]
        return result

    result["Domain_Rule_Version"] = rules["rule_version"]
    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            result.loc[indices, "Domain_Reasons"] = "No cluster assignment"
            continue

        fingerprints = {idx: _fingerprint(result.at[idx, smiles_col]) for idx in indices}
        valid = [idx for idx in indices if fingerprints[idx] is not None]
        invalid = [idx for idx in indices if fingerprints[idx] is None]
        if invalid:
            result.loc[invalid, "Domain_Status"] = rules["missing_structure_status"]
            result.loc[invalid, "Domain_Reasons"] = "Standardized structure missing or invalid"
        if not valid:
            continue

        scaffolds = [str(result.at[idx, "Scaffold"]).strip() for idx in valid if "Scaffold" in result.columns and pd.notna(result.at[idx, "Scaffold"]) and str(result.at[idx, "Scaffold"]).strip()]
        scaffold_counts = Counter(scaffolds)
        if scaffold_counts:
            highest_scaffold_count = max(scaffold_counts.values())
            representative_scaffold = sorted(scaffold for scaffold, count in scaffold_counts.items() if count == highest_scaffold_count)[0]
            scaffold_count = highest_scaffold_count
        else:
            representative_scaffold, scaffold_count = "Not supplied", 0
        coverage = scaffold_count / len(valid)
        result.loc[valid, "Cluster_Representative_Scaffold"] = representative_scaffold
        result.loc[valid, "Cluster_Scaffold_Coverage"] = round(coverage, 3)

        pairwise, nearest = [], {idx: [] for idx in valid}
        for position, left in enumerate(valid):
            for right in valid[:position]:
                similarity = DataStructs.TanimotoSimilarity(fingerprints[left], fingerprints[right])
                pairwise.append(similarity)
                nearest[left].append(similarity)
                nearest[right].append(similarity)
        if pairwise:
            result.loc[valid, "Cluster_Min_Pairwise_Tanimoto"] = round(min(pairwise), 3)
            result.loc[valid, "Cluster_Median_Pairwise_Tanimoto"] = round(float(np.median(pairwise)), 3)
        for idx in valid:
            if nearest[idx]:
                result.at[idx, "Nearest_Neighbour_Tanimoto"] = round(max(nearest[idx]), 3)

        property_outliers: dict[Any, list[str]] = {idx: [] for idx in valid}
        property_not_assessed: list[str] = []
        for column in rules["descriptor_columns"]:
            series = _column_series(result, column)
            if series is None:
                property_not_assessed.append(column)
                continue
            valid_series = series.loc[valid]
            if pd.to_numeric(valid_series, errors="coerce").notna().sum() < 4:
                property_not_assessed.append(column)
                continue
            for idx in _descriptor_outlier_values(valid_series, rules["descriptor_outlier_threshold"]):
                property_outliers[idx].append(column)

        profiles = {idx: _profile(_column_value(result, idx, "Toxicophore_Profile")) for idx in valid}
        common_profile_tuple = _unique_mode([tuple(sorted(profile)) for profile in profiles.values()])
        common_profile = set(common_profile_tuple) if common_profile_tuple is not None else None
        ionisations = {idx: str(_column_value(result, idx, "Ionisation_Indicator", "Unknown")) for idx in valid}
        common_ionisation = _unique_mode(list(ionisations.values()))

        for idx in valid:
            identity_status = _column_value(result, idx, rules["identity_status_column"])
            reasons = [
                f"Representative scaffold {representative_scaffold}; coverage {coverage:.0%}",
                f"Nearest-neighbour Tanimoto {result.at[idx, 'Nearest_Neighbour_Tanimoto']:.3f}" if pd.notna(result.at[idx, "Nearest_Neighbour_Tanimoto"]) else "Singleton cluster; no pairwise similarity",
            ]
            review_flags = []
            if len(valid) == 1:
                review_flags.append("Singleton cluster requires analogue review")
            if len(valid) > 1 and result.at[idx, "Nearest_Neighbour_Tanimoto"] < rules["minimum_nearest_neighbour_similarity"]:
                review_flags.append("Nearest-neighbour similarity below review threshold")
            if coverage < rules["scaffold_coverage_threshold"]:
                review_flags.append("Shared-scaffold coverage below review threshold")
            if property_outliers[idx]:
                review_flags.append("Property outlier: " + ", ".join(property_outliers[idx]))
            if common_profile is None:
                review_flags.append("Toxicophore profile has no unique cluster majority")
            elif profiles[idx] != common_profile:
                review_flags.append("Toxicophore profile differs from cluster majority")
            if common_ionisation is None:
                review_flags.append("Ionisation indicator has no unique cluster majority")
            elif ionisations[idx] != common_ionisation:
                review_flags.append("Ionisation indicator differs from cluster majority")

            property_status = "; ".join(property_outliers[idx]) if property_outliers[idx] else "None"
            if property_not_assessed:
                property_status += "; Not assessed (fewer than 4 valid values or unavailable): " + ", ".join(property_not_assessed)
            result.at[idx, "Property_Outlier_Flags"] = property_status
            result.at[idx, "Ionisation_Consistency"] = "No unique cluster majority" if common_ionisation is None else ("Consistent" if ionisations[idx] == common_ionisation else "Different from cluster majority")
            result.at[idx, "Toxicophore_Consistency"] = "No unique cluster majority" if common_profile is None else ("Consistent" if profiles[idx] == common_profile else "Different from cluster majority")
            identity_domain_status, identity_reason = _identity_evidence(identity_status, rules)
            status = identity_domain_status or (rules["singleton_status"] if len(valid) == 1 else ("Borderline" if review_flags else "Inside"))
            result.at[idx, "Domain_Status"] = status
            result.at[idx, "Domain_Reasons"] = "; ".join(reasons + review_flags + ([identity_reason] if identity_reason else []))
    return result
