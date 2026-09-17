"""Review-only descriptor extension for clustering.

This layer expands the clustering evidence with the additional physicochemical
descriptors already calculated elsewhere in the app. It remains soft and
auditable: it can influence the similarity matrix, but it never hard-gates a
cluster assignment.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_DESCRIPTOR_LAYER_RULES_PATH = Path(__file__).resolve().parent / "config" / "clustering_descriptor_layer_rules.json"
_MISSING_TEXT = {"", "none", "n/a", "not assessable", "unknown", "not implemented"}


def _rule_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_DESCRIPTOR_LAYER_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "clustering_descriptor_layer_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_descriptor_layer_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _rule_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                rules = json.load(handle)
            if rules.get("mode") != "review_only":
                raise ValueError("Descriptor layer rules must remain review_only until separately approved.")
            if not str(rules.get("rule_version", "")).strip():
                raise ValueError("Descriptor layer rules must define rule_version.")
            if not rules.get("descriptor_columns"):
                raise ValueError("Descriptor layer rules must define descriptor_columns.")
            return rules
    raise FileNotFoundError("clustering_descriptor_layer_rules.json not found")


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.casefold() in _MISSING_TEXT:
        return ""
    return text


def _available_columns(df: pd.DataFrame, rules: dict[str, Any]) -> list[str]:
    return [column for column in rules["descriptor_columns"] if column in df.columns]


def _normalize_numeric(values: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    valid = numeric.dropna()
    if valid.empty:
        return np.zeros(len(numeric), dtype=float)
    minimum = float(valid.min())
    maximum = float(valid.max())
    if np.isclose(minimum, maximum):
        return np.zeros(len(numeric), dtype=float)
    return ((numeric.fillna(minimum) - minimum) / (maximum - minimum)).to_numpy(dtype=float)


def _pairwise_numeric_distance(values: pd.Series) -> np.ndarray:
    normalized = _normalize_numeric(values)
    return np.abs(normalized[:, np.newaxis] - normalized[np.newaxis, :])


def _pairwise_ionisation_distance(values: pd.Series, rules: dict[str, Any]) -> np.ndarray:
    params = rules["ionisation_distance"]
    cleaned = [_clean_text(value) or "Unknown" for value in values.tolist()]
    matrix = np.zeros((len(cleaned), len(cleaned)), dtype=float)
    for i in range(len(cleaned)):
        for j in range(i):
            left = cleaned[i]
            right = cleaned[j]
            if left == right and left not in {"Unknown", ""}:
                distance = float(params["same_known"])
            elif left == right == "Unknown":
                distance = float(params["both_unknown"])
            elif "Unknown" in {left, right}:
                distance = float(params["one_unknown"])
            else:
                distance = float(params["different_known"])
            matrix[i, j] = distance
            matrix[j, i] = distance
    return matrix


def build_descriptor_distance_matrix(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> tuple[np.ndarray, dict[str, Any]]:
    """Build a review-only distance matrix from the expanded descriptor layer."""
    rules = rules or load_descriptor_layer_rules()
    available_columns = _available_columns(df, rules)
    if len(df) < 2 or len(available_columns) < 1:
        metadata = {
            "rule_version": rules["rule_version"],
            "available_columns": available_columns,
            "coverage_ratio": 0.0,
            "status": "Not assessable",
        }
        return np.zeros((len(df), len(df)), dtype=float), metadata

    matrix = np.zeros((len(df), len(df)), dtype=float)
    used_weights: dict[str, float] = {}
    for column in available_columns:
        weight = float(rules["descriptor_weights"].get(column, 0.0))
        if weight <= 0:
            continue
        used_weights[column] = weight
        if column in rules.get("numeric_descriptor_columns", []):
            column_distance = _pairwise_numeric_distance(df[column])
        else:
            column_distance = _pairwise_ionisation_distance(df[column], rules)
        matrix += weight * column_distance

    total_weight = sum(used_weights.values())
    if total_weight > 0:
        matrix = matrix / total_weight

    row_coverage = []
    for _, row in df.iterrows():
        available = sum(1 for column in available_columns if _clean_text(row.get(column, "")) != "")
        row_coverage.append(round(available / len(available_columns), 3) if available_columns else 0.0)

    metadata = {
        "rule_version": rules["rule_version"],
        "available_columns": available_columns,
        "coverage_ratio": round(len(available_columns) / len(rules["descriptor_columns"]), 3),
        "row_coverage": row_coverage,
        "status": "Calculated" if available_columns else "Not assessable",
    }
    return matrix, metadata


def annotate_descriptor_layer_evidence(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    """Add cluster-level summary evidence for the expanded descriptor layer."""
    result = df.copy()
    rules = rules or load_descriptor_layer_rules()
    result["Descriptor_Layer_Status"] = "Not assessable"
    result["Descriptor_Layer_Cohesion"] = "Not assessable"
    result["Descriptor_Layer_Reasons"] = ""
    result["Descriptor_Layer_Coverage"] = np.nan

    if "Cluster ID" not in result.columns:
        return result

    matrix, metadata = build_descriptor_distance_matrix(result, rules)
    result["Descriptor_Layer_Coverage"] = metadata.get("row_coverage", [0.0] * len(result))

    thresholds = rules.get("cluster_cohesion_thresholds", {})
    perfect_threshold = float(thresholds.get("perfect_mean_pairwise_distance", 0.05))
    aligned_threshold = float(thresholds.get("aligned_mean_pairwise_distance", 0.25))

    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            continue

        if len(indices) == 1 or not metadata.get("available_columns"):
            result.loc[indices, "Descriptor_Layer_Status"] = "Not assessable"
            result.loc[indices, "Descriptor_Layer_Cohesion"] = "Not assessable"
            result.loc[indices, "Descriptor_Layer_Reasons"] = "Singleton cluster or no usable descriptor columns"
            continue

        cluster_matrix = matrix[np.ix_(indices, indices)]
        pairwise = cluster_matrix[np.tril_indices(len(indices), k=-1)]
        pairwise = pairwise[np.isfinite(pairwise)]
        if pairwise.size == 0:
            mean_pairwise = 0.0
        else:
            mean_pairwise = float(np.mean(pairwise))

        mean_coverage = float(np.mean(result.loc[indices, "Descriptor_Layer_Coverage"].fillna(0.0)))
        if mean_pairwise <= perfect_threshold:
            cohesion = "Chemically coherent"
        elif mean_pairwise <= aligned_threshold:
            cohesion = "Chemically aligned"
        else:
            cohesion = "Chemically mixed"

        result.loc[indices, "Descriptor_Layer_Status"] = "Calculated"
        result.loc[indices, "Descriptor_Layer_Cohesion"] = cohesion
        result.loc[indices, "Descriptor_Layer_Reasons"] = (
            f"Mean pairwise descriptor distance {mean_pairwise:.3f}; "
            f"coverage {mean_coverage:.0%}; "
            f"available descriptors: {', '.join(metadata['available_columns'])}"
        )

    return result
