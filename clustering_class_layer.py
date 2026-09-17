"""Review-only chemical-class soft layer for clustering.

This module stays lightweight and avoids RDKit so the class-aware layer can be
tested independently from the full clustering engine.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import numpy as np


DEFAULT_CLASS_LAYER_RULES_PATH = Path(__file__).resolve().parent / "config" / "stage2_class_layer_rules.json"
_MISSING_MARKERS = {"", "none", "n/a", "not assessable", "unclassified", "not further subclassified", "unknown"}


def _rule_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_CLASS_LAYER_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "stage2_class_layer_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_class_layer_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _rule_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                rules = json.load(handle)
            if rules.get("mode") != "review_only":
                raise ValueError("Chemical class clustering rules must remain review_only until separately approved.")
            if not str(rules.get("rule_version", "")).strip():
                raise ValueError("Chemical class clustering rules must define rule_version.")
            weights = rules.get("hierarchy_distance_weights", {})
            required = {
                "same_structural_type",
                "same_subclass",
                "same_functional_group",
                "same_parent_class",
                "same_corrected_class",
                "missing_class_information",
                "different_family",
            }
            if not required.issubset(weights):
                raise ValueError("Chemical class clustering rules are missing hierarchy distance weights.")
            soft_layer_weight = rules.get("soft_layer_weight")
            if not isinstance(soft_layer_weight, (int, float)) or not 0 <= soft_layer_weight <= 1:
                raise ValueError("Chemical class clustering soft_layer_weight must be between 0 and 1.")
            return rules
    raise FileNotFoundError("stage2_class_layer_rules.json not found")


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.casefold() in _MISSING_MARKERS:
        return ""
    return re.sub(r"\s+", " ", text)


def _parse_hierarchy_string(value: object) -> dict[str, str]:
    text = _clean_text(value)
    if not text:
        return {}
    parsed: dict[str, str] = {}
    for part in text.split("|"):
        if ":" not in part:
            continue
        key, raw_value = part.split(":", 1)
        key = _clean_text(key)
        raw_value = _clean_text(raw_value)
        if key and raw_value:
            parsed[key] = raw_value
    return parsed


def extract_class_signature(row: pd.Series | dict[str, Any]) -> dict[str, str]:
    """Extract a normalised taxonomy signature from a row."""
    if isinstance(row, pd.Series):
        getter = row.get
    else:
        getter = row.get

    hierarchy: dict[str, str] = {}
    hierarchy_value = getter("Taxonomy Hierarchy", None)
    if isinstance(hierarchy_value, dict):
        hierarchy = {str(key): _clean_text(value) for key, value in hierarchy_value.items()}
    else:
        hierarchy = _parse_hierarchy_string(hierarchy_value)

    tax_path = _clean_text(getter("Taxonomy Path", None))
    path_parts = [part.strip() for part in tax_path.split("→") if part.strip()] if tax_path else []

    parent_class = _clean_text(
        hierarchy.get("Parent Class")
        or getter("Parent Class", None)
        or (path_parts[0] if len(path_parts) >= 1 else "")
    )
    functional_group = _clean_text(
        hierarchy.get("Functional Group")
        or getter("Functional Group", None)
        or getter("Corrected Chemical Class", None)
        or getter("Chemical Class", None)
        or (path_parts[1] if len(path_parts) >= 2 else "")
    )
    subclass = _clean_text(
        hierarchy.get("Subclass")
        or getter("Subclass", None)
        or (path_parts[2] if len(path_parts) >= 3 else "")
    )
    structural_type = _clean_text(
        hierarchy.get("Structural Type")
        or getter("Structural Type", None)
        or (path_parts[3] if len(path_parts) >= 4 else "")
    )
    corrected_class = _clean_text(
        getter("Corrected Chemical Class", None)
        or getter("Chemical Class", None)
        or functional_group
    )
    if not parent_class and corrected_class:
        parent_class = corrected_class

    return {
        "parent_class": parent_class,
        "functional_group": functional_group,
        "subclass": subclass,
        "structural_type": structural_type,
        "corrected_class": corrected_class,
        "hierarchy_source": "taxonomy_hierarchy" if hierarchy else ("taxonomy_path" if path_parts else "class_columns"),
    }


def class_distance(left: dict[str, str], right: dict[str, str], rules: dict[str, Any] | None = None) -> float:
    """Return a review-only class distance in the 0-1 interval."""
    rules = rules or load_class_layer_rules()
    weights = rules["hierarchy_distance_weights"]

    left_known = any(left.get(key) for key in ("parent_class", "functional_group", "subclass", "structural_type", "corrected_class"))
    right_known = any(right.get(key) for key in ("parent_class", "functional_group", "subclass", "structural_type", "corrected_class"))
    if not left_known or not right_known:
        return float(weights["missing_class_information"])

    if left.get("structural_type") and left["structural_type"] == right.get("structural_type"):
        return float(weights["same_structural_type"])
    if left.get("subclass") and left["subclass"] == right.get("subclass"):
        return float(weights["same_subclass"])
    if left.get("functional_group") and left["functional_group"] == right.get("functional_group"):
        return float(weights["same_functional_group"])
    if left.get("parent_class") and left["parent_class"] == right.get("parent_class"):
        return float(weights["same_parent_class"])
    if left.get("corrected_class") and left["corrected_class"] == right.get("corrected_class"):
        return float(weights["same_corrected_class"])
    return float(weights["different_family"])


def build_class_distance_matrix(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> tuple[np.ndarray, list[dict[str, str]]]:
    """Build a pairwise chemical-class distance matrix and the extracted signatures."""
    rules = rules or load_class_layer_rules()
    signatures = [extract_class_signature(row) for _, row in df.iterrows()]
    n = len(signatures)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i):
            value = class_distance(signatures[i], signatures[j], rules)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix, signatures


def annotate_cluster_class_evidence(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    """Add cluster-level chemical class cohesion evidence without altering the cluster logic."""
    result = df.copy()
    rules = rules or load_class_layer_rules()
    labels = rules.get("cohesion_status_labels", {})
    result["Cluster_Class_Key"] = ""
    result["Cluster_Class_Cohesion"] = labels.get("not_assessable", "Not assessable")
    result["Cluster_Class_Majority"] = ""
    result["Cluster_Class_Reasons"] = ""

    if "Cluster ID" not in result.columns:
        return result

    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            continue

        signatures = {idx: extract_class_signature(result.loc[idx]) for idx in indices}
        class_keys = {
            idx: " | ".join(
                [
                    signatures[idx].get("parent_class", ""),
                    signatures[idx].get("functional_group", ""),
                    signatures[idx].get("subclass", ""),
                    signatures[idx].get("structural_type", ""),
                ]
            ).strip(" |")
            for idx in indices
        }
        for idx in indices:
            result.at[idx, "Cluster_Class_Key"] = class_keys[idx]

        non_empty = [key for key in class_keys.values() if key.strip()]
        if not non_empty:
            result.loc[indices, "Cluster_Class_Cohesion"] = labels.get("not_assessable", "Not assessable")
            result.loc[indices, "Cluster_Class_Reasons"] = "No class hierarchy available"
            continue

        counts = {}
        for key in non_empty:
            counts[key] = counts.get(key, 0) + 1
        majority_key = max(counts, key=counts.get)
        majority_count = counts[majority_key]
        total = len(non_empty)
        cohesion_ratio = majority_count / total if total else 0.0
        if cohesion_ratio >= 1.0:
            cohesion = labels.get("perfect", "Chemically coherent")
        elif cohesion_ratio >= rules.get("cluster_cohesion_threshold", 0.65):
            cohesion = labels.get("strong", "Chemically aligned")
        else:
            cohesion = labels.get("mixed", "Chemically mixed")

        if total == 1:
            cohesion = labels.get("not_assessable", "Not assessable")
            reason = "Singleton cluster; class cohesion not assessable"
        elif cohesion == labels.get("mixed", "Chemically mixed"):
            reason = f"Multiple class signatures present; majority signature {majority_key}"
        elif cohesion == labels.get("strong", "Chemically aligned"):
            reason = f"Majority class signature {majority_key} covers {cohesion_ratio:.0%} of members"
        else:
            reason = f"All members share class signature {majority_key}"

        result.loc[indices, "Cluster_Class_Cohesion"] = cohesion
        result.loc[indices, "Cluster_Class_Majority"] = majority_key
        result.loc[indices, "Cluster_Class_Reasons"] = reason

    return result
