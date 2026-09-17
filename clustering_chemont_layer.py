"""Review-only ChemOnt/taxonomy soft layer for clustering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from clustering_class_layer import extract_class_signature


DEFAULT_CHEMONT_LAYER_RULES_PATH = Path(__file__).resolve().parent / "config" / "clustering_chemont_layer_rules.json"
_MISSING_TEXT = {"", "none", "n/a", "not assessable", "unknown", "not implemented", "unclassified"}


def _rule_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_CHEMONT_LAYER_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "clustering_chemont_layer_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_chemont_layer_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _rule_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                rules = json.load(handle)
            if rules.get("mode") != "review_only":
                raise ValueError("ChemOnt layer rules must remain review_only until separately approved.")
            if not str(rules.get("rule_version", "")).strip():
                raise ValueError("ChemOnt layer rules must define rule_version.")
            if not rules.get("hierarchy_columns"):
                raise ValueError("ChemOnt layer rules must define hierarchy_columns.")
            return rules
    raise FileNotFoundError("clustering_chemont_layer_rules.json not found")


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.casefold() in _MISSING_TEXT:
        return ""
    return text


def _split_tokens(value: object) -> set[str]:
    text = _clean_text(value)
    if not text:
        return set()
    return {part.strip() for part in text.split(";") if part.strip() and part.strip().casefold() not in _MISSING_TEXT}


def _parse_taxonomy_hierarchy(value: object) -> dict[str, str]:
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


def _parse_taxonomy_path(value: object) -> str:
    text = _clean_text(value)
    return text


def extract_chemont_signature(row: pd.Series | dict[str, Any]) -> dict[str, Any]:
    getter = row.get
    chemont_values = {
        "ChemOnt_Kingdom": _clean_text(getter("ChemOnt_Kingdom", "")),
        "ChemOnt_Superclass": _clean_text(getter("ChemOnt_Superclass", "")),
        "ChemOnt_Class": _clean_text(getter("ChemOnt_Class", "")),
        "ChemOnt_Subclass": _clean_text(getter("ChemOnt_Subclass", "")),
        "ChemOnt_Direct_Parent": _clean_text(getter("ChemOnt_Direct_Parent", "")),
        "ChemOnt_Molecular_Framework": _clean_text(getter("ChemOnt_Molecular_Framework", "")),
        "ChemOnt_Substituents": _split_tokens(getter("ChemOnt_Substituents", "")),
    }
    taxonomy_hierarchy = _parse_taxonomy_hierarchy(getter("Taxonomy Hierarchy", ""))
    taxonomy_path = _parse_taxonomy_path(getter("Taxonomy Path", ""))
    corrected_class = _clean_text(getter("Corrected Chemical Class", "")) or _clean_text(getter("Chemical Class", ""))
    source = "chemont" if any(chemont_values[key] for key in (
        "ChemOnt_Kingdom", "ChemOnt_Superclass", "ChemOnt_Class", "ChemOnt_Subclass",
        "ChemOnt_Direct_Parent", "ChemOnt_Molecular_Framework"
    )) else ("taxonomy" if taxonomy_hierarchy or taxonomy_path else ("class" if corrected_class else "unassessable"))
    class_signature = extract_class_signature(row)
    return {
        "chemont": chemont_values,
        "taxonomy_hierarchy": taxonomy_hierarchy,
        "taxonomy_path": taxonomy_path,
        "class_signature": class_signature,
        "corrected_class": corrected_class,
        "source": source,
    }


def _same_or_mismatch(left: str, right: str, missing_penalty: float = 0.5) -> float:
    left_known = bool(left)
    right_known = bool(right)
    if left_known and right_known:
        return 0.0 if left == right else 1.0
    if not left_known and not right_known:
        return 0.0
    return missing_penalty


def _set_distance(left: set[str], right: set[str], rules: dict[str, Any]) -> float:
    params = rules["substituent_distance"]
    if not left and not right:
        return float(params["both_unknown"])
    if not left or not right:
        return float(params["one_unknown"])
    if left == right:
        return float(params["same"])
    if left.issubset(right) or right.issubset(left):
        return float(params["subset"])
    overlap = len(left & right)
    if overlap:
        return float(params["partial_overlap"])
    return float(params["different"])


def chemont_distance(left: dict[str, Any], right: dict[str, Any], rules: dict[str, Any] | None = None) -> float:
    rules = rules or load_chemont_layer_rules()
    weights = rules["level_weights"]

    distance = 0.0
    distance += weights["ChemOnt_Kingdom"] * _same_or_mismatch(left["chemont"].get("ChemOnt_Kingdom", ""), right["chemont"].get("ChemOnt_Kingdom", ""))
    distance += weights["ChemOnt_Superclass"] * _same_or_mismatch(left["chemont"].get("ChemOnt_Superclass", ""), right["chemont"].get("ChemOnt_Superclass", ""))
    distance += weights["ChemOnt_Class"] * _same_or_mismatch(left["chemont"].get("ChemOnt_Class", ""), right["chemont"].get("ChemOnt_Class", ""))
    distance += weights["ChemOnt_Subclass"] * _same_or_mismatch(left["chemont"].get("ChemOnt_Subclass", ""), right["chemont"].get("ChemOnt_Subclass", ""))
    distance += weights["ChemOnt_Direct_Parent"] * _same_or_mismatch(left["chemont"].get("ChemOnt_Direct_Parent", ""), right["chemont"].get("ChemOnt_Direct_Parent", ""))
    distance += weights["ChemOnt_Molecular_Framework"] * _same_or_mismatch(left["chemont"].get("ChemOnt_Molecular_Framework", ""), right["chemont"].get("ChemOnt_Molecular_Framework", ""))
    distance += weights["ChemOnt_Substituents"] * _set_distance(left["chemont"].get("ChemOnt_Substituents", set()), right["chemont"].get("ChemOnt_Substituents", set()), rules)

    # Taxonomy fallback, used when ChemOnt is absent or sparse.
    taxonomy_left = left.get("taxonomy_hierarchy", {})
    taxonomy_right = right.get("taxonomy_hierarchy", {})
    if taxonomy_left or taxonomy_right:
        distance += weights["Taxonomy_Hierarchy"] * (
            0.0 if taxonomy_left == taxonomy_right and taxonomy_left else 1.0 if taxonomy_left and taxonomy_right else 0.5
        )
    else:
        distance += weights["Taxonomy_Hierarchy"] * 0.5

    distance += weights["Taxonomy_Path"] * (
        0.0 if left.get("taxonomy_path") and left.get("taxonomy_path") == right.get("taxonomy_path") else 1.0 if left.get("taxonomy_path") and right.get("taxonomy_path") else 0.5
    )

    # When ChemOnt is absent, gently fall back to the existing class hierarchy.
    if left.get("source") != "chemont" or right.get("source") != "chemont":
        class_left = left.get("class_signature", {})
        class_right = right.get("class_signature", {})
        distance += 0.06 * (
            _same_or_mismatch(class_left.get("parent_class", ""), class_right.get("parent_class", "")) * 0.25
            + _same_or_mismatch(class_left.get("functional_group", ""), class_right.get("functional_group", "")) * 0.25
            + _same_or_mismatch(class_left.get("subclass", ""), class_right.get("subclass", "")) * 0.25
            + _same_or_mismatch(class_left.get("structural_type", ""), class_right.get("structural_type", "")) * 0.25
        )

    return float(np.clip(distance, 0.0, 1.0))


def build_chemont_distance_matrix(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> tuple[np.ndarray, list[dict[str, Any]]]:
    rules = rules or load_chemont_layer_rules()
    signatures = [extract_chemont_signature(row) for _, row in df.iterrows()]
    n = len(signatures)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i):
            value = chemont_distance(signatures[i], signatures[j], rules)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix, signatures


def annotate_chemont_layer_evidence(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    result = df.copy()
    rules = rules or load_chemont_layer_rules()
    labels = rules.get("status_labels", {})
    result["ChemOnt_Layer_Status"] = "Not assessable"
    result["ChemOnt_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
    result["ChemOnt_Layer_Key"] = ""
    result["ChemOnt_Layer_Reasons"] = ""
    result["ChemOnt_Layer_Source"] = ""
    result["ChemOnt_Layer_Unique_Classes"] = np.nan
    result["ChemOnt_Layer_Unique_Frameworks"] = np.nan

    if "Cluster ID" not in result.columns:
        return result

    matrix, signatures = build_chemont_distance_matrix(result, rules)
    result["ChemOnt_Layer_Source"] = [sig.get("source", "unassessable") for sig in signatures]
    result["ChemOnt_Layer_Key"] = [
        " | ".join([
            sig["chemont"].get("ChemOnt_Kingdom", "") or "No kingdom",
            sig["chemont"].get("ChemOnt_Superclass", "") or "No superclass",
            sig["chemont"].get("ChemOnt_Class", "") or "No class",
            sig["chemont"].get("ChemOnt_Subclass", "") or "No subclass",
            sig["chemont"].get("ChemOnt_Direct_Parent", "") or "No direct parent",
            sig["chemont"].get("ChemOnt_Molecular_Framework", "") or "No framework",
        ])
        for sig in signatures
    ]

    perfect_threshold = float(rules["cohesion_thresholds"]["perfect_mean_distance"])
    aligned_threshold = float(rules["cohesion_thresholds"]["aligned_mean_distance"])

    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            continue

        if len(indices) == 1:
            result.loc[indices, "ChemOnt_Layer_Status"] = "Not assessable"
            result.loc[indices, "ChemOnt_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
            result.loc[indices, "ChemOnt_Layer_Reasons"] = "Singleton cluster; ChemOnt cohesion not assessable"
            continue

        cluster_matrix = matrix[np.ix_(indices, indices)]
        pairwise = cluster_matrix[np.tril_indices(len(indices), k=-1)]
        pairwise = pairwise[np.isfinite(pairwise)]
        mean_pairwise = float(np.mean(pairwise)) if pairwise.size else 0.0

        unique_classes = sorted({sig["chemont"].get("ChemOnt_Class", "") or sig["class_signature"].get("corrected_class", "") or "Unknown" for sig in [signatures[idx] for idx in indices]})
        unique_frameworks = sorted({sig["chemont"].get("ChemOnt_Molecular_Framework", "") or "Unknown" for sig in [signatures[idx] for idx in indices]})

        class_counts = {}
        for sig in [signatures[idx] for idx in indices]:
            key = sig["chemont"].get("ChemOnt_Class", "") or sig["class_signature"].get("corrected_class", "") or "Unknown"
            class_counts[key] = class_counts.get(key, 0) + 1
        majority_class = max(class_counts, key=class_counts.get)
        majority_ratio = class_counts[majority_class] / len(indices)

        if mean_pairwise <= perfect_threshold:
            cohesion = labels.get("perfect", "Taxonomically coherent")
        elif mean_pairwise <= aligned_threshold:
            cohesion = labels.get("strong", "Taxonomically aligned")
        else:
            cohesion = labels.get("mixed", "Taxonomically mixed")

        result.loc[indices, "ChemOnt_Layer_Status"] = "Calculated"
        result.loc[indices, "ChemOnt_Layer_Cohesion"] = cohesion
        result.loc[indices, "ChemOnt_Layer_Unique_Classes"] = len(unique_classes)
        result.loc[indices, "ChemOnt_Layer_Unique_Frameworks"] = len(unique_frameworks)
        result.loc[indices, "ChemOnt_Layer_Reasons"] = (
            f"Mean pairwise taxonomy distance {mean_pairwise:.3f}; "
            f"majority class {majority_class} ({majority_ratio:.0%}); "
            f"unique classes {len(unique_classes)}; unique frameworks {len(unique_frameworks)}; "
            f"sources: {', '.join(sorted(set(result.loc[indices, 'ChemOnt_Layer_Source'].tolist())))}"
        )

    return result
