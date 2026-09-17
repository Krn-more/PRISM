"""Review-only topology/ring soft layer for clustering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_TOPOLOGY_LAYER_RULES_PATH = Path(__file__).resolve().parent / "config" / "clustering_topology_layer_rules.json"
_MISSING_TEXT = {"", "none", "n/a", "not assessable", "unknown", "not implemented", "unclassified"}


def _rule_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_TOPOLOGY_LAYER_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "clustering_topology_layer_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_topology_layer_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _rule_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                rules = json.load(handle)
            if rules.get("mode") != "review_only":
                raise ValueError("Topology layer rules must remain review_only until separately approved.")
            if not str(rules.get("rule_version", "")).strip():
                raise ValueError("Topology layer rules must define rule_version.")
            if not rules.get("topology_columns"):
                raise ValueError("Topology layer rules must define topology_columns.")
            return rules
    raise FileNotFoundError("clustering_topology_layer_rules.json not found")


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.casefold() in _MISSING_TEXT:
        return ""
    return text


def _split_modifiers(value: object) -> set[str]:
    text = _clean_text(value)
    if not text:
        return set()
    return {part.strip() for part in text.split(";") if part.strip() and part.strip().casefold() not in _MISSING_TEXT}


def _parse_int(value: object) -> int | None:
    text = _clean_text(value)
    if not text:
        return None
    try:
        return int(float(text))
    except Exception:
        return None


def extract_topology_signature(row: pd.Series | dict[str, Any]) -> dict[str, Any]:
    getter = row.get
    topology_class = _clean_text(getter("Topology Class", ""))
    modifiers = _split_modifiers(getter("Topology Modifiers", ""))
    ring_system = _clean_text(getter("Ring System", ""))
    ring_profile = _clean_text(getter("Ring Profile", ""))
    alcohol_context = _clean_text(getter("Alcohol Context", ""))
    ring_count = _parse_int(getter("Ring Count", None))
    aromatic_ring_count = _parse_int(getter("Aromatic Ring Count", None))

    if not topology_class and ring_count is not None:
        if ring_count == 0:
            topology_class = "Acyclic"
            ring_system = ring_system or "Open-chain"
        elif ring_count == 1:
            topology_class = "Monocyclic"
        elif ring_count == 2:
            topology_class = "Bicyclic"
        else:
            topology_class = "Polycyclic"

    if not ring_profile and ring_count is not None:
        aromatic_count = aromatic_ring_count if aromatic_ring_count is not None else 0
        ring_profile = f"{ring_count} ring(s), {aromatic_count} aromatic ring(s)"

    if not ring_system and topology_class:
        ring_system = topology_class

    if not alcohol_context and topology_class:
        alcohol_context = "Not applicable"

    return {
        "topology_class": topology_class,
        "modifiers": modifiers,
        "ring_system": ring_system,
        "ring_profile": ring_profile,
        "alcohol_context": alcohol_context,
        "ring_count": ring_count,
        "aromatic_ring_count": aromatic_ring_count,
    }


def topology_distance(left: dict[str, Any], right: dict[str, Any], rules: dict[str, Any] | None = None) -> float:
    rules = rules or load_topology_layer_rules()
    cat_weights = rules["categorical_distance_weights"]
    num_weights = rules["numeric_distance_weights"]

    left_known = any(
        value not in (None, "", set())
        for value in (
            left.get("topology_class"),
            left.get("ring_system"),
            left.get("ring_profile"),
            left.get("ring_count"),
            left.get("aromatic_ring_count"),
            left.get("modifiers"),
        )
    )
    right_known = any(
        value not in (None, "", set())
        for value in (
            right.get("topology_class"),
            right.get("ring_system"),
            right.get("ring_profile"),
            right.get("ring_count"),
            right.get("aromatic_ring_count"),
            right.get("modifiers"),
        )
    )
    if not left_known or not right_known:
        return 0.5

    distance = 0.0

    if left.get("topology_class") and right.get("topology_class"):
        distance += cat_weights["Topology Class"] * (0.0 if left["topology_class"] == right["topology_class"] else 1.0)
    else:
        distance += 0.5 * cat_weights["Topology Class"]

    left_mods = left.get("modifiers", set()) or set()
    right_mods = right.get("modifiers", set()) or set()
    if left_mods or right_mods:
        union = len(left_mods | right_mods)
        inter = len(left_mods & right_mods)
        mod_distance = 1.0 if union == 0 else 1.0 - (inter / union)
        distance += cat_weights["Topology Modifiers"] * mod_distance
    else:
        distance += 0.5 * cat_weights["Topology Modifiers"]

    if left.get("ring_system") and right.get("ring_system"):
        distance += cat_weights["Ring System"] * (0.0 if left["ring_system"] == right["ring_system"] else 1.0)
    else:
        distance += 0.5 * cat_weights["Ring System"]

    if left.get("ring_profile") and right.get("ring_profile"):
        distance += cat_weights["Ring Profile"] * (0.0 if left["ring_profile"] == right["ring_profile"] else 1.0)
    else:
        distance += 0.5 * cat_weights["Ring Profile"]

    if left.get("alcohol_context") and right.get("alcohol_context"):
        distance += cat_weights["Alcohol Context"] * (0.0 if left["alcohol_context"] == right["alcohol_context"] else 1.0)
    else:
        distance += 0.5 * cat_weights["Alcohol Context"]

    ring_left = left.get("ring_count")
    ring_right = right.get("ring_count")
    if ring_left is not None and ring_right is not None:
        max_ring = max(ring_left, ring_right, 1)
        distance += num_weights["Ring Count"] * abs(ring_left - ring_right) / max_ring
    else:
        distance += 0.5 * num_weights["Ring Count"]

    aromatic_left = left.get("aromatic_ring_count")
    aromatic_right = right.get("aromatic_ring_count")
    if aromatic_left is not None and aromatic_right is not None:
        max_aromatic = max(aromatic_left, aromatic_right, 1)
        distance += num_weights["Aromatic Ring Count"] * abs(aromatic_left - aromatic_right) / max_aromatic
    else:
        distance += 0.5 * num_weights["Aromatic Ring Count"]

    return float(np.clip(distance, 0.0, 1.0))


def build_topology_distance_matrix(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> tuple[np.ndarray, list[dict[str, Any]]]:
    rules = rules or load_topology_layer_rules()
    signatures = [extract_topology_signature(row) for _, row in df.iterrows()]
    n = len(signatures)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i):
            value = topology_distance(signatures[i], signatures[j], rules)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix, signatures


def annotate_topology_layer_evidence(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    result = df.copy()
    rules = rules or load_topology_layer_rules()
    labels = rules.get("status_labels", {})
    result["Topology_Layer_Status"] = "Not assessable"
    result["Topology_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
    result["Topology_Layer_Key"] = ""
    result["Topology_Layer_Reasons"] = ""

    if "Cluster ID" not in result.columns:
        return result

    matrix, signatures = build_topology_distance_matrix(result, rules)
    result["Topology_Layer_Key"] = [
        " | ".join(
            [
                sig.get("topology_class", ""),
                ",".join(sorted(sig.get("modifiers", set()))) if sig.get("modifiers") else "",
                sig.get("ring_system", ""),
                sig.get("ring_profile", ""),
                sig.get("alcohol_context", ""),
            ]
        ).strip(" |")
        for sig in signatures
    ]

    perfect_threshold = 0.08
    aligned_threshold = 0.25

    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            continue

        if len(indices) == 1:
            result.loc[indices, "Topology_Layer_Status"] = "Not assessable"
            result.loc[indices, "Topology_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
            result.loc[indices, "Topology_Layer_Reasons"] = "Singleton cluster; topology cohesion not assessable"
            continue

        cluster_matrix = matrix[np.ix_(indices, indices)]
        pairwise = cluster_matrix[np.tril_indices(len(indices), k=-1)]
        pairwise = pairwise[np.isfinite(pairwise)]
        mean_pairwise = float(np.mean(pairwise)) if pairwise.size else 0.0

        signatures_subset = [signatures[idx] for idx in indices]
        class_counts = {}
        for sig in signatures_subset:
            key = sig.get("topology_class", "") or "Unknown"
            class_counts[key] = class_counts.get(key, 0) + 1
        majority_class = max(class_counts, key=class_counts.get)
        majority_ratio = class_counts[majority_class] / len(indices)

        if mean_pairwise <= perfect_threshold:
            cohesion = labels.get("perfect", "Topologically coherent")
        elif mean_pairwise <= aligned_threshold:
            cohesion = labels.get("strong", "Topologically aligned")
        else:
            cohesion = labels.get("mixed", "Topologically mixed")

        result.loc[indices, "Topology_Layer_Status"] = "Calculated"
        result.loc[indices, "Topology_Layer_Cohesion"] = cohesion
        result.loc[indices, "Topology_Layer_Reasons"] = (
            f"Mean pairwise topology distance {mean_pairwise:.3f}; "
            f"majority class {majority_class} ({majority_ratio:.0%}); "
            f"topology keys: {', '.join(sorted(set(result.loc[indices, 'Topology_Layer_Key'].tolist())))}"
        )

    return result
