"""Review-only ionisation soft layer for clustering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_IONISATION_LAYER_RULES_PATH = Path(__file__).resolve().parent / "config" / "clustering_ionisation_layer_rules.json"
_MISSING_TEXT = {"", "none", "n/a", "not assessable", "unknown", "not implemented"}


def _rule_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_IONISATION_LAYER_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "clustering_ionisation_layer_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_ionisation_layer_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _rule_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                rules = json.load(handle)
            if rules.get("mode") != "review_only":
                raise ValueError("Ionisation layer rules must remain review_only until separately approved.")
            if not str(rules.get("rule_version", "")).strip():
                raise ValueError("Ionisation layer rules must define rule_version.")
            if not rules.get("ionisation_columns"):
                raise ValueError("Ionisation layer rules must define ionisation_columns.")
            return rules
    raise FileNotFoundError("clustering_ionisation_layer_rules.json not found")


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.casefold() in _MISSING_TEXT:
        return ""
    return text


def _parse_charge(value: object) -> int | None:
    text = _clean_text(value)
    if not text:
        return None
    try:
        return int(float(text))
    except Exception:
        return None


def extract_ionisation_signature(row: pd.Series | dict[str, Any]) -> dict[str, Any]:
    getter = row.get
    indicator = _clean_text(getter("Ionisation_Indicator", "")) or "Unknown"
    charge = _parse_charge(getter("Formal_Charge", None))
    if charge is None:
        charge = 0 if indicator == "Neutral" else None
    charge_band = "Unknown"
    if charge is not None:
        if charge > 0:
            charge_band = "Positive"
        elif charge < 0:
            charge_band = "Negative"
        else:
            charge_band = "Neutral"
    return {
        "indicator": indicator,
        "formal_charge": charge,
        "charge_band": charge_band,
    }


def _indicator_distance(left: str, right: str, rules: dict[str, Any]) -> float:
    table = rules["indicator_distance"]
    if left == right and left != "Unknown":
        return float(table["same_known"])
    if left == right == "Unknown":
        return float(table["same_unknown"])
    if "Unknown" in {left, right}:
        return float(table["one_unknown"])
    charged_labels = {"Permanently charged", "Likely acidic", "Likely basic", "Amphoteric"}
    if left in charged_labels and right in charged_labels:
        if {left, right} in ({"Likely acidic", "Likely basic"}, {"Likely basic", "Likely acidic"}):
            return float(table["acid_base_opposition"])
        if "Amphoteric" in {left, right}:
            return 0.55
        if "Permanently charged" in {left, right} and (("Neutral" in {left, right}) or ("Unknown" in {left, right})):
            return float(table["charged_vs_unknown"])
        if "Permanently charged" in {left, right}:
            return 0.80
        return float(table["different_known"])
    if {left, right} == {"Neutral", "Likely acidic"} or {left, right} == {"Neutral", "Likely basic"}:
        return float(table["charged_vs_neutral"])
    if "Neutral" in {left, right} and "Permanently charged" in {left, right}:
        return float(table["charged_vs_neutral"])
    return float(table["different_known"])


def ionisation_distance(left: dict[str, Any], right: dict[str, Any], rules: dict[str, Any] | None = None) -> float:
    rules = rules or load_ionisation_layer_rules()
    indicator_distance = _indicator_distance(left.get("indicator", "Unknown"), right.get("indicator", "Unknown"), rules)

    max_charge_difference = max(int(rules.get("max_charge_difference", 4)), 1)
    left_charge = left.get("formal_charge")
    right_charge = right.get("formal_charge")
    if left_charge is None and right_charge is None:
        charge_distance = 0.25
    elif left_charge is None or right_charge is None:
        charge_distance = 0.50
    else:
        charge_distance = min(abs(int(left_charge) - int(right_charge)) / max_charge_difference, 1.0)

    return float(np.clip(
        float(rules["indicator_weight"]) * indicator_distance
        + float(rules["formal_charge_weight"]) * charge_distance,
        0.0,
        1.0,
    ))


def build_ionisation_distance_matrix(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> tuple[np.ndarray, list[dict[str, Any]]]:
    rules = rules or load_ionisation_layer_rules()
    signatures = [extract_ionisation_signature(row) for _, row in df.iterrows()]
    n = len(signatures)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i):
            value = ionisation_distance(signatures[i], signatures[j], rules)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix, signatures


def annotate_ionisation_layer_evidence(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    result = df.copy()
    rules = rules or load_ionisation_layer_rules()
    labels = rules.get("status_labels", {})
    result["Ionisation_Layer_Status"] = "Not assessable"
    result["Ionisation_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
    result["Ionisation_Layer_Key"] = ""
    result["Ionisation_Layer_Reasons"] = ""
    result["Ionisation_Layer_Charge_Span"] = np.nan

    if "Cluster ID" not in result.columns:
        return result

    matrix, signatures = build_ionisation_distance_matrix(result, rules)
    result["Ionisation_Layer_Key"] = [
        f"{sig.get('indicator', 'Unknown')} | charge {sig.get('formal_charge', 'NA')} | band {sig.get('charge_band', 'Unknown')}"
        for sig in signatures
    ]

    perfect_threshold = 0.08
    aligned_threshold = 0.28

    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            continue

        if len(indices) == 1:
            result.loc[indices, "Ionisation_Layer_Status"] = "Not assessable"
            result.loc[indices, "Ionisation_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
            result.loc[indices, "Ionisation_Layer_Reasons"] = "Singleton cluster; ionisation cohesion not assessable"
            continue

        cluster_matrix = matrix[np.ix_(indices, indices)]
        pairwise = cluster_matrix[np.tril_indices(len(indices), k=-1)]
        pairwise = pairwise[np.isfinite(pairwise)]
        mean_pairwise = float(np.mean(pairwise)) if pairwise.size else 0.0

        charges = [signatures[idx].get("formal_charge") for idx in indices if signatures[idx].get("formal_charge") is not None]
        if charges:
            charge_span = max(charges) - min(charges)
        else:
            charge_span = np.nan

        indicator_counts = {}
        for idx in indices:
            key = signatures[idx].get("indicator", "Unknown") or "Unknown"
            indicator_counts[key] = indicator_counts.get(key, 0) + 1
        majority_indicator = max(indicator_counts, key=indicator_counts.get)
        majority_ratio = indicator_counts[majority_indicator] / len(indices)

        if mean_pairwise <= perfect_threshold:
            cohesion = labels.get("perfect", "Ionisation coherent")
        elif mean_pairwise <= aligned_threshold:
            cohesion = labels.get("strong", "Ionisation aligned")
        else:
            cohesion = labels.get("mixed", "Ionisation mixed")

        result.loc[indices, "Ionisation_Layer_Status"] = "Calculated"
        result.loc[indices, "Ionisation_Layer_Cohesion"] = cohesion
        result.loc[indices, "Ionisation_Layer_Charge_Span"] = charge_span
        result.loc[indices, "Ionisation_Layer_Reasons"] = (
            f"Mean pairwise ionisation distance {mean_pairwise:.3f}; "
            f"majority indicator {majority_indicator} ({majority_ratio:.0%}); "
            f"charge span {charge_span if pd.notna(charge_span) else 'NA'}; "
            f"ionisation keys: {', '.join(sorted(set(result.loc[indices, 'Ionisation_Layer_Key'].tolist())))}"
        )

    return result
