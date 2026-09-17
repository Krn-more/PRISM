"""Review-only toxicophore/alert soft layer for clustering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_TOXICOPHORE_LAYER_RULES_PATH = Path(__file__).resolve().parent / "config" / "clustering_toxicophore_layer_rules.json"
_MISSING_TEXT = {"", "none", "n/a", "not assessable", "unknown", "not implemented"}


def _rule_candidates(rules_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if rules_path:
        candidates.append(Path(rules_path))
    candidates.append(DEFAULT_TOXICOPHORE_LAYER_RULES_PATH)
    candidates.append(Path.cwd() / "config" / "clustering_toxicophore_layer_rules.json")
    return candidates


@lru_cache(maxsize=4)
def load_toxicophore_layer_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    for candidate in _rule_candidates(rules_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                rules = json.load(handle)
            if rules.get("mode") != "review_only":
                raise ValueError("Toxicophore layer rules must remain review_only until separately approved.")
            if not str(rules.get("rule_version", "")).strip():
                raise ValueError("Toxicophore layer rules must define rule_version.")
            if not rules.get("toxicophore_columns"):
                raise ValueError("Toxicophore layer rules must define toxicophore_columns.")
            return rules
    raise FileNotFoundError("clustering_toxicophore_layer_rules.json not found")


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.casefold() in _MISSING_TEXT:
        return ""
    return text


def _tokenize_set(value: object) -> set[str]:
    text = _clean_text(value)
    if not text:
        return set()
    return {part.strip() for part in text.split(";") if part.strip() and part.strip().casefold() not in _MISSING_TEXT}


def extract_toxicophore_signature(row: pd.Series | dict[str, Any]) -> dict[str, Any]:
    getter = row.get
    profile = _tokenize_set(getter("Toxicophore_Profile", ""))
    alerts = _tokenize_set(getter("All_Structural_Alerts", ""))
    first_alert = _clean_text(getter("Alerts", "")) or "None"
    return {
        "profile": profile,
        "alerts": alerts,
        "first_alert": first_alert,
    }


def _set_distance(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 0.0
    if not left or not right:
        return 0.55
    union = left | right
    intersection = left & right
    if not union:
        return 0.0
    return 1.0 - (len(intersection) / len(union))


def _first_alert_distance(left: str, right: str) -> float:
    if left == right and left not in {"", "None", "N/A", "Unknown"}:
        return 0.0
    if left in {"", "None", "N/A", "Unknown"} and right in {"", "None", "N/A", "Unknown"}:
        return 0.0
    if left in {"", "None", "N/A", "Unknown"} or right in {"", "None", "N/A", "Unknown"}:
        return 0.45
    return 1.0


def toxicophore_distance(left: dict[str, Any], right: dict[str, Any], rules: dict[str, Any] | None = None) -> float:
    rules = rules or load_toxicophore_layer_rules()
    profile_distance = _set_distance(left.get("profile", set()), right.get("profile", set()))
    alerts_distance = _set_distance(left.get("alerts", set()), right.get("alerts", set()))
    first_alert_distance = _first_alert_distance(left.get("first_alert", "None"), right.get("first_alert", "None"))
    return float(np.clip(
        float(rules["profile_weight"]) * profile_distance
        + float(rules["all_alerts_weight"]) * alerts_distance
        + float(rules["first_alert_weight"]) * first_alert_distance,
        0.0,
        1.0,
    ))


def build_toxicophore_distance_matrix(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> tuple[np.ndarray, list[dict[str, Any]]]:
    rules = rules or load_toxicophore_layer_rules()
    signatures = [extract_toxicophore_signature(row) for _, row in df.iterrows()]
    n = len(signatures)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i):
            value = toxicophore_distance(signatures[i], signatures[j], rules)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix, signatures


def annotate_toxicophore_layer_evidence(df: pd.DataFrame, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    result = df.copy()
    rules = rules or load_toxicophore_layer_rules()
    labels = rules.get("status_labels", {})
    result["Toxicophore_Layer_Status"] = "Not assessable"
    result["Toxicophore_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
    result["Toxicophore_Layer_Key"] = ""
    result["Toxicophore_Layer_Reasons"] = ""
    result["Toxicophore_Layer_Unique_Profiles"] = np.nan
    result["Toxicophore_Layer_Unique_Alerts"] = np.nan

    if "Cluster ID" not in result.columns:
        return result

    matrix, signatures = build_toxicophore_distance_matrix(result, rules)
    result["Toxicophore_Layer_Key"] = [
        f"profile={'; '.join(sorted(sig.get('profile', set()))) or 'None'} | alerts={'; '.join(sorted(sig.get('alerts', set()))) or 'None'} | first={sig.get('first_alert', 'None')}"
        for sig in signatures
    ]

    perfect_threshold = 0.06
    aligned_threshold = 0.22

    for cluster_id, group in result.groupby("Cluster ID", dropna=False, sort=False):
        indices = list(group.index)
        if pd.isna(cluster_id):
            continue

        if len(indices) == 1:
            result.loc[indices, "Toxicophore_Layer_Status"] = "Not assessable"
            result.loc[indices, "Toxicophore_Layer_Cohesion"] = labels.get("not_assessable", "Not assessable")
            result.loc[indices, "Toxicophore_Layer_Reasons"] = "Singleton cluster; toxicophore cohesion not assessable"
            continue

        cluster_matrix = matrix[np.ix_(indices, indices)]
        pairwise = cluster_matrix[np.tril_indices(len(indices), k=-1)]
        pairwise = pairwise[np.isfinite(pairwise)]
        mean_pairwise = float(np.mean(pairwise)) if pairwise.size else 0.0

        profile_sets = [signatures[idx].get("profile", set()) for idx in indices]
        alert_sets = [signatures[idx].get("alerts", set()) for idx in indices]
        unique_profiles = sorted({"; ".join(sorted(profile)) or "None" for profile in profile_sets})
        unique_alerts = sorted({"; ".join(sorted(alerts)) or "None" for alerts in alert_sets})

        profile_counts = {}
        for profile in unique_profiles:
            profile_counts[profile] = sum(1 for item in profile_sets if ("; ".join(sorted(item)) or "None") == profile)
        majority_profile = max(profile_counts, key=profile_counts.get) if profile_counts else "None"
        majority_ratio = profile_counts.get(majority_profile, 0) / len(indices) if indices else 0.0

        if mean_pairwise <= perfect_threshold:
            cohesion = labels.get("perfect", "Toxicophore coherent")
        elif mean_pairwise <= aligned_threshold:
            cohesion = labels.get("strong", "Toxicophore aligned")
        else:
            cohesion = labels.get("mixed", "Toxicophore mixed")

        result.loc[indices, "Toxicophore_Layer_Status"] = "Calculated"
        result.loc[indices, "Toxicophore_Layer_Cohesion"] = cohesion
        result.loc[indices, "Toxicophore_Layer_Unique_Profiles"] = len(unique_profiles)
        result.loc[indices, "Toxicophore_Layer_Unique_Alerts"] = len(unique_alerts)
        result.loc[indices, "Toxicophore_Layer_Reasons"] = (
            f"Mean pairwise toxicophore distance {mean_pairwise:.3f}; "
            f"majority profile {majority_profile} ({majority_ratio:.0%}); "
            f"unique profiles {len(unique_profiles)}; unique alerts {len(unique_alerts)}; "
            f"keys: {', '.join(sorted(set(result.loc[indices, 'Toxicophore_Layer_Key'].tolist())))}"
        )

    return result
