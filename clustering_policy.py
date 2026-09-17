"""Frozen clustering policy metadata.

This module intentionally stays lightweight so the policy can be validated
without importing the full RDKit-based clustering stack.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_CLUSTERING_POLICY_PATH = Path(__file__).resolve().parent / "config" / "stage1_clustering_policy.json"


def _policy_candidates(policy_path: str | Path | None = None) -> list[Path]:
    candidates: list[Path] = []
    if policy_path:
        candidates.append(Path(policy_path))
    candidates.append(DEFAULT_CLUSTERING_POLICY_PATH)
    candidates.append(Path.cwd() / "config" / "stage1_clustering_policy.json")
    return candidates


@lru_cache(maxsize=4)
def load_clustering_policy(policy_path: str | Path | None = None) -> dict[str, Any]:
    """Load the frozen clustering policy used to document the current baseline."""
    for candidate in _policy_candidates(policy_path):
        if candidate.exists():
            with open(candidate, "r", encoding="utf-8") as handle:
                policy = json.load(handle)
            if policy.get("mode") != "review_only":
                raise ValueError("Clustering policy must remain review_only until separately approved.")
            if "rule_version" not in policy or not str(policy["rule_version"]).strip():
                raise ValueError("Clustering policy must define a rule_version.")
            components = policy.get("distance_components", [])
            if not isinstance(components, list) or not components:
                raise ValueError("Clustering policy must define distance_components.")
            total_weight = 0.0
            for component in components:
                weight = component.get("weight")
                if not isinstance(weight, (int, float)) or weight < 0:
                    raise ValueError("Clustering policy weights must be non-negative numbers.")
                total_weight += float(weight)
            if not np.isclose(total_weight, 1.0, atol=1e-6):
                raise ValueError("Clustering policy distance weights must sum to 1.0.")
            return policy
    raise FileNotFoundError("stage1_clustering_policy.json not found")
