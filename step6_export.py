from __future__ import annotations

from typing import Any, Mapping


STEP6_EXPORT_FIELDS = [
    "Structural Evidence Status",
    "Structural Evidence Version",
    "Toxicophore Profile",
    "All Structural Alerts",
    "Alerts",
    "Scaffold",
]


def build_step6_export_fields(
    structural_evidence: Mapping[str, Any] | None,
    scaffold: str = "",
) -> dict[str, Any]:
    """Return the Step 6 workbook mirror fields from the shared evidence source."""
    evidence = dict(structural_evidence or {})
    return {
        "Structural Evidence Status": evidence.get("Structural_Evidence_Status", ""),
        "Structural Evidence Version": evidence.get("Structural_Evidence_Version", ""),
        "Toxicophore Profile": evidence.get("Toxicophore_Profile", "None"),
        "All Structural Alerts": evidence.get("All_Structural_Alerts", "None"),
        "Alerts": evidence.get("Alerts", "None"),
        "Scaffold": scaffold or "Not available",
    }
