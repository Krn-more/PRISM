"""Stage 3 deterministic, review-only uncertainty register."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_RULES_PATH = Path(__file__).with_name("config") / "stage3_uncertainty_rules.json"
LEVEL_ORDER = {"Low": 0, "Medium": 1, "High": 2}


def load_uncertainty_rules(path: str | Path | None = None) -> dict[str, Any]:
    """Load the versioned Stage 3 review-only rules."""
    with open(path or DEFAULT_RULES_PATH, encoding="utf-8") as handle:
        rules = json.load(handle)
    if rules.get("mode") != "review_only":
        raise ValueError("Stage 3 rules must remain review_only until separately approved.")
    if rules.get("categories") != ["Identity", "Structure", "Reactivity", "Property", "Biological", "Exposure", "Domain"]:
        raise ValueError("Stage 3 rules must define the seven controlled uncertainty categories.")
    for category in rules["categories"]:
        if category.casefold() not in rules:
            raise ValueError(f"Stage 3 configuration is missing {category} rules.")
    return rules


def _text(value: object) -> str:
    if isinstance(value, pd.Series):
        non_null = value.dropna()
        if non_null.empty:
            return ""
        value = non_null.iloc[0]
    elif isinstance(value, (list, tuple)) and value:
        value = value[0]
    elif isinstance(value, dict):
        value = next(iter(value.values()), "")

    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _has_value(value: object) -> bool:
    return bool(_text(value))


def _outcome(rule: dict[str, str], reason: str) -> tuple[str, str, str]:
    """Build a level/reason/action tuple from the versioned rule configuration."""
    return rule["level"], reason, rule["action"]


def _record(scope: str, scope_id: str, category: str, level: str, reason: str, required_action: str, evidence_fields: str, version: str) -> dict[str, str]:
    if level in {"Medium", "High"} and not required_action:
        raise ValueError(f"{category} {level} uncertainty requires an action.")
    return {
        "Scope": scope,
        "Scope_ID": scope_id,
        "Category": category,
        "Level": level,
        "Reason": reason,
        "Required_Action": required_action,
        "Evidence_Fields": evidence_fields,
        "Rule_Version": version,
    }


def _member_records(row: pd.Series, row_id: str, rules: dict[str, Any]) -> list[dict[str, str]]:
    version = rules["rule_version"]
    identity = _text(row.get("Identity Status"))
    identity_lower = identity.casefold()
    structure = _text(row.get("Structural_Evidence_Status"))
    alerts = _text(row.get("All_Structural_Alerts"))
    toxicophores = _text(row.get("Toxicophore_Profile"))
    property_flags = _text(row.get("Property_Outlier_Flags"))
    domain = _text(row.get("Domain_Status"))
    exposure_fields = [column for column in rules["exposure_context_columns"] if _has_value(row.get(column))]

    identity_rules = rules["identity"]
    if not identity:
        identity_record = _outcome(identity_rules["missing"], "Identity status is not supplied.")
    elif any(term.casefold() in identity_lower for term in identity_rules["unresolved_terms"]):
        identity_record = _outcome(identity_rules["unresolved"], f"Identity status is unresolved: {identity}")
    elif any(term.casefold() in identity_lower for term in identity_rules["review_terms"]):
        identity_record = _outcome(identity_rules["review"], f"Identity status requires review: {identity}")
    elif any(term.casefold() in identity_lower for term in identity_rules["resolved_terms"]):
        identity_record = _outcome(identity_rules["resolved"], f"Identity status is resolved: {identity}.")
    else:
        identity_record = _outcome(identity_rules["other"], f"Identity status is not a configured resolved status: {identity}.")

    structure_rules = rules["structure"]
    if structure != "Calculated":
        structure_record = _outcome(structure_rules["unavailable"], f"Structural evidence status is {structure or 'not supplied'}.")
    else:
        structure_record = _outcome(structure_rules["calculated"], "Local structural evidence was calculated from standardized SMILES.")

    reactivity_rules = rules["reactivity"]
    if alerts not in {"", "None", "Not assessable"} or toxicophores not in {"", "None", "Not assessable"}:
        reactivity_record = _outcome(reactivity_rules["flagged"], f"Structural alert/toxicophore review required: alerts={alerts or 'None'}; toxicophores={toxicophores or 'None'}.")
    elif structure != "Calculated":
        reactivity_record = _outcome(reactivity_rules["unavailable"], "Reactivity evidence cannot be assessed without a valid structure.")
    else:
        reactivity_record = _outcome(reactivity_rules["clear"], "No Stage 1 structural alerts or configured toxicophores were identified.")

    property_rules = rules["property"]
    if "Not assessed" in property_flags:
        property_record = _outcome(property_rules["incomplete"], f"Property outlier screening is incomplete: {property_flags}.")
    elif property_flags not in {"", "None"}:
        property_record = _outcome(property_rules["outlier"], f"Property outlier flag(s): {property_flags}.")
    else:
        property_record = _outcome(property_rules["clear"], "No Stage 2 property outlier flags were reported.")

    biological_record = _outcome(rules["biological"], rules["biological"]["reason"] + ".")
    exposure_rules = rules["exposure"]
    if exposure_fields:
        exposure_record = _outcome(exposure_rules["provided"], "Analytical concentration field(s) supplied: " + ", ".join(exposure_fields) + "; exposure context remains incomplete.")
    else:
        exposure_record = _outcome(exposure_rules["missing"], "No configured exposure-context field is supplied.")

    domain_rules = rules["domain"]
    if not domain:
        domain_record = _outcome(domain_rules["missing"], "Stage 2 domain status is not supplied.")
    elif domain in domain_rules["high_statuses"]:
        domain_record = _outcome(domain_rules["high"], f"Stage 2 domain status is {domain}.")
    elif domain in domain_rules["medium_statuses"]:
        domain_record = _outcome(domain_rules["medium"], f"Stage 2 domain status is {domain}.")
    elif domain in domain_rules["low_statuses"]:
        domain_record = _outcome(domain_rules["low"], f"Stage 2 domain status is {domain}.")
    else:
        domain_record = _outcome(domain_rules["other"], f"Stage 2 domain status is unrecognised: {domain}.")

    values = {
        "Identity": (identity_record, "Identity Status"),
        "Structure": (structure_record, "Structural_Evidence_Status; Standardized SMILES"),
        "Reactivity": (reactivity_record, "All_Structural_Alerts; Toxicophore_Profile"),
        "Property": (property_record, "Property_Outlier_Flags; MW; LogP; TPSA"),
        "Biological": (biological_record, "No external biological evidence in Stage 3"),
        "Exposure": (exposure_record, "; ".join(rules["exposure_context_columns"])),
        "Domain": (domain_record, "Domain_Status; Domain_Reasons; Nearest_Neighbour_Tanimoto"),
    }
    return [_record("Compound", row_id, category, *values[category][0], values[category][1], version) for category in rules["categories"]]


def _cluster_records(cluster_id: str, member_records: list[dict[str, str]], rules: dict[str, Any]) -> list[dict[str, str]]:
    records = []
    for category in rules["categories"]:
        matching = [record for record in member_records if record["Category"] == category]
        highest_level = max(LEVEL_ORDER[record["Level"]] for record in matching)
        highest = [record for record in matching if LEVEL_ORDER[record["Level"]] == highest_level]
        level = highest[0]["Level"]
        reasons = sorted({record["Reason"] for record in highest})
        actions = sorted({record["Required_Action"] for record in highest if record["Required_Action"]})
        evidence_fields = "; ".join(sorted({record["Evidence_Fields"] for record in highest}))
        reason = f"{len(highest)} of {len(matching)} member(s) have {level} {category.lower()} uncertainty. " + " | ".join(reasons)
        required_action = " | ".join(actions)
        records.append(_record("Cluster", cluster_id, category, level, reason, required_action, evidence_fields, rules["rule_version"]))
    return records


def add_uncertainty_register(df: pd.DataFrame, identity_col: str = "Compound Name", rules_path: str | Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Append compact member summaries and return normalized compound/cluster uncertainty rows."""
    rules = load_uncertainty_rules(rules_path)
    result = df.copy()
    result["_stage3_row_position"] = range(1, len(result) + 1)
    all_records: list[dict[str, str]] = []
    summaries: dict[Any, str] = {}
    groups = result.groupby("Cluster ID", dropna=False, sort=False) if "Cluster ID" in result.columns else [("Unassigned", result)]
    for cluster_id, group in groups:
        cluster_name = _text(cluster_id) or "Unassigned"
        cluster_member_records: list[dict[str, str]] = []
        for index, row in group.iterrows():
            identifier = _text(row.get("InChIKey")) or _text(row.get("CAS")) or _text(row.get(identity_col)) or _text(row.get("Standardized SMILES")) or "Unidentified compound"
            source_position = int(row["_stage3_row_position"])
            row_id = f"{cluster_name} | source row {source_position} | {identifier}"
            member = _member_records(row, row_id, rules)
            cluster_member_records.extend(member)
            all_records.extend(member)
            grouped = {level: [record["Category"] for record in member if record["Level"] == level] for level in LEVEL_ORDER}
            summaries[source_position] = " | ".join(f"{level}: {', '.join(grouped[level])}" for level in ("High", "Medium") if grouped[level]) or "Low uncertainty across Stage 3 categories"
        all_records.extend(_cluster_records(cluster_name, cluster_member_records, rules))
    register = pd.DataFrame(all_records, columns=["Scope", "Scope_ID", "Category", "Level", "Reason", "Required_Action", "Evidence_Fields", "Rule_Version"])
    result["Uncertainty_Summary"] = result["_stage3_row_position"].map(summaries)
    result["Uncertainty_Rule_Version"] = rules["rule_version"]
    result = result.drop(columns=["_stage3_row_position"])
    result.attrs["uncertainty_register"] = register
    return result, register
