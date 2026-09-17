from __future__ import annotations

from typing import Iterable

import pandas as pd
from rdkit import Chem


WORKBOOK_CONTRACT_SECTIONS: list[dict[str, object]] = [
    {
        "section": "Input capture and identity resolution",
        "status": "computed",
        "notes": "Source row, identity, and lookup results that seed the workbook.",
        "columns": [
            "Compound Name",
            "CASRN",
            "Concentration",
            "Source",
            "Cleaned Compound Name",
            "Conc_Float",
            "Resolved Name",
            "Identity Status",
        ],
    },
    {
        "section": "Structure ingestion and normalization",
        "status": "computed",
        "notes": "Structure fields used for all downstream evidence layers.",
        "columns": [
            "Original SMILES",
            "Standardized SMILES",
            "InChIKey",
            "InChI",
            "Classification Input Source",
            "Classification Structure SMILES",
            "Classification Standardization Version",
            "Classification Stereochemistry Status",
            "Structural_Evidence_Status",
            "Structural_Evidence_Version",
            "3D_Evidence_Status",
        ],
    },
    {
        "section": "Primary classification decision",
        "status": "computed",
        "notes": "Core class assignment and audit trail.",
        "columns": [
            "Chemical Class",
            "Corrected Chemical Class",
            "PRISM Structural Class",
            "Final Classification Record",
            "Classification Status",
            "Classification Rule ID",
            "Classification Rule Version",
            "Classification Scope",
            "Classification Review Recommended",
            "Classification Review Recommendation",
            "Classification Suggested Action",
            "Manual Review Flag",
            "Manual Review Reason",
            "Detected Feature Profile",
            "Confidence",
            "Tanimoto Analogs Found",
        ],
    },
    {
        "section": "Functional group and topology interpretation",
        "status": "computed",
        "notes": "Chemistry hierarchy and topology explanation for toxicology review.",
        "columns": [
            "Primary Functional Group",
            "Secondary Functional Groups",
            "Secondary Functional Groups Basis",
            "Taxonomy Path",
            "Taxonomy Path Steps",
            "Matched Categories",
            "Direct Parent",
            "Taxonomy Dictionary Version",
            "Taxonomy Hierarchy",
            "Topology Profile",
            "Reference Taxonomy Path",
        ],
    },
    {
        "section": "Physicochemical descriptor generation",
        "status": "computed",
        "notes": "Calculated descriptor space used for similarity, clustering, and domain logic.",
        "columns": [
            "MW",
            "LogP",
            "TPSA",
            "Min_EState",
            "Max_EState",
            "3D_Asphericity",
            "3D_PMI1",
            "3D_PMI2",
            "3D_PMI3",
            "3D_RadiusOfGyration",
            "HBD",
            "HBA",
            "Rotatable_Bonds",
            "Formal_Charge",
            "Ring_Count",
            "Aromatic_Ring_Count",
            "Fraction_CSP3",
            "Heavy_Atom_Count",
            "Molecular_Refractivity",
            "Ionisation_Indicator",
        ],
    },
    {
        "section": "Toxicophore and structural alert review",
        "status": "computed",
        "notes": "Hazard-relevant substructures and screening context.",
        "columns": [
            "Toxicophore_Profile",
            "All_Structural_Alerts",
            "Alerts",
            "Scaffold",
        ],
    },
    {
        "section": "Cluster and domain preview",
        "status": "preview-only",
        "notes": "Batch clustering uses deduplicated standardized structures, compatibility-gated Morgan/ECFP4 subclusters, and an SME-reviewable Cluster Summary sheet.",
        "columns": [
            "Cluster ID",
            "Compatibility Group",
            "Cluster Size",
            "Cluster Unique Structure Count",
            "Cluster_Status",
            "Domain",
            "Decision",
            "Compatibility Gate",
            "Compatibility Gate Reason",
            "Duplicate Structure Count",
            "Cluster Property Compatibility",
            "Cluster Property Compatibility Reason",
            "Identity Structure Consistency",
            "Identity Structure Consistency Reason",
            "Domain_Status",
            "Domain_Reasons",
            "Domain_Rule_Version",
            "Nearest_Neighbour_Tanimoto",
            "Cluster_Min_Pairwise_Tanimoto",
            "Cluster_Median_Pairwise_Tanimoto",
            "Cluster_Scaffold_Coverage",
            "Cluster_Representative_Scaffold",
            "Property_Outlier_Flags",
            "Ionisation_Consistency",
            "Toxicophore_Consistency",
            "Uncertainty_Summary",
            "Uncertainty_Rule_Version",
            "Cluster_Reasoning_Status",
            "Cluster_Nearest_Members",
            "Cluster_Nearest_Tanimoto",
            "Cluster_Consensus_Method",
            "Cluster_Membership_Parameters",
            "Cluster_Membership_Rationale",
            "SME_Review_Boundary",
            "Cluster_Reasoning_Rule_Version",
        ],
    },
    {
        "section": "Endpoint evidence and external taxonomy",
        "status": "computed",
        "notes": "Supporting evidence layers mirrored from exported workbook sources.",
        "columns": [
            "EPA_CTX_Status",
            "EPA_CTX_DTXSID",
            "EPA_CTX_Source_Version",
            "ChemOnt_Kingdom",
            "ChemOnt_Superclass",
            "ChemOnt_Class",
            "ChemOnt_Subclass",
            "ChemOnt_Direct_Parent",
            "ChemOnt_Molecular_Framework",
            "ChemOnt_Substituents",
            "ChemOnt_Classification_Version",
            "ChemOnt_Source",
            "ChemOnt_Retrieval_Status",
        ],
    },
    {
        "section": "Read-across and final export support",
        "status": "computed",
        "notes": "Workbook-only support columns for downstream review.",
        "columns": [
            "ChEMBL_Max_Phase",
            "ChEMBL_Targets",
            "XAI_Status",
            "AI_Cluster_Reasoning",
            "AI_Evidence_Fields",
            "AI_Input_Hash",
            "AI_Prompt_Version",
            "AI_Response_Timestamp",
        ],
    },
]

# Retired workbook fields are tracked for contract review, but excluded from the
# final exported workbook so the schema does not drift back to the discarded flow.
RETIRED_WORKBOOK_COLUMNS = [
    "Cramer Class",
    "TTC Limit",
    # External services intentionally removed from the current local workflow.
    "EPA_CTX_Status",
    "EPA_CTX_DTXSID",
    "EPA_CTX_Source_Version",
    "ChemOnt_Kingdom",
    "ChemOnt_Superclass",
    "ChemOnt_Class",
    "ChemOnt_Subclass",
    "ChemOnt_Direct_Parent",
    "ChemOnt_Molecular_Framework",
    "ChemOnt_Substituents",
    "ChemOnt_Classification_Version",
    "ChemOnt_Source",
    "ChemOnt_Retrieval_Status",
]


# Reader-facing definitions are emitted at the top of Workbook Contract.  They
# describe the lifecycle of a field, not its toxicological acceptability.
WORKBOOK_STATUS_DEFINITIONS = (
    (
        "Computed",
        "Generated during the current PRISM run from the available input, local rules, and deterministic calculations. "
        "Computed values remain subject to the field-specific uncertainty and SME-review boundaries; they do not establish endpoint equivalence or safety.",
    ),
    (
        "Preview only",
        "A review-support value shown for transparency before an SME decision. It may guide inspection, but it is not an approved conclusion or a release decision.",
    ),
    (
        "Retired",
        "A legacy or intentionally disabled field. It is documented for traceability but excluded from the reader-facing workbook and is not calculated for the current run.",
    ),
    (
        "Cluster identifier convention",
        "A Cluster ID (for example, Cluster_001) identifies membership in one deterministic structural grouping. "
        "Its number is assigned only for stable display and sorting; it is not a similarity value, score, rank, priority, or toxicological measurement. "
        "Where a legacy cluster letter is present, the letter only distinguishes a partition within a broader family and has no measurement meaning.",
    ),
)


def _unique_columns(columns: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for column in columns:
        if column in seen:
            continue
        seen.add(column)
        ordered.append(column)
    return ordered


WORKBOOK_COLUMN_ORDER = _unique_columns(
    column
    for section in WORKBOOK_CONTRACT_SECTIONS
    for column in section["columns"]
    if column not in RETIRED_WORKBOOK_COLUMNS
)


def apply_workbook_contract(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder the workbook export, enforce review-only statuses, and drop retired fields."""
    if df is None or df.empty:
        return df.copy() if isinstance(df, pd.DataFrame) else pd.DataFrame()

    contract_df = df.copy()
    if "Cluster Family" in contract_df.columns:
        if "Compatibility Group" not in contract_df.columns:
            contract_df = contract_df.rename(columns={"Cluster Family": "Compatibility Group"})
        else:
            contract_df = contract_df.drop(columns=["Cluster Family"])
    for retired in RETIRED_WORKBOOK_COLUMNS:
        if retired in contract_df.columns:
            contract_df = contract_df.drop(columns=[retired])

    for column in WORKBOOK_COLUMN_ORDER:
        if column not in contract_df.columns:
            contract_df[column] = ""

    def _valid_structure(value: object) -> bool:
        if value is None or pd.isna(value):
            return False
        text = str(value).strip()
        if not text:
            return False
        try:
            return Chem.MolFromSmiles(text) is not None
        except Exception:
            return False

    # Enforce the structural eligibility invariant at the export boundary as
    # well as in the clustering engine. This prevents stale/legacy rows from
    # being exported as clustered or passed when their structure is missing.
    if "Standardized SMILES" in contract_df.columns:
        # Legacy/merged workbook payloads can contain duplicate header names.
        # Pandas returns a DataFrame (not a Series) for such a selection, and
        # boolean expressions on it raise "truth value is ambiguous".  Use the
        # first occurrence for contract decisions while retaining the existing
        # output column ordering.
        def _contract_series(column):
            selected = contract_df[column]
            return selected.iloc[:, 0] if isinstance(selected, pd.DataFrame) else selected

        def _ensure_object(column):
            selected = contract_df[column]
            if isinstance(selected, pd.DataFrame):
                # Convert each duplicate occurrence independently; assigning
                # ``contract_df[column]`` would itself be a DataFrame and can
                # trigger pandas' ambiguous truth-value path.
                for idx in [i for i, name in enumerate(contract_df.columns) if name == column]:
                    contract_df.iloc[:, idx] = contract_df.iloc[:, idx].astype(object)
            elif not pd.api.types.is_object_dtype(selected):
                contract_df[column] = selected.astype(object)

        valid = _contract_series("Standardized SMILES").map(_valid_structure)
        invalid_indices = contract_df.index[~valid]
        if len(invalid_indices):
            explicit_review = {
                "Cluster ID": "Not assessable",
                "Cluster Size": "Not assessable",
                "Cluster_Status": "Review Required",
                "Domain": "Not assessable",
                "Decision": "REVIEW",
                "Domain_Status": "Not assessable",
                "Domain_Reasons": "Standardized structure missing or invalid",
                "Nearest_Neighbour_Tanimoto": "Not assessable",
                "Cluster_Min_Pairwise_Tanimoto": "Not assessable",
                "Cluster_Median_Pairwise_Tanimoto": "Not assessable",
                "Cluster_Scaffold_Coverage": "Not assessable",
                "Cluster_Representative_Scaffold": "Not assessable",
                "Property_Outlier_Flags": "Not assessable",
                "Ionisation_Consistency": "Not assessable",
                "Toxicophore_Consistency": "Not assessable",
                "Toxicophore_Profile": "Not assessable",
                "All_Structural_Alerts": "Not assessable",
                "Alerts": "Not assessable",
                "Scaffold": "Not assessable",
                "ChEMBL_Max_Phase": "Not assessable",
                "ChEMBL_Targets": "Not assessable",
            }
            for column, value in explicit_review.items():
                if column in contract_df.columns:
                    # Some legacy workbooks infer numeric dtypes (for example
                    # Cluster Size or Tanimoto fields). Explicit review
                    # statuses must still be representable without a pandas
                    # lossy-setitem failure.
                    _ensure_object(column)
                    contract_df.loc[invalid_indices, column] = value

        # Reconcile stale review metadata from earlier identity-resolution
        # passes.  A parseable surrogate is computationally usable, but its
        # identity remains review-worthy; only rows with no parseable structure
        # may use the no-structure reason.
        if "Manual Review Reason" in contract_df.columns:
            reason_series = _contract_series("Manual Review Reason").astype(str).str.strip()
            identity_series = _contract_series("Identity Status").astype(str).str.strip() if "Identity Status" in contract_df.columns else pd.Series("", index=contract_df.index)
            surrogate = identity_series.str.casefold().str.startswith("surrogate")
            old_no_structure = reason_series.str.casefold().eq("no valid standardized or original smiles could be parsed.")
            if not pd.api.types.is_object_dtype(_contract_series("Manual Review Reason")):
                _ensure_object("Manual Review Reason")
            if "Manual Review Flag" in contract_df.columns:
                _ensure_object("Manual Review Flag")
            if len(invalid_indices):
                contract_df.loc[invalid_indices, "Manual Review Reason"] = "No valid standardized or original SMILES could be parsed."
                if "Manual Review Flag" in contract_df.columns:
                    contract_df.loc[invalid_indices, "Manual Review Flag"] = 1
            valid_surrogate = valid & surrogate
            contract_df.loc[valid_surrogate, "Manual Review Reason"] = "Surrogate structure used; identity uncertainty requires toxicologist review."
            if "Manual Review Flag" in contract_df.columns:
                contract_df.loc[valid_surrogate, "Manual Review Flag"] = 1
            # Clear only the stale no-structure reason on otherwise valid,
            # non-surrogate rows; preserve legitimate classification reasons.
            stale_valid = valid & (~surrogate) & old_no_structure
            contract_df.loc[stale_valid, "Manual Review Reason"] = ""
            if "Manual Review Flag" in contract_df.columns:
                contract_df.loc[stale_valid, "Manual Review Flag"] = 0

        # Make no-finding values explicit for valid structures. Numeric fields
        # are intentionally untouched so zero remains a real numeric zero.
        no_finding_columns = ["Toxicophore_Profile", "All_Structural_Alerts", "Alerts"]
        for column in no_finding_columns:
            if column not in contract_df.columns:
                continue
            selected = _contract_series(column)
            missing = selected.isna() | selected.astype(str).str.strip().isin({"", "nan", "None"})
            _ensure_object(column)
            contract_df.loc[valid & missing, column] = "None"
        for column in ["ChEMBL_Max_Phase", "ChEMBL_Targets"]:
            if column not in contract_df.columns:
                continue
            selected = _contract_series(column)
            missing = selected.isna() | selected.astype(str).str.strip().isin({"", "nan"})
            _ensure_object(column)
            contract_df.loc[valid & missing, column] = "None"

    ordered_columns = WORKBOOK_COLUMN_ORDER + [
        column for column in contract_df.columns
        if column not in WORKBOOK_COLUMN_ORDER and column not in RETIRED_WORKBOOK_COLUMNS
    ]
    return contract_df[ordered_columns]


def build_workbook_contract_frame() -> pd.DataFrame:
    rows: list[dict[str, object]] = [
        {
            "Workbook Order": "",
            "Section": "Status definitions",
            "Column Headers - Detailed Analysis sheet": label,
            "Status": "Definition",
            "Notes": definition,
        }
        for label, definition in WORKBOOK_STATUS_DEFINITIONS
    ]
    order = 1
    for section in WORKBOOK_CONTRACT_SECTIONS:
        for column in section["columns"]:
            if column in RETIRED_WORKBOOK_COLUMNS:
                continue
            rows.append(
                {
                    "Workbook Order": order,
                    "Section": section["section"],
                    "Column Headers - Detailed Analysis sheet": column,
                    "Status": section["status"],
                    "Notes": section["notes"],
                }
            )
            order += 1
    for retired in RETIRED_WORKBOOK_COLUMNS:
        rows.append(
            {
                "Workbook Order": "",
                "Section": "Retired workbook fields",
                "Column Headers - Detailed Analysis sheet": retired,
                "Status": "retired",
                "Notes": "Excluded from the final workbook contract.",
            }
        )
    return pd.DataFrame(rows, columns=[
        "Workbook Order", "Section", "Column Headers - Detailed Analysis sheet", "Status", "Notes",
    ])
