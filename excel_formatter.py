import pandas as pd
import xlsxwriter
from io import BytesIO
from datetime import datetime

from rdkit import Chem

from chemont_export import chemont_display_frame
from cluster_summary import build_cluster_summary
from workbook_contract import build_workbook_contract_frame, RETIRED_WORKBOOK_COLUMNS


# ChemOnt and EPA CTX are intentionally out of the current local workflow.
# Keep their backend data available for future reactivation, but do not emit
# empty/placeholder tabs in the user-facing workbook.
EXPORT_EXTERNAL_EVIDENCE = False


# Reader-facing export view.  Keep this explicit and ordered: the Summary
# worksheet is intentionally a compact projection of Detailed Analysis, not a
# second calculation path or an independently maintained dataset.
SUMMARY_COLUMNS = [
    "Compound Name",
    "CASRN",
    "Source",
    "Conc_Float",
    "Original SMILES",
    "Standardized SMILES",
    "InChIKey",
    "InChI",
    "Detected Feature Profile",
    "Primary Functional Group",
    "Secondary Functional Groups",
    "Taxonomy Path",
    "Taxonomy Hierarchy",
    "MW",
    "LogP",
    "TPSA",
    "Ionisation_Indicator",
    "Toxicophore_Profile",
    "All_Structural_Alerts",
    "Alerts",
    "Cluster ID",
    "Cluster Size",
    "Uncertainty_Summary",
    "Cluster_Class_Reasons",
]


# Detailed Analysis is an auditable reader-facing export, not a debug dump of
# every transient pipeline field.  Keep this allow-list explicit so disabled
# integrations and future implementation fields cannot create blank columns in
# the workbook.  Summary remains its compact, stable projection.
DETAILED_ANALYSIS_COLUMNS = [
    "Compound Name", "CASRN", "Concentration", "Source", "Cleaned Compound Name", "Conc_Float", "Resolved Name", "Identity Status",
    "Original SMILES", "Standardized SMILES", "InChIKey", "InChI",
    "Classification Input Source", "Classification Structure SMILES", "Classification Standardization Version", "Classification Stereochemistry Status",
    "Structural_Evidence_Status", "Structural_Evidence_Version",
    "Chemical Class", "Corrected Chemical Class", "PRISM Structural Class", "Final Classification Record", "Classification Status",
    "Classification Rule ID", "Classification Rule Version", "Classification Scope", "Classification Review Recommended",
    "Classification Review Recommendation", "Classification Suggested Action", "Manual Review Flag", "Manual Review Reason",
    "Detected Feature Profile", "Confidence", "Tanimoto Analogs Found", "Primary Functional Group", "Secondary Functional Groups",
    "Secondary Functional Groups Basis", "Taxonomy Path", "Taxonomy Path Steps", "Matched Categories", "Direct Parent",
    "Taxonomy Dictionary Version", "Taxonomy Hierarchy", "Topology Profile",
    "MW", "LogP", "TPSA", "Min_EState", "Max_EState", "HBD", "HBA", "Rotatable_Bonds", "Formal_Charge",
    "Ring_Count", "Aromatic_Ring_Count", "Fraction_CSP3", "Heavy_Atom_Count", "Molecular_Refractivity", "Ionisation_Indicator",
    "Toxicophore_Profile", "All_Structural_Alerts", "Alerts", "Scaffold",
    "Cluster ID", "Compatibility Group", "Cluster Size", "Cluster Unique Structure Count", "Cluster_Status", "Domain", "Decision",
    "Compatibility Gate", "Compatibility Gate Reason", "Cluster Property Compatibility", "Cluster Property Compatibility Reason",
    "Duplicate Structure Count", "Identity Structure Consistency", "Identity Structure Consistency Reason",
    "Domain_Status", "Domain_Reasons", "Domain_Rule_Version", "Nearest_Neighbour_Tanimoto", "Cluster_Min_Pairwise_Tanimoto",
    "Cluster_Median_Pairwise_Tanimoto", "Cluster_Scaffold_Coverage", "Cluster_Representative_Scaffold", "Property_Outlier_Flags",
    "Ionisation_Consistency", "Toxicophore_Consistency", "Uncertainty_Summary", "Uncertainty_Rule_Version", "Cluster_Reasoning_Status",
    "Cluster_Nearest_Members", "Cluster_Nearest_Tanimoto", "Cluster_Consensus_Method", "Cluster_Membership_Parameters",
    "Cluster_Membership_Rationale", "SME_Review_Boundary", "Cluster_Reasoning_Rule_Version",
    "XAI_Status", "AI_Cluster_Reasoning", "AI_Evidence_Fields", "AI_Input_Hash", "AI_Prompt_Version", "AI_Response_Timestamp",
    "AI_Validation_Status", "AI_Response_Hash",
    # Keep only the compact ontology crosswalk needed for SME traceability.
    # The removed hierarchy/source fields are derivable duplicates or static
    # reference metadata and remain available in the underlying pipeline.
    "Classification Rationale", "Ontology Mapping Status", "Ontology Mapping Version", "Ontology Mapping Evidence",
    "Local Taxonomy Match Status", "Local Taxonomy Record ID", "Reference Taxonomy Path", "Local Taxonomy Reference",
]


def _cluster_class_reason(row: pd.Series) -> str:
    """Create the requested compact Summary explanation from deterministic evidence."""
    raw_existing = row.get("Cluster_Class_Reasons", "")
    existing = "" if pd.isna(raw_existing) else str(raw_existing or "").strip()
    if existing and existing.casefold() not in {"nan", "none", "not available"}:
        return existing
    cluster_id = str(row.get("Cluster ID", "") or "").strip()
    if not cluster_id or cluster_id.casefold() == "not assessable":
        return "Not assessable: no valid standardized structure is available for structural clustering."
    gate = str(row.get("Compatibility Gate", "") or "").strip()
    compatibility = str(row.get("Cluster Property Compatibility Reason", "") or "").strip()
    if gate and compatibility:
        return f"{gate}. {compatibility}"
    return gate or "Cluster assignment requires SME review."


def _first_matching_series(frame: pd.DataFrame, column: str) -> pd.Series | None:
    if column not in frame.columns:
        return None
    subset = frame.loc[:, frame.columns == column]
    if isinstance(subset, pd.Series):
        return subset
    if subset.shape[1] == 1:
        return subset.iloc[:, 0]
    return subset.apply(
        lambda row: next(
            (value for value in row.tolist() if not pd.isna(value) and str(value).strip()),
            next((value for value in row.tolist() if not pd.isna(value)), "")
        ),
        axis=1,
    )


def _safe_column_width(frame: pd.DataFrame, column: str, cap: int = 40, floor: int = 12) -> int:
    series = _first_matching_series(frame, column)
    if series is None or series.empty:
        return max(floor, min(cap, len(str(column)) + 2))
    max_len = series.map(lambda value: len(str(value))).max()
    return min(cap, max(floor, max(len(str(column)), int(max_len)) + 2))


def _first_column_index(frame: pd.DataFrame, column: str) -> int | None:
    matches = [idx for idx, name in enumerate(frame.columns) if name == column]
    return matches[0] if matches else None


def _flatten_metadata(value, prefix=""):
    rows = []
    if isinstance(value, dict):
        for key, nested_value in value.items():
            next_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(_flatten_metadata(nested_value, next_prefix))
    elif isinstance(value, list):
        if not value:
            rows.append({"Path": prefix, "Value": ""})
        else:
            for index, nested_value in enumerate(value):
                next_prefix = f"{prefix}[{index}]"
                rows.extend(_flatten_metadata(nested_value, next_prefix))
    else:
        rows.append({"Path": prefix, "Value": "" if value is None else value})
    return rows


def _usable_structure_mask(frame: pd.DataFrame) -> pd.Series:
    """Return a conservative structure-eligibility mask for dashboard counts."""
    smiles = _first_matching_series(frame, "Standardized SMILES")
    if smiles is None:
        return pd.Series(False, index=frame.index, dtype=bool)

    def _is_usable(value: object) -> bool:
        if value is None or pd.isna(value):
            return False
        text = str(value).strip()
        if not text or text.casefold() in {"nan", "none", "not available", "not assessable"}:
            return False
        try:
            return Chem.MolFromSmiles(text) is not None
        except Exception:
            return False

    return smiles.map(_is_usable)


def _nonempty_unique_count(series: pd.Series | None) -> int:
    if series is None:
        return 0
    values = series.astype(str).str.strip()
    return int(values.loc[~values.str.casefold().isin({"", "nan", "none", "not available", "not assessable"})].nunique())


def _build_dashboard_data(df: pd.DataFrame, raw_df: pd.DataFrame | None, cluster_summary_df: pd.DataFrame) -> tuple[list[tuple[str, int]], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build display-only dashboard counts from the same frames exported below."""
    usable_structure = _usable_structure_mask(df)
    casrn = _first_matching_series(df, "CASRN")
    cluster_id = _first_matching_series(df, "Cluster ID")
    compatibility = _first_matching_series(df, "Compatibility Group")
    if compatibility is None:
        compatibility = _first_matching_series(df, "Cluster Family")
    manual_flag = _first_matching_series(df, "Manual Review Flag")
    manual_reason = _first_matching_series(df, "Manual Review Reason")
    alerts = _first_matching_series(df, "Alerts")

    if cluster_id is None:
        clustered = pd.Series(False, index=df.index, dtype=bool)
        cluster_values = pd.Series("", index=df.index, dtype=object)
    else:
        cluster_values = cluster_id.astype(str).str.strip()
        clustered = usable_structure & ~cluster_values.str.casefold().isin({"", "nan", "none", "not assessable", "not available"})

    if manual_flag is None:
        manual = pd.Series(False, index=df.index, dtype=bool)
    else:
        manual = manual_flag.astype(str).str.strip().str.casefold().isin({"1", "true", "yes", "y", "review", "required"})
    if manual_reason is not None:
        reason_text = manual_reason.astype(str).str.strip().str.casefold()
        manual = manual | ~reason_text.isin({"", "nan", "none", "not available"})

    if alerts is None:
        alert_mask = pd.Series(False, index=df.index, dtype=bool)
    else:
        alert_mask = ~alerts.astype(str).str.strip().str.casefold().isin({"", "nan", "none", "not available", "not assessed", "not assessable"})

    valid_cluster_ids = cluster_values.loc[clustered]
    cluster_sizes = valid_cluster_ids.value_counts()
    if "Cluster Size" in df.columns:
        exported_sizes = pd.to_numeric(_first_matching_series(df, "Cluster Size"), errors="coerce")
        singleton_clusters = int(valid_cluster_ids.loc[exported_sizes.loc[clustered].eq(1)].nunique())
    else:
        singleton_clusters = int((cluster_sizes == 1).sum())

    metrics = [
        ("Raw extracted source rows", int(len(raw_df)) if raw_df is not None else 0),
        ("Final deduplicated chemical records", int(len(df))),
        ("Unique reported CASRNs", _nonempty_unique_count(casrn)),
        ("Unique resolved standardised structures", _nonempty_unique_count(_first_matching_series(df.loc[usable_structure], "Standardized SMILES"))),
        ("Structure-resolved records", int(usable_structure.sum())),
        ("Unresolved / non-structural records", int((~usable_structure).sum())),
        ("Grouping-eligible records", int(usable_structure.sum())),
        ("Compatibility groups", _nonempty_unique_count(compatibility.loc[usable_structure] if compatibility is not None else None)),
        ("Structural clusters", int(valid_cluster_ids.nunique())),
        ("Clustered records", int(clustered.sum())),
        ("Singleton clusters", singleton_clusters),
        ("SME review required", int(manual.sum())),
        ("Records with structural alerts", int(alert_mask.sum())),
    ]

    def _distribution(series: pd.Series | None, label: str) -> pd.DataFrame:
        if series is None:
            return pd.DataFrame(columns=[label, "Record Count"])
        values = series.astype(str).str.strip().replace({"": "Not available", "nan": "Not available", "None": "Not available"})
        return values.value_counts().rename_axis(label).reset_index(name="Record Count").head(12)

    compatibility_distribution = _distribution(compatibility.loc[usable_structure] if compatibility is not None else None, "Compatibility Group")
    class_distribution = _distribution(_first_matching_series(df, "Corrected Chemical Class"), "Corrected Chemical Class")
    resolution_distribution = pd.DataFrame([
        {"Identity / Structure Status": "Structure resolved", "Record Count": int(usable_structure.sum())},
        {"Identity / Structure Status": "Unresolved or non-structural", "Record Count": int((~usable_structure).sum())},
        {"Identity / Structure Status": "SME review required", "Record Count": int(manual.sum())},
    ])
    return metrics, compatibility_distribution, class_distribution, resolution_distribution


def _write_sme_review_guide(workbook, header_format, cell_format, alt_row_format):
    """Create the concise, workbook-native SME instruction sheet."""
    worksheet = workbook.add_worksheet("SME Review Guide")
    worksheet.hide_gridlines(2)
    title_format = workbook.add_format({'bold': True, 'font_size': 16, 'font_color': '#203764'})
    subtitle_format = workbook.add_format({'italic': True, 'font_color': '#595959', 'text_wrap': True, 'valign': 'top'})
    warning_format = workbook.add_format({'bold': True, 'font_color': '#9C0006', 'bg_color': '#FCE4D6', 'text_wrap': True, 'valign': 'top', 'border': 1, 'border_color': '#F4B183'})
    worksheet.merge_range('A1:D1', 'PRISM SME Workbook Review Guide', title_format)
    worksheet.merge_range('A2:D3', 'Use this guide to review the exported chemical records, deterministic classifications, and structural grouping. Detailed Analysis and Cluster Summary remain the authoritative evidence sheets.', subtitle_format)
    worksheet.merge_range('A5:D6', 'Review boundary: PRISM supports structured chemical grouping. It does not establish endpoint equivalence, read-across acceptability, exposure equivalence, safety, or a final toxicological conclusion. SME review remains required.', warning_format)

    rows = [
        ('1', 'Read the Workbook Contract', 'Confirm field definitions, status values, and calculation boundaries before interpreting results.'),
        ('2', 'Verify Raw Extraction (Pre-Dedup)', 'Check priority compound names, CASRNs, sources, pages/tables, concentrations, and duplicate/alias handling against the analytical report.'),
        ('3', 'Triage Summary', 'Filter Manual Review Flag first, then Uncertainty Summary, Cluster ID, Cluster Size, and Compound Name.'),
        ('4', 'Confirm identity and structure', 'Review Compound Name, CASRN, Standardized SMILES, InChIKey, identity status, and source evidence. Do not use missing/invalid structures for analogue evidence.'),
        ('5', 'Review classification', 'Assess Corrected Chemical Class, functional groups, taxonomy, toxicophore profile, and structural alerts for chemical plausibility.'),
        ('6', 'Review Cluster Summary', 'Check Compatibility Group, cluster cohesion, scaffold/ionisation/alert compatibility, property compatibility, membership rationale, and stability.'),
        ('7', 'Record the SME outcome', 'Use: Accept for analogue review; Accept with limitations; Do not use for analogue review; or SME review/follow-up required.'),
    ]
    worksheet.write_row(8, 0, ['Step', 'Review activity', 'What to do'], header_format)
    for offset, row in enumerate(rows, start=9):
        worksheet.write(offset, 0, row[0], alt_row_format if offset % 2 else cell_format)
        worksheet.write(offset, 1, row[1], alt_row_format if offset % 2 else cell_format)
        worksheet.write(offset, 2, row[2], alt_row_format if offset % 2 else cell_format)

    worksheet.write_row(18, 0, ['Status', 'Meaning for review'], header_format)
    definitions = [
        ('Computed', 'Calculated in the current PRISM run; still subject to field-specific uncertainty and SME review.'),
        ('Preview only', 'Review-support information, not an approved conclusion or release decision.'),
        ('Retired', 'Historical/disabled field retained only for traceability; do not use in a new decision.'),
        ('Not assessable', 'A defensible calculation could not be made from available information; it is not a negative result.'),
        ('Unclassified', 'No sufficiently supported structural class is available; review identity/source evidence.'),
        ('SME Review Required', 'A limitation or boundary requires documented SME judgement before use.'),
    ]
    for offset, row in enumerate(definitions, start=19):
        worksheet.write(offset, 0, row[0], alt_row_format if offset % 2 else cell_format)
        worksheet.merge_range(offset, 1, offset, 3, row[1], alt_row_format if offset % 2 else cell_format)

    worksheet.write_row(27, 0, ['SME decision record', 'Complete for priority compound or cluster review'], header_format)
    worksheet.write_row(28, 0, ['Reviewer', 'Review date', 'Endpoint / intended use', 'Outcome', 'Rationale and follow-up'], header_format)
    for row in range(29, 34):
        worksheet.write_row(row, 0, ['', '', '', '', ''], alt_row_format if row % 2 else cell_format)
    worksheet.set_column('A:A', 14)
    worksheet.set_column('B:B', 30)
    worksheet.set_column('C:C', 62)
    worksheet.set_column('D:D', 24)
    worksheet.freeze_panes(8, 0)
    return worksheet


def _write_dashboard(workbook, metrics, compatibility_distribution, class_distribution, resolution_distribution):
    """Create a display-only, auditable overview driven by export data."""
    worksheet = workbook.add_worksheet("Dashboard")
    worksheet.hide_gridlines(2)
    title_format = workbook.add_format({'bold': True, 'font_size': 16, 'font_color': '#203764'})
    subtitle_format = workbook.add_format({'italic': True, 'font_color': '#595959'})
    metric_label_format = workbook.add_format({'bold': True, 'font_color': '#203764', 'bg_color': '#D9EAF7', 'border': 1, 'border_color': '#8EA9DB', 'text_wrap': True, 'valign': 'vcenter'})
    metric_value_format = workbook.add_format({'bold': True, 'font_size': 14, 'align': 'center', 'bg_color': '#F7FBFF', 'border': 1, 'border_color': '#8EA9DB', 'num_format': '#,##0'})
    note_format = workbook.add_format({'italic': True, 'font_color': '#595959', 'text_wrap': True, 'valign': 'top'})
    header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'align': 'center', 'fg_color': '#203764', 'font_color': 'white', 'border': 1, 'border_color': '#8EA9DB'})
    cell_format = workbook.add_format({'valign': 'top', 'text_wrap': True, 'border': 1, 'border_color': '#A6A6A6'})

    worksheet.merge_range('A1:L1', 'PRISM Workbook Dashboard', title_format)
    worksheet.merge_range('A2:L2', 'Overview generated from the same records exported in this workbook. Detailed Analysis and Cluster Summary remain authoritative.', subtitle_format)
    worksheet.write('J3', 'Generated', note_format)
    worksheet.write('K3', datetime.now().strftime('%Y-%m-%d %H:%M:%S'), note_format)

    for index, (label, value) in enumerate(metrics):
        row = 4 + (index // 4) * 2
        col = (index % 4) * 3
        worksheet.merge_range(row, col, row, col + 1, label, metric_label_format)
        worksheet.merge_range(row + 1, col, row + 1, col + 1, value, metric_value_format)

    worksheet.merge_range('A13:H14', 'Definitions: Raw data = all extracted source rows before duplicate handling. Unique chemical records = final deduplicated output records; unique CASRNs and standardised structures are shown separately. Grouping-eligible = usable standardised structures. Clustered = grouping-eligible records assigned to a structural cluster; singleton clusters are counted separately.', note_format)

    def _write_distribution(start_row, start_col, title, frame):
        worksheet.write(start_row, start_col, title, header_format)
        worksheet.write_row(start_row + 1, start_col, list(frame.columns), header_format)
        for offset, values in enumerate(frame.itertuples(index=False, name=None), start=start_row + 2):
            worksheet.write_row(offset, start_col, list(values), cell_format)
        return max(start_row + 2, start_row + 1 + len(frame))

    compatibility_end = _write_distribution(16, 0, 'Compounds by compatibility group', compatibility_distribution)
    class_end = _write_distribution(16, 4, 'Compounds by corrected chemical class', class_distribution)
    resolution_end = _write_distribution(16, 8, 'Identity and review status', resolution_distribution)

    if not compatibility_distribution.empty:
        chart = workbook.add_chart({'type': 'bar'})
        chart.add_series({
            'name': 'Records',
            'categories': ['Dashboard', 18, 0, compatibility_end, 0],
            'values': ['Dashboard', 18, 1, compatibility_end, 1],
            'fill': {'color': '#4472C4'},
        })
        chart.set_title({'name': 'Compounds by compatibility group'})
        chart.set_legend({'none': True})
        chart.set_x_axis({'name': 'Record count', 'min': 0})
        chart.set_y_axis({'reverse': True})
        worksheet.insert_chart('A32', chart, {'x_scale': 1.15, 'y_scale': 1.15})
    if not resolution_distribution.empty:
        chart = workbook.add_chart({'type': 'doughnut'})
        chart.add_series({
            'name': 'Records',
            'categories': ['Dashboard', 18, 8, resolution_end, 8],
            'values': ['Dashboard', 18, 9, resolution_end, 9],
            'data_labels': {'value': True},
        })
        chart.set_title({'name': 'Identity and review status'})
        chart.set_legend({'position': 'bottom'})
        worksheet.insert_chart('H32', chart, {'x_scale': 1.05, 'y_scale': 1.15})

    worksheet.set_column('A:A', 26)
    worksheet.set_column('B:B', 14)
    worksheet.set_column('C:C', 3)
    worksheet.set_column('D:D', 3)
    worksheet.set_column('E:E', 30)
    worksheet.set_column('F:F', 14)
    worksheet.set_column('G:G', 3)
    worksheet.set_column('H:H', 3)
    worksheet.set_column('I:I', 32)
    worksheet.set_column('J:J', 14)
    worksheet.set_column('K:L', 16)
    return worksheet


def format_excel_output(df: pd.DataFrame, out_buffer: BytesIO, raw_tables: list = None):
    """Format the Phase B workbook with one deterministic compound sheet."""
    source_attrs = dict(getattr(df, 'attrs', {}))
    uncertainty_register = source_attrs.get('uncertainty_register')
    endpoint_evidence_ledger = source_attrs.get('endpoint_evidence_ledger')
    cluster_evidence_matrix = source_attrs.get('cluster_evidence_matrix')
    chemont_classifications = source_attrs.get('chemont_classifications')
    workbook_contract = source_attrs.get('workbook_contract')
    clustering_policy = source_attrs.get('clustering_policy')
    class_layer_policy = source_attrs.get('class_layer_policy')
    descriptor_layer_policy = source_attrs.get('descriptor_layer_policy')
    topology_layer_policy = source_attrs.get('topology_layer_policy')
    ionisation_layer_policy = source_attrs.get('ionisation_layer_policy')
    toxicophore_layer_policy = source_attrs.get('toxicophore_layer_policy')
    chemont_layer_policy = source_attrs.get('chemont_layer_policy')
    distance_layer_metadata = source_attrs.get('distance_layer_metadata')
    df = df.copy()
    # Enforce the public export schema even when callers invoke the formatter
    # directly without first applying the workbook contract.
    drop_columns = [column for column in RETIRED_WORKBOOK_COLUMNS if column in df.columns]
    if drop_columns:
        df = df.drop(columns=drop_columns)
    df["Cluster_Class_Reasons"] = df.apply(_cluster_class_reason, axis=1)
    df.attrs.update(source_attrs)
    
    # Deterministic SME-review ordering for the reader-facing sheets.
    # Compatibility groups and cluster IDs such as ``Family_004`` and
    # ``Cluster_004`` sort numerically; unresolved/non-clustered rows remain
    # at the end.  Retain the legacy column as a fallback for old callers.
    compatibility_column = (
        'Compatibility Group' if 'Compatibility Group' in df.columns
        else 'Cluster Family' if 'Cluster Family' in df.columns
        else None
    )
    if compatibility_column:
        compatibility_text = df[compatibility_column].astype(str)
        df['_sort_compatibility_group'] = pd.to_numeric(
            compatibility_text.str.extract(r'(\d+)', expand=False), errors='coerce'
        ).fillna(999999)
    else:
        df['_sort_compatibility_group'] = 999999
    if 'Cluster ID' in df.columns:
        cluster_text = df['Cluster ID'].astype(str)
        df['_sort_cluster'] = pd.to_numeric(
            cluster_text.str.extract(r'(\d+)', expand=False), errors='coerce'
        ).fillna(999999)
    else:
        df['_sort_cluster'] = 999999
    df['_sort_clustered'] = (df['_sort_cluster'] < 999999).astype(int)
    cluster_sizes = df['Cluster Size'] if 'Cluster Size' in df.columns else pd.Series(0, index=df.index)
    df['_sort_cluster_size'] = pd.to_numeric(cluster_sizes, errors='coerce').fillna(0)
    sort_columns = [
        '_sort_clustered', '_sort_compatibility_group', '_sort_cluster',
        '_sort_cluster_size',
    ]
    sort_ascending = [False, True, True, False]
    if 'Compound Name' in df.columns:
        sort_columns.append('Compound Name')
        sort_ascending.append(True)
    df = df.sort_values(by=sort_columns, ascending=sort_ascending, kind='mergesort')
    df = df.drop(columns=[
        '_sort_clustered', '_sort_compatibility_group', '_sort_cluster',
        '_sort_cluster_size',
    ])

    # The Summary is a direct, ordered projection of the complete detailed
    # sheet.  reindex keeps the public schema stable and presents blank cells
    # when an optional upstream field was not calculated in a particular run.
    summary_df = df.reindex(columns=SUMMARY_COLUMNS).copy()
    # A single auditable detailed sheet keeps high-confidence, analogue, and
    # unresolved rows together so users can filter/sort without losing context.
    # Cluster Summary is a review layer only; it does not independently assign
    # membership and always traces back to Detailed Analysis.
    cluster_summary_df = build_cluster_summary(df)
    detailed_df = df.reindex(columns=DETAILED_ANALYSIS_COLUMNS).copy()
    sheets_to_write = (
        [('Summary', summary_df), ('Cluster Summary', cluster_summary_df), ('Detailed Analysis', detailed_df)]
        if not df.empty else []
    )

    writer = pd.ExcelWriter(out_buffer, engine='xlsxwriter')
    workbook = writer.book
    
    # Define stunning formats
    header_format = workbook.add_format({
        'bold': True,
        'text_wrap': True,
        'valign': 'top',
        'align': 'center',
        'fg_color': '#203764',   # Dark blue
        'font_color': 'white',
        'border': 1,
        'border_color': '#8EA9DB'
    })
    
    cell_format = workbook.add_format({
        'valign': 'top',
        'text_wrap': True,
        'border': 1,
        'border_color': '#A6A6A6'
    })
    
    alt_row_format = workbook.add_format({
        'valign': 'top',
        'text_wrap': True,
        'border': 1,
        'border_color': '#A6A6A6',
        'bg_color': '#F2F2F2'
    })
    
    green_format = workbook.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})
    red_format = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
    yellow_format = workbook.add_format({'bg_color': '#FFEB9C', 'font_color': '#9C6500'})

    # Reader-facing guide and overview are placed first.  They are generated
    # solely from the same source, detailed, and cluster frames exported below;
    # neither worksheet feeds any calculation or modifies a scientific result.
    combined_raw = pd.concat(raw_tables, ignore_index=True) if raw_tables else None
    metrics, compatibility_distribution, class_distribution, resolution_distribution = _build_dashboard_data(
        df, combined_raw, cluster_summary_df
    )
    _write_sme_review_guide(workbook, header_format, cell_format, alt_row_format)
    _write_dashboard(workbook, metrics, compatibility_distribution, class_distribution, resolution_distribution)

    # The contract precedes the source and analysis sheets so field definitions
    # are immediately available to reviewers.
    if not isinstance(workbook_contract, pd.DataFrame) or workbook_contract.empty:
        workbook_contract = build_workbook_contract_frame()
    elif (
        "Column" in workbook_contract.columns
        and "Column Headers - Detailed Analysis sheet" not in workbook_contract.columns
    ):
        # Older saved run-state payloads may carry the former generic header.
        # Normalize it at export so every final workbook uses the same contract.
        workbook_contract = workbook_contract.rename(columns={
            "Column": "Column Headers - Detailed Analysis sheet"
        })
    if isinstance(workbook_contract, pd.DataFrame) and not workbook_contract.empty:
        workbook_contract.to_excel(writer, index=False, sheet_name='Workbook Contract')
        contract_worksheet = writer.sheets['Workbook Contract']
        for col_num, value in enumerate(workbook_contract.columns.values):
            contract_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(workbook_contract)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(workbook_contract.iloc[row_num]):
                contract_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        for col_num, column in enumerate(workbook_contract.columns):
            contract_worksheet.set_column(col_num, col_num, min(80, _safe_column_width(workbook_contract, column, cap=80, floor=14)))
        contract_worksheet.freeze_panes(1, 0)
        contract_worksheet.autofilter(0, 0, len(workbook_contract), len(workbook_contract.columns) - 1)

    # The raw PDF extraction is an immutable source view. It is displayed
    # after the review guide, dashboard, and contract, and never feeds
    # downstream calculations in the exported workbook.
    if combined_raw is not None:
        raw_sheet_name = 'Raw Extraction (Pre-Dedup)'
        combined_raw.to_excel(writer, index=False, sheet_name=raw_sheet_name)
        raw_worksheet = writer.sheets[raw_sheet_name]
        for col_num, value in enumerate(combined_raw.columns.values):
            raw_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(combined_raw)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num in range(len(combined_raw.columns)):
                value = combined_raw.iloc[row_num, col_num]
                try:
                    is_na = pd.isna(value)
                    if isinstance(is_na, (list, tuple)) or type(is_na).__name__ == 'ndarray':
                        is_na = False
                except Exception:
                    is_na = False
                if is_na:
                    value = ""
                elif isinstance(value, (list, dict, tuple)) or type(value).__name__ == 'ndarray':
                    value = str(value)
                raw_worksheet.write(row_num + 1, col_num, value, fmt)
        for col_num, col_name in enumerate(combined_raw.columns):
            raw_worksheet.set_column(col_num, col_num, _safe_column_width(combined_raw, col_name))
        raw_worksheet.freeze_panes(1, 0)
        raw_worksheet.autofilter(0, 0, len(combined_raw), len(combined_raw.columns) - 1)

    for sheet_name, sheet_df in sheets_to_write:
        sheet_df.to_excel(writer, index=False, sheet_name=sheet_name)
        worksheet = writer.sheets[sheet_name]
        
        # Headers
        for col_num, value in enumerate(sheet_df.columns.values):
            worksheet.write(0, col_num, value, header_format)
            
        # Data
        for row_num in range(len(sheet_df)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num in range(len(sheet_df.columns)):
                val = sheet_df.iloc[row_num, col_num]
                
                is_na = False
                try:
                    is_na = pd.isna(val)
                    if isinstance(is_na, (list, tuple)) or type(is_na).__name__ == 'ndarray':
                        is_na = False
                except:
                    pass
                    
                if is_na:
                    val = ""
                elif isinstance(val, (list, dict, tuple)) or type(val).__name__ == 'ndarray':
                    val = str(val)
                    
                worksheet.write(row_num + 1, col_num, val, fmt)
        
        # Column widths
        for col_num, col_name in enumerate(sheet_df.columns):
            worksheet.set_column(col_num, col_num, _safe_column_width(sheet_df, col_name))
            
        worksheet.freeze_panes(1, 0)
        worksheet.autofilter(0, 0, len(sheet_df), len(sheet_df.columns) - 1)
        
        # Conditional formatting
        if 'Confidence' in sheet_df.columns:
            conf_col = _first_column_index(sheet_df, 'Confidence')
            if conf_col is not None:
                worksheet.conditional_format(1, conf_col, len(sheet_df), conf_col,
                                             {'type': 'cell', 'criteria': '==', 'value': '"High"', 'format': green_format})
                worksheet.conditional_format(1, conf_col, len(sheet_df), conf_col,
                                             {'type': 'cell', 'criteria': '==', 'value': '"Medium"', 'format': yellow_format})
                worksheet.conditional_format(1, conf_col, len(sheet_df), conf_col,
                                             {'type': 'cell', 'criteria': '==', 'value': '"Low"', 'format': red_format})
                                         
        if 'Cluster ID' in sheet_df.columns:
            cluster_col = _first_column_index(sheet_df, 'Cluster ID')
            if cluster_col is not None:
                worksheet.conditional_format(1, cluster_col, len(sheet_df), cluster_col,
                                             {'type': '3_color_scale',
                                              'mid_color': '#FFEB84',
                                             'max_color': '#63BE7B'})

        if 'SME Review Status' in sheet_df.columns:
            review_col = _first_column_index(sheet_df, 'SME Review Status')
            if review_col is not None:
                worksheet.conditional_format(1, review_col, len(sheet_df), review_col,
                                             {'type': 'cell', 'criteria': '==', 'value': '"Pending"', 'format': yellow_format})

    # Stage 3 normalized, review-only uncertainty evidence. This is deliberately
    # separate from the processed compound sheets so every action remains auditable.
    if isinstance(uncertainty_register, pd.DataFrame) and not uncertainty_register.empty:
        uncertainty_register.to_excel(writer, index=False, sheet_name='Uncertainty Register')
        uncertainty_worksheet = writer.sheets['Uncertainty Register']
        for col_num, value in enumerate(uncertainty_register.columns.values):
            uncertainty_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(uncertainty_register)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(uncertainty_register.iloc[row_num]):
                uncertainty_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        for col_num, column in enumerate(uncertainty_register.columns):
            uncertainty_worksheet.set_column(col_num, col_num, min(60, _safe_column_width(uncertainty_register, column, cap=60, floor=14)))
        uncertainty_worksheet.freeze_panes(1, 0)
        uncertainty_worksheet.autofilter(0, 0, len(uncertainty_register), len(uncertainty_register.columns) - 1)
        if 'Level' in uncertainty_register.columns:
            level_col = _first_column_index(uncertainty_register, 'Level')
            if level_col is not None:
                uncertainty_worksheet.conditional_format(1, level_col, len(uncertainty_register), level_col, {'type': 'cell', 'criteria': '==', 'value': '"High"', 'format': red_format})
                uncertainty_worksheet.conditional_format(1, level_col, len(uncertainty_register), level_col, {'type': 'cell', 'criteria': '==', 'value': '"Medium"', 'format': yellow_format})
                uncertainty_worksheet.conditional_format(1, level_col, len(uncertainty_register), level_col, {'type': 'cell', 'criteria': '==', 'value': '"Low"', 'format': green_format})

    # Stage 6 deterministic packet used for grounded AI drafting. It is always
    # a separate audit sheet; AI availability must never hide this evidence.
    matrix_columns = ['Cluster_ID', 'Evidence_ID', 'Uncertainty_ID', 'Evidence_Domain', 'Evidence_Type', 'Evidence_Text', 'Source_Fields', 'Matrix_Version', 'Matrix_Status']
    if isinstance(cluster_evidence_matrix, pd.DataFrame) and not cluster_evidence_matrix.empty:
        matrix_for_export = cluster_evidence_matrix.copy()
        matrix_for_export['Matrix_Status'] = 'Deterministic evidence row'
    else:
        matrix_for_export = pd.DataFrame([{'Matrix_Status': 'Stage not run: no cluster evidence matrix available'}], columns=matrix_columns)
    matrix_for_export.to_excel(writer, index=False, sheet_name='Cluster Evidence Matrix')
    matrix_worksheet = writer.sheets['Cluster Evidence Matrix']
    for col_num, value in enumerate(matrix_for_export.columns.values):
        matrix_worksheet.write(0, col_num, value, header_format)
    for row_num in range(len(matrix_for_export)):
        fmt = alt_row_format if row_num % 2 == 1 else cell_format
        for col_num, value in enumerate(matrix_for_export.iloc[row_num]):
            matrix_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
    for col_num, column in enumerate(matrix_for_export.columns):
        matrix_worksheet.set_column(col_num, col_num, min(80, _safe_column_width(matrix_for_export, column, cap=80, floor=14)))
    matrix_worksheet.freeze_panes(1, 0)
    matrix_worksheet.autofilter(0, 0, len(matrix_for_export), len(matrix_for_export.columns) - 1)

    # Stage 7 structural taxonomy is exported separately from the processed
    # sheets, so hierarchy provenance remains reviewable even when live lookup
    # is disabled or unavailable.
    chemont_columns = ['InChIKey', 'Standardized SMILES', 'ChemOnt_Kingdom', 'ChemOnt_Superclass', 'ChemOnt_Class', 'ChemOnt_Subclass', 'ChemOnt_Direct_Parent', 'ChemOnt_Molecular_Framework', 'ChemOnt_Substituents', 'ChemOnt_Classification_Version', 'ChemOnt_Source', 'ChemOnt_Retrieval_Status']
    if EXPORT_EXTERNAL_EVIDENCE and isinstance(chemont_classifications, pd.DataFrame) and not chemont_classifications.empty:
        chemont_for_export = chemont_display_frame(chemont_classifications.copy())
    else:
        chemont_for_export = chemont_display_frame(pd.DataFrame([{'ChemOnt_Retrieval_Status': 'Stage not run: no ChemOnt classification available'}], columns=chemont_columns))
    if EXPORT_EXTERNAL_EVIDENCE:
        chemont_for_export.to_excel(writer, index=False, sheet_name='ChemOnt Classification')
        chemont_worksheet = writer.sheets['ChemOnt Classification']
        for col_num, value in enumerate(chemont_for_export.columns.values):
            chemont_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(chemont_for_export)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(chemont_for_export.iloc[row_num]):
                chemont_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        for col_num, column in enumerate(chemont_for_export.columns):
            chemont_worksheet.set_column(col_num, col_num, min(70, _safe_column_width(chemont_for_export, column, cap=70, floor=14)))
        chemont_worksheet.freeze_panes(1, 0)
        chemont_worksheet.autofilter(0, 0, len(chemont_for_export), len(chemont_for_export.columns) - 1)

    evidence_columns = ['evidence_id', 'identity_key', 'compound_mapping_method', 'endpoint', 'result', 'study_type', 'species', 'route', 'dose', 'dose_unit', 'duration', 'source', 'source_record', 'source_url', 'source_database_version', 'reliability', 'evidence_status', 'raw_record_hash', 'vocabulary_version', 'retrieved_at', 'imported_at', 'imported_by', 'reviewer_status', 'supersedes_evidence_id', 'notes', 'Ledger_Status']
    if isinstance(endpoint_evidence_ledger, pd.DataFrame) and not endpoint_evidence_ledger.empty:
        ledger_for_export = endpoint_evidence_ledger.copy()
        ledger_for_export['Ledger_Status'] = 'Matched local evidence record'
    else:
        ledger_for_export = pd.DataFrame([{'Ledger_Status': 'Stage not run: No matching local evidence records'}], columns=evidence_columns)
    if EXPORT_EXTERNAL_EVIDENCE and not ledger_for_export.empty:
        ledger_for_export.to_excel(writer, index=False, sheet_name='Endpoint Evidence Ledger')
        evidence_worksheet = writer.sheets['Endpoint Evidence Ledger']
        for col_num, value in enumerate(ledger_for_export.columns.values):
            evidence_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(ledger_for_export)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(ledger_for_export.iloc[row_num]):
                evidence_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        for col_num, column in enumerate(ledger_for_export.columns):
            evidence_worksheet.set_column(col_num, col_num, min(60, _safe_column_width(ledger_for_export, column, cap=60, floor=14)))
        evidence_worksheet.freeze_panes(1, 0)
        evidence_worksheet.autofilter(0, 0, len(ledger_for_export), len(ledger_for_export.columns) - 1)

    if isinstance(clustering_policy, dict) and clustering_policy:
        policy_rows = _flatten_metadata(clustering_policy)
        policy_frame = pd.DataFrame(policy_rows, columns=["Path", "Value"])
        policy_frame.to_excel(writer, index=False, sheet_name='Clustering Policy')
        policy_worksheet = writer.sheets['Clustering Policy']
        for col_num, value in enumerate(policy_frame.columns.values):
            policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(policy_frame.iloc[row_num]):
                policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        policy_worksheet.set_column(0, 0, 42)
        policy_worksheet.set_column(1, 1, 90)
        policy_worksheet.freeze_panes(1, 0)
        policy_worksheet.autofilter(0, 0, len(policy_frame), len(policy_frame.columns) - 1)
        policy_worksheet.hide()

    if isinstance(class_layer_policy, dict) and class_layer_policy:
        class_policy_rows = _flatten_metadata(class_layer_policy)
        class_policy_frame = pd.DataFrame(class_policy_rows, columns=["Path", "Value"])
        class_policy_frame.to_excel(writer, index=False, sheet_name='Class Layer Policy')
        class_policy_worksheet = writer.sheets['Class Layer Policy']
        for col_num, value in enumerate(class_policy_frame.columns.values):
            class_policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(class_policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(class_policy_frame.iloc[row_num]):
                class_policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        class_policy_worksheet.set_column(0, 0, 42)
        class_policy_worksheet.set_column(1, 1, 90)
        class_policy_worksheet.freeze_panes(1, 0)
        class_policy_worksheet.autofilter(0, 0, len(class_policy_frame), len(class_policy_frame.columns) - 1)
        class_policy_worksheet.hide()

    if isinstance(descriptor_layer_policy, dict) and descriptor_layer_policy:
        descriptor_policy_rows = _flatten_metadata(descriptor_layer_policy)
        descriptor_policy_frame = pd.DataFrame(descriptor_policy_rows, columns=["Path", "Value"])
        descriptor_policy_frame.to_excel(writer, index=False, sheet_name='Descriptor Layer Policy')
        descriptor_policy_worksheet = writer.sheets['Descriptor Layer Policy']
        for col_num, value in enumerate(descriptor_policy_frame.columns.values):
            descriptor_policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(descriptor_policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(descriptor_policy_frame.iloc[row_num]):
                descriptor_policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        descriptor_policy_worksheet.set_column(0, 0, 42)
        descriptor_policy_worksheet.set_column(1, 1, 90)
        descriptor_policy_worksheet.freeze_panes(1, 0)
        descriptor_policy_worksheet.autofilter(0, 0, len(descriptor_policy_frame), len(descriptor_policy_frame.columns) - 1)
        descriptor_policy_worksheet.hide()

    if isinstance(topology_layer_policy, dict) and topology_layer_policy:
        topology_policy_rows = _flatten_metadata(topology_layer_policy)
        topology_policy_frame = pd.DataFrame(topology_policy_rows, columns=["Path", "Value"])
        topology_policy_frame.to_excel(writer, index=False, sheet_name='Topology Layer Policy')
        topology_policy_worksheet = writer.sheets['Topology Layer Policy']
        for col_num, value in enumerate(topology_policy_frame.columns.values):
            topology_policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(topology_policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(topology_policy_frame.iloc[row_num]):
                topology_policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        topology_policy_worksheet.set_column(0, 0, 42)
        topology_policy_worksheet.set_column(1, 1, 90)
        topology_policy_worksheet.freeze_panes(1, 0)
        topology_policy_worksheet.autofilter(0, 0, len(topology_policy_frame), len(topology_policy_frame.columns) - 1)
        topology_policy_worksheet.hide()

    if isinstance(ionisation_layer_policy, dict) and ionisation_layer_policy:
        ionisation_policy_rows = _flatten_metadata(ionisation_layer_policy)
        ionisation_policy_frame = pd.DataFrame(ionisation_policy_rows, columns=["Path", "Value"])
        ionisation_policy_frame.to_excel(writer, index=False, sheet_name='Ionisation Layer Policy')
        ionisation_policy_worksheet = writer.sheets['Ionisation Layer Policy']
        for col_num, value in enumerate(ionisation_policy_frame.columns.values):
            ionisation_policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(ionisation_policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(ionisation_policy_frame.iloc[row_num]):
                ionisation_policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        ionisation_policy_worksheet.set_column(0, 0, 42)
        ionisation_policy_worksheet.set_column(1, 1, 90)
        ionisation_policy_worksheet.freeze_panes(1, 0)
        ionisation_policy_worksheet.autofilter(0, 0, len(ionisation_policy_frame), len(ionisation_policy_frame.columns) - 1)
        ionisation_policy_worksheet.hide()

    if isinstance(toxicophore_layer_policy, dict) and toxicophore_layer_policy:
        toxicophore_policy_rows = _flatten_metadata(toxicophore_layer_policy)
        toxicophore_policy_frame = pd.DataFrame(toxicophore_policy_rows, columns=["Path", "Value"])
        toxicophore_policy_frame.to_excel(writer, index=False, sheet_name='Toxicophore Layer Policy')
        toxicophore_policy_worksheet = writer.sheets['Toxicophore Layer Policy']
        for col_num, value in enumerate(toxicophore_policy_frame.columns.values):
            toxicophore_policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(toxicophore_policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(toxicophore_policy_frame.iloc[row_num]):
                toxicophore_policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        toxicophore_policy_worksheet.set_column(0, 0, 42)
        toxicophore_policy_worksheet.set_column(1, 1, 90)
        toxicophore_policy_worksheet.freeze_panes(1, 0)
        toxicophore_policy_worksheet.autofilter(0, 0, len(toxicophore_policy_frame), len(toxicophore_policy_frame.columns) - 1)
        toxicophore_policy_worksheet.hide()

    if EXPORT_EXTERNAL_EVIDENCE and isinstance(chemont_layer_policy, dict) and chemont_layer_policy:
        chemont_policy_rows = _flatten_metadata(chemont_layer_policy)
        chemont_policy_frame = pd.DataFrame(chemont_policy_rows, columns=["Path", "Value"])
        chemont_policy_frame.to_excel(writer, index=False, sheet_name='ChemOnt Layer Policy')
        chemont_policy_worksheet = writer.sheets['ChemOnt Layer Policy']
        for col_num, value in enumerate(chemont_policy_frame.columns.values):
            chemont_policy_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(chemont_policy_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(chemont_policy_frame.iloc[row_num]):
                chemont_policy_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        chemont_policy_worksheet.set_column(0, 0, 42)
        chemont_policy_worksheet.set_column(1, 1, 90)
        chemont_policy_worksheet.freeze_panes(1, 0)
        chemont_policy_worksheet.autofilter(0, 0, len(chemont_policy_frame), len(chemont_policy_frame.columns) - 1)
        chemont_policy_worksheet.hide()

    if isinstance(distance_layer_metadata, dict) and distance_layer_metadata:
        distance_frame = pd.DataFrame(_flatten_metadata(distance_layer_metadata), columns=["Path", "Value"])
        distance_frame.to_excel(writer, index=False, sheet_name='Distance Layer Metadata')
        distance_worksheet = writer.sheets['Distance Layer Metadata']
        for col_num, value in enumerate(distance_frame.columns.values):
            distance_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(distance_frame)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(distance_frame.iloc[row_num]):
                distance_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        distance_worksheet.set_column(0, 0, 42)
        distance_worksheet.set_column(1, 1, 90)
        distance_worksheet.freeze_panes(1, 0)
        distance_worksheet.autofilter(0, 0, len(distance_frame), len(distance_frame.columns) - 1)
        distance_worksheet.hide()

    classification_columns = [
        'Resolved Name',
        'Original SMILES',
        'Standardized SMILES',
        'Corrected Chemical Class',
        'Taxonomy Path',
        'Taxonomy Path Steps',
        'Taxonomy Hierarchy',
        'Final Classification Record',
        'Classification Input Source',
        'Classification Status',
        'Classification Rule ID',
        'Classification Rule Version',
        'Manual Review Flag',
        'Manual Review Reason',
        'Detected Feature Profile'
    ]
    if all(column in df.columns for column in ['Corrected Chemical Class', 'Final Classification Record']):
        classification_for_export = df[[column for column in classification_columns if column in df.columns]].copy()
        classification_for_export.to_excel(writer, index=False, sheet_name='Structure Classification Audit')
        classification_worksheet = writer.sheets['Structure Classification Audit']
        for col_num, value in enumerate(classification_for_export.columns.values):
            classification_worksheet.write(0, col_num, value, header_format)
        for row_num in range(len(classification_for_export)):
            fmt = alt_row_format if row_num % 2 == 1 else cell_format
            for col_num, value in enumerate(classification_for_export.iloc[row_num]):
                if isinstance(value, (list, dict, tuple)) or type(value).__name__ == 'ndarray':
                    value = str(value)
                classification_worksheet.write(row_num + 1, col_num, '' if pd.isna(value) else value, fmt)
        for col_num, column in enumerate(classification_for_export.columns):
            classification_worksheet.set_column(col_num, col_num, min(80, _safe_column_width(classification_for_export, column, cap=80, floor=14)))
        classification_worksheet.freeze_panes(1, 0)
        classification_worksheet.autofilter(0, 0, len(classification_for_export), len(classification_for_export.columns) - 1)

    writer.close()
