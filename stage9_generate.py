from pathlib import Path
from unittest.mock import patch
import os
import pandas as pd
from clustering_engine import process_clustering_and_adme
from excel_formatter import format_excel_output

fixture = Path('phase_a_D61148-1 C33511 Draft 1 (2).xlsx')
os.environ['DETERMINISTIC_LOCAL_MODE'] = 'true'
os.environ['PHASE_B_STAGE_DEADLINE_SECONDS'] = '30'
source = pd.read_excel(fixture, sheet_name='High Confidence')
frames = [source]
for sheet in ('Medium-Low Analogs', 'Unresolved (Review)'):
    frames.append(pd.read_excel(fixture, sheet_name=sheet))
source = pd.concat(frames, ignore_index=True).fillna('')
with patch('evidence_ledger.attach_evidence_ledger', side_effect=lambda df: df), \
     patch('epa_ctx_adapter.enrich_with_epa_ctx', side_effect=lambda df, session: df), \
     patch('external_api.fetch_chembl_for_dataframe', side_effect=lambda df, smiles_col: df), \
     patch('graph_exporter.generate_cypher_queries', return_value=None), \
     patch('xai_engine.generate_justifications', side_effect=lambda df, cluster_col='Cluster ID': df):
    result = process_clustering_and_adme(source)
out = Path('phase_a_stage9_final.xlsx')
format_excel_output(result, out)
print({'output': str(out), 'rows': len(result), 'structured': int(result['Standardized SMILES'].astype(str).str.strip().ne('').sum()), 'unresolved': int(result['Identity Status'].astype(str).str.contains('unresolved|review', case=False, regex=True).sum()), 'sheets': pd.ExcelFile(out).sheet_names})
