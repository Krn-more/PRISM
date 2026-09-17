import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
from zipfile import ZipFile

import pandas as pd
from rdkit import Chem

from clustering_engine import process_clustering_and_adme, run_consensus_clustering
from excel_formatter import format_excel_output


FIXTURE_PATH = Path("fixtures") / "stage1_to_3_representative_benchmark.csv"


class EndToEndSmokeTests(unittest.TestCase):
    def test_consensus_copyback_accepts_excel_string_dtype_enrichment_columns(self):
        source = pd.DataFrame({
            "Standardized SMILES": pd.Series(["CCO", "not-a-smiles"], dtype="string"),
            "Cluster ID": pd.Series(["", ""], dtype="string"),
            "Cluster Size": pd.Series(["", ""], dtype="string"),
        })

        result = run_consensus_clustering(source, "Standardized SMILES")

        self.assertEqual(result.loc[1, "Cluster ID"], "Not assessable")
        self.assertNotEqual(result.loc[0, "Cluster ID"], "Not assessable")
        self.assertFalse(isinstance(result["Cluster Size"].dtype, pd.StringDtype))

    def test_frozen_fixture_runs_through_batch_pipeline_and_workbook_export(self):
        source = pd.read_csv(FIXTURE_PATH, keep_default_na=False)
        source["Domain"] = "LEGACY DOMAIN"
        source["Decision"] = "LEGACY DECISION"

        with patch("evidence_ledger.attach_evidence_ledger", side_effect=lambda df: df), \
            patch("epa_ctx_adapter.enrich_with_epa_ctx", side_effect=lambda df, session: df), \
            patch("external_api.fetch_chembl_for_dataframe", side_effect=lambda df, smiles_col: df), \
            patch("graph_exporter.generate_cypher_queries", return_value=None), \
            patch(
                "xai_engine.generate_justifications",
                side_effect=lambda df, cluster_col="Cluster ID": df.assign(
                    XAI_Status="Unavailable",
                    Read_Across_Justification="AI read-across summary unavailable.",
                ),
            ):
            result = process_clustering_and_adme(source)

        self.assertEqual(len(result), len(source))
        self.assertIn("Cluster ID", result.columns)
        self.assertIn("Domain_Status", result.columns)
        self.assertIn("Uncertainty_Summary", result.columns)
        self.assertIn("ChemOnt_Retrieval_Status", result.columns)
        self.assertFalse(result.columns.duplicated().any())
        # Invalid structures are intentionally excluded from clustering.  The
        # fixture includes two malformed-SMILES rows to verify that the batch
        # process marks them as not assessable instead of assigning analogues.
        valid_structure = source["Standardized SMILES"].map(
            lambda smiles: bool(str(smiles).strip()) and Chem.MolFromSmiles(str(smiles)) is not None
        )
        self.assertTrue(
            result.loc[valid_structure, "Cluster ID"].astype(str).str.startswith("Cluster_").all()
        )
        self.assertTrue(
            result.loc[~valid_structure, "Cluster ID"].eq("Not assessable").all()
        )

        workbook = BytesIO()
        format_excel_output(result, workbook, raw_tables=[source])
        sheet_names = pd.ExcelFile(BytesIO(workbook.getvalue())).sheet_names
        self.assertEqual(
            sheet_names[:7],
            [
                "SME Review Guide",
                "Dashboard",
                "Workbook Contract",
                "Raw Extraction (Pre-Dedup)",
                "Summary",
                "Cluster Summary",
                "Detailed Analysis",
            ],
        )
        contract_headers = pd.read_excel(BytesIO(workbook.getvalue()), sheet_name="Workbook Contract", nrows=0).columns.tolist()
        self.assertIn("Column Headers - Detailed Analysis sheet", contract_headers)
        with ZipFile(BytesIO(workbook.getvalue())) as archive:
            workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")

        for sheet_name in (
            "Raw Extraction (Pre-Dedup)",
            "Cluster Summary",
            "Workbook Contract",
            "Uncertainty Register",
            "Cluster Evidence Matrix",
        ):
            self.assertIn(sheet_name, workbook_xml)

        self.assertNotIn("ChemOnt Classification", workbook_xml)
        self.assertNotIn("Endpoint Evidence Ledger", workbook_xml)


if __name__ == "__main__":
    unittest.main()
