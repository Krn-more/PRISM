import unittest
from io import BytesIO
from unittest.mock import patch
from zipfile import ZipFile

import pandas as pd

from clustering_engine import process_clustering_and_adme
from excel_formatter import format_excel_output


class Stage6ClusterDomainIntegrationTests(unittest.TestCase):
    def test_batch_cluster_and_domain_fields_are_present_and_exported(self):
        frame = pd.DataFrame({
            "Compound Name": ["member_a", "member_b", "member_c"],
            "Standardized SMILES": ["CCO", "CCCO", "CCCCO"],
        })

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
            result = process_clustering_and_adme(frame)

        self.assertFalse(result.columns.duplicated().any())
        for column in (
            "Cluster ID",
            "Cluster Size",
            "Cluster_Status",
            "Domain",
            "Decision",
            "Domain_Status",
            "Domain_Reasons",
            "Nearest_Neighbour_Tanimoto",
            "Cluster_Min_Pairwise_Tanimoto",
            "Cluster_Median_Pairwise_Tanimoto",
            "Cluster_Scaffold_Coverage",
            "Cluster_Representative_Scaffold",
            "Property_Outlier_Flags",
            "Ionisation_Consistency",
            "Toxicophore_Consistency",
            "Uncertainty_Summary",
            "XAI_Status",
            "Read_Across_Justification",
        ):
            self.assertIn(column, result.columns)

        self.assertTrue(result["Cluster ID"].astype(str).str.startswith("Cluster_").all())
        self.assertTrue(result["Cluster_Status"].astype(str).str.len().gt(0).all())
        self.assertTrue(result["Domain_Status"].astype(str).isin({"Inside", "Borderline", "Outside", "Not assessable"}).all())
        self.assertIn("cluster_evidence_matrix", result.attrs)
        self.assertIn("uncertainty_register", result.attrs)

        workbook = BytesIO()
        format_excel_output(result, workbook, raw_tables=[])
        with ZipFile(BytesIO(workbook.getvalue())) as archive:
            workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
        self.assertIn("Cluster Evidence Matrix", workbook_xml)
        self.assertIn("Uncertainty Register", workbook_xml)


if __name__ == "__main__":
    unittest.main()
