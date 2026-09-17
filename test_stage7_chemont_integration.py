import unittest
from io import BytesIO
from unittest.mock import patch
from zipfile import ZipFile

import pandas as pd

from clustering_engine import process_clustering_and_adme
from excel_formatter import format_excel_output


class FakeSession:
    def close(self):
        pass


class Stage7ChemOntIntegrationTests(unittest.TestCase):
    def test_batch_chemont_enrichment_surfaces_hierarchy_without_changing_clustering(self):
        frame = pd.DataFrame({
            "Compound Name": ["member_a", "member_b", "member_c"],
            "Standardized SMILES": ["CCO", "CCCO", "CCCCO"],
            "InChIKey": ["AAAA", "BBBB", "CCCC"],
        })

        def fake_enrich(df, session):
            result = df.copy()
            result["ChemOnt_Kingdom"] = "Organic compounds"
            result["ChemOnt_Superclass"] = "Organic oxygen compounds"
            result["ChemOnt_Class"] = "Organooxygen compounds"
            result["ChemOnt_Subclass"] = "Ethers"
            result["ChemOnt_Direct_Parent"] = "Dialkyl ethers"
            result["ChemOnt_Molecular_Framework"] = "Aliphatic acyclic compounds"
            result["ChemOnt_Substituents"] = "Ether; Alkyl"
            result["ChemOnt_Classification_Version"] = "2.1"
            result["ChemOnt_Source"] = "ClassyFire ChemOnt API"
            result["ChemOnt_Retrieval_Status"] = "Retrieved"
            result.attrs["chemont_classifications"] = result[
                [
                    "InChIKey",
                    "Standardized SMILES",
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
            ].copy()
            return result

        with patch("database.SessionLocal", return_value=FakeSession()), \
            patch("chemont_adapter.chemont_enabled", return_value=True), \
            patch("chemont_adapter.fetch_chemont_by_inchikey", side_effect=lambda inchikey, db_session, rules=None, get=None, now=None: ({
                "classification_version": "2.1",
                "kingdom": {"name": "Organic compounds"},
                "superclass": {"name": "Organic oxygen compounds"},
                "class": {"name": "Organooxygen compounds"},
                "subclass": {"name": "Ethers"},
                "direct_parent": {"name": "Dialkyl ethers"},
                "molecular_framework": "Aliphatic acyclic compounds",
                "substituents": [{"name": "Ether"}, {"name": "Alkyl"}],
            }, "Retrieved")), \
            patch("chemont_adapter.enrich_with_chemont", side_effect=fake_enrich), \
            patch("evidence_ledger.attach_evidence_ledger", side_effect=lambda df: df), \
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

        self.assertIn("ChemOnt_Layer_Status", result.columns)
        self.assertIn("ChemOnt_Layer_Cohesion", result.columns)
        self.assertEqual(result["ChemOnt_Retrieval_Status"].iloc[0], "Retrieved")
        self.assertEqual(result["ChemOnt_Class"].iloc[0], "Organooxygen compounds")
        self.assertFalse(result.columns.duplicated().any())
        self.assertIn("chemont_classifications", result.attrs)

        workbook = BytesIO()
        format_excel_output(result, workbook, raw_tables=[])
        with ZipFile(BytesIO(workbook.getvalue())) as archive:
            workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
        self.assertIn("ChemOnt Classification", workbook_xml)
        self.assertIn("Cluster Evidence Matrix", workbook_xml)


if __name__ == "__main__":
    unittest.main()
