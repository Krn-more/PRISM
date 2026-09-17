import unittest
from io import BytesIO
from zipfile import ZipFile

import pandas as pd

from cluster_summary import build_cluster_summary
from clustering_engine import run_revalidated_clustering
from excel_formatter import format_excel_output


class RevalidatedClusteringTests(unittest.TestCase):
    def setUp(self):
        # Two source records share one structure.  They must remain visible in
        # the export, but may only contribute one structure to clustering.
        self.source = pd.DataFrame(
            {
                "Compound Name": ["D4", "D4 alias", "Ethyl benzoate", "Di-tert-butyl peroxide"],
                "Standardized SMILES": [
                    "C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1",
                    "C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1",
                    "CCOC(=O)c1ccccc1",
                    "CC(C)(C)OOC(C)(C)C",
                ],
                "Primary Functional Group": ["Siloxane", "Siloxane", "Ester", "Peroxide"],
                "Detected Feature Profile": ["Siloxane", "Siloxane", "Ester", "Peroxide"],
            }
        )

    def test_clustering_accepts_duplicate_source_index_labels(self):
        frame = pd.DataFrame(
            [
                {"Compound Name": "ethanol", "Standardized SMILES": "CCO"},
                {"Compound Name": "propanol", "Standardized SMILES": "CCCO"},
            ],
            index=[7, 7],
        )

        result = run_revalidated_clustering(frame)

        self.assertEqual(len(result), 2)
        self.assertTrue(result["Cluster ID"].astype(str).str.startswith("Cluster_").all())

    def test_duplicates_are_collapsed_and_incompatible_families_are_separated(self):
        result = run_revalidated_clustering(self.source)

        self.assertEqual(result.loc[0, "Subcluster ID"], result.loc[1, "Subcluster ID"])
        self.assertEqual(result.loc[0, "Cluster Size"], 2)
        self.assertEqual(result.loc[0, "Cluster Unique Structure Count"], 1)
        self.assertEqual(result.loc[0, "Duplicate Structure Count"], 1)
        self.assertNotEqual(result.loc[0, "Subcluster ID"], result.loc[2, "Subcluster ID"])
        self.assertIn("Ester", result.loc[2, "Compatibility Gate"])
        self.assertNotIn("Siloxane", result.loc[2, "Compatibility Gate"])
        self.assertIn("separate review domain", result.loc[3, "Compatibility Gate"].casefold())
        self.assertTrue(result["Decision"].eq("REVIEW").all())
        self.assertTrue(result["Cluster ID"].str.fullmatch(r"Cluster_\d{3}").all())
        self.assertTrue(result["Compatibility Group"].str.fullmatch(r"Family_\d{3}").all())

    def test_summary_is_auditable_and_exported(self):
        result = run_revalidated_clustering(self.source)
        summary = build_cluster_summary(result)

        d4_summary = summary.loc[summary["Cluster ID"] == result.loc[0, "Cluster ID"]].iloc[0]
        self.assertEqual(d4_summary["Source Row Count"], 2)
        self.assertEqual(d4_summary["Unique Structure Count"], 1)
        self.assertEqual(d4_summary["Duplicate or Alias Row Count"], 1)
        self.assertEqual(d4_summary["SME Review Status"], "Pending")
        self.assertIn("SME review is required", d4_summary["Assessment Boundary"])

        workbook = BytesIO()
        format_excel_output(result, workbook)
        with ZipFile(BytesIO(workbook.getvalue())) as archive:
            workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
        self.assertIn("Cluster Summary", workbook_xml)


if __name__ == "__main__":
    unittest.main()
