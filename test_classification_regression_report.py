import unittest
from pathlib import Path

from classification_regression_report import build_regression_report


class ClassificationRegressionReportTests(unittest.TestCase):
    def test_approved_corpus_has_no_changed_existing_decisions(self):
        report, summary = build_regression_report(
            Path("List of Inchikey and SMILES.xlsx"),
            Path("outputs/classification_corpus_assessed.xlsx"),
        )

        self.assertEqual(summary["source_rows"], 219)
        self.assertFalse(report.empty)
        self.assertTrue(report["Change Type"].eq("Added field").all())
        self.assertFalse(report["Field"].isin({
            "Corrected Chemical Class",
            "Classification Rule ID",
            "Classification Rule Version",
            "Taxonomy Path",
            "Manual Review Flag",
            "Manual Review Reason",
        }).any())


if __name__ == "__main__":
    unittest.main()
