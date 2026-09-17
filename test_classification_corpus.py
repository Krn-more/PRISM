import unittest
import json
from pathlib import Path


class ClassificationCorpusTests(unittest.TestCase):
    def test_committed_corpus_baseline_matches_the_source_contract(self):
        source = Path(__file__).resolve().parent / "List of Inchikey and SMILES.xlsx"
        fixture_dir = source.parent / "fixtures"
        summary = json.loads((fixture_dir / "classification_corpus_baseline.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["row_count"], 219)
        self.assertEqual(summary["classified_rows"], 190)
        self.assertEqual(summary["manual_review_rows"], 29)
        self.assertEqual(summary["missing_supplied_inchikey_rows"], 52)
        self.assertEqual(summary["derived_inchikey_candidate_rows"], 23)
        # Version 1.5.1 resolves every valid corpus structure with a
        # structure-specific rule; the remaining 29 rows lack usable structure.
        self.assertEqual(summary["rule_review_queue_rows"], 0)
        for filename in summary["outputs"].values():
            self.assertTrue((fixture_dir / filename).exists())


if __name__ == "__main__":
    unittest.main()
