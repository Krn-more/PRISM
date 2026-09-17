import csv
import hashlib
import unittest
from pathlib import Path


MANIFEST = Path("fixtures") / "stage6_parallel_benchmark_manifest.csv"
REVIEW_LOG = Path("fixtures") / "stage6_parallel_reviewer_log.csv"
HASH_MANIFEST = Path("fixtures") / "stage6_parallel_benchmark_manifest.sha256"


class Stage6ParallelBenchmarkManifestTests(unittest.TestCase):
    def test_manifest_hash_is_frozen(self):
        expected_hash = HASH_MANIFEST.read_text(encoding="utf-8").split()[0]
        self.assertEqual(hashlib.sha256(MANIFEST.read_bytes()).hexdigest().upper(), expected_hash)

    def test_manifest_has_twenty_unique_review_cases(self):
        with open(MANIFEST, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 20)
        self.assertEqual(len({row["Case_ID"] for row in rows}), 20)
        self.assertTrue(all(row["Expected_Grounded_Eligibility"] == "Eligible" for row in rows))
        focus = " ".join(" ".join((row["Cluster_Profile"], row["Primary_Evidence_Variation"], row["Reviewer_Focus"])) for row in rows).casefold()
        for required in ("endpoint", "exposure", "domain", "identity", "property", "reactivity"):
            self.assertIn(required, focus)

    def test_reviewer_log_has_all_required_traceability_columns(self):
        with open(REVIEW_LOG, newline="", encoding="utf-8") as handle:
            columns = next(csv.reader(handle))
        for required in ("Case_ID", "Grounded_XAI_Status", "Matrix_Hash", "Prompt_Version", "Matrix_Version", "Deployment", "Response_Timestamp", "Claim_Evidence_Links_Valid", "Unsupported_Claims", "Rendered_Length", "Reviewer_Disposition"):
            self.assertIn(required, columns)


if __name__ == "__main__":
    unittest.main()
