import hashlib
import unittest
from pathlib import Path

import pandas as pd

from clustering_engine import extract_scaffolds, flag_alerts
from domain_evidence import add_domain_evidence
from structural_evidence import add_structural_evidence
from uncertainty_register import add_uncertainty_register


FIXTURE_DIR = Path("fixtures")
FIXTURE_PATH = FIXTURE_DIR / "stage1_to_3_representative_benchmark.csv"
MANIFEST_PATH = FIXTURE_DIR / "stage1_to_3_representative_benchmark.sha256"


def run_benchmark() -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_csv(FIXTURE_PATH, keep_default_na=False)
    # Legacy fields are deliberately present to assert review-only invariance.
    source["Domain"] = "LEGACY DOMAIN"
    source["Decision"] = "LEGACY DECISION"
    result = add_structural_evidence(source)
    result = flag_alerts(result)
    result = extract_scaffolds(result)
    result = add_domain_evidence(result)
    result, register = add_uncertainty_register(result)
    return result.set_index("Fixture_ID"), register


class RepresentativeBenchmarkTests(unittest.TestCase):
    def test_fixture_hash_is_frozen(self):
        expected_hash = MANIFEST_PATH.read_text(encoding="utf-8").split()[0]
        actual_hash = hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest().upper()
        self.assertEqual(actual_hash, expected_hash)

    def test_stage1_to_3_expected_review_outcomes_and_legacy_invariance(self):
        result, register = run_benchmark()
        self.assertEqual(len(result), 13)
        self.assertEqual(result["Cluster ID"].to_dict(), pd.read_csv(FIXTURE_PATH, keep_default_na=False).set_index("Fixture_ID")["Cluster ID"].to_dict())
        self.assertTrue((result["Domain"] == "LEGACY DOMAIN").all())
        self.assertTrue((result["Decision"] == "LEGACY DECISION").all())

        self.assertEqual(result.at["B07", "Ionisation_Indicator"], "Likely acidic")
        self.assertEqual(result.at["B08", "Ionisation_Indicator"], "Likely basic")
        self.assertEqual(result.at["B09", "Ionisation_Indicator"], "Amphoteric")
        self.assertEqual(result.at["B10", "Ionisation_Indicator"], "Permanently charged")
        self.assertEqual(result.at["B12", "Structural_Evidence_Status"], "Missing structure")
        self.assertEqual(result.at["B13", "Structural_Evidence_Status"], "Invalid structure")

        self.assertIn("Rotatable_Bonds", result.at["B04", "Property_Outlier_Flags"])
        self.assertEqual(result.at["B06", "Toxicophore_Profile"], "Aldehyde")
        self.assertEqual(result.at["B06", "Domain_Status"], "Borderline")
        for fixture_id in ("B07", "B08", "B09", "B10", "B11"):
            self.assertEqual(result.at[fixture_id, "Domain_Status"], "Borderline")
        for fixture_id in ("B12", "B13"):
            self.assertEqual(result.at[fixture_id, "Domain_Status"], "Not assessable")

        compound_rows = register[register["Scope"] == "Compound"]
        b11_identity = compound_rows[compound_rows["Scope_ID"].str.contains("Surrogate identity phenyl ether", regex=False) & (compound_rows["Category"] == "Identity")].iloc[0]
        b12_identity = compound_rows[compound_rows["Scope_ID"].str.contains("Unresolved identity", regex=False) & (compound_rows["Category"] == "Identity")].iloc[0]
        self.assertEqual(b11_identity["Level"], "Medium")
        self.assertEqual(b12_identity["Level"], "High")
        self.assertIn("Document member-specific differences", compound_rows[(compound_rows["Scope_ID"].str.contains("Acetic acid", regex=False)) & (compound_rows["Category"] == "Domain")].iloc[0]["Required_Action"])

    def test_review_outcomes_are_stable_when_fixture_rows_are_reordered(self):
        first, _ = run_benchmark()
        source = pd.read_csv(FIXTURE_PATH, keep_default_na=False).iloc[::-1].reset_index(drop=True)
        source["Domain"] = "LEGACY DOMAIN"
        source["Decision"] = "LEGACY DECISION"
        second = add_structural_evidence(source)
        second = flag_alerts(second)
        second = extract_scaffolds(second)
        second = add_domain_evidence(second).set_index("Fixture_ID")
        columns = ["Ionisation_Indicator", "Structural_Evidence_Status", "Domain_Status", "Property_Outlier_Flags", "Domain_Reasons"]
        self.assertTrue(first[columns].sort_index().equals(second[columns].sort_index()))


if __name__ == "__main__":
    unittest.main()
