import unittest

import pandas as pd

from main import _validate_phase_b_export


def _frame(decision="PASS", domain="INSIDE DOMAIN"):
    return pd.DataFrame([{
        "Compound Name": "test",
        "Standardized SMILES": "C1CCOC1",
        "Identity Status": "Resolved",
        "Classification Status": "Calculated",
        "Cluster ID": "Cluster_1",
        "Cluster Size": 1,
        "Cluster_Status": "Valid Group",
        "Domain": domain,
        "Decision": decision,
        "Domain_Status": "Inside",
        "Structural_Evidence_Status": "Calculated",
    }])


class B2ValidationTests(unittest.TestCase):
    def test_required_fields_allow_valid_export_without_run_state(self):
        self.assertEqual(_validate_phase_b_export(_frame(), None), [])

    def test_unresolved_pass_is_blocked(self):
        frame = _frame()
        frame.loc[0, "Standardized SMILES"] = ""
        self.assertTrue(any("Unresolved rows" in error for error in _validate_phase_b_export(frame, None)))

    def test_retired_cramer_fields_are_blocked(self):
        frame = _frame()
        frame["Cramer Class"] = "Class II"
        self.assertTrue(any("Retired Cramer/TTC" in error for error in _validate_phase_b_export(frame, None)))


if __name__ == "__main__":
    unittest.main()
