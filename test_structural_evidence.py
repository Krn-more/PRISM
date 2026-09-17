import os
import unittest
from unittest.mock import patch

import pandas as pd

from structural_evidence import calculate_structural_evidence
from structural_evidence import add_structural_evidence


class StructuralEvidenceTests(unittest.TestCase):
    def test_fast_deterministic_mode_skips_only_optional_3d_embedding(self):
        with patch.dict(os.environ, {"PRISM_ENABLE_3D_EVIDENCE": "false"}, clear=False):
            result = calculate_structural_evidence("CC(C)(C)OOC(C)(C)C")

        self.assertEqual(result["Structural_Evidence_Status"], "Calculated")
        self.assertEqual(result["3D_Evidence_Status"], "Skipped in fast deterministic mode")
        self.assertIsNone(result["3D_Asphericity"])
        self.assertGreater(result["MW"], 0)
        self.assertGreaterEqual(result["TPSA"], 0)

    def test_refresh_replaces_existing_evidence_columns_without_duplicates(self):
        source = pd.DataFrame({
            "Standardized SMILES": ["CCO"],
            "MW": [999.0],
            "3D_Evidence_Status": ["stale"],
        })

        result = add_structural_evidence(source)

        self.assertFalse(result.columns.duplicated().any())
        self.assertNotEqual(result.loc[0, "MW"], 999.0)
        self.assertEqual(result.loc[0, "3D_Evidence_Status"], "Skipped in fast deterministic mode")


if __name__ == "__main__":
    unittest.main()
