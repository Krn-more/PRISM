import unittest

from step6_export import build_step6_export_fields
from structural_evidence import calculate_structural_evidence


class Step6ExportTests(unittest.TestCase):
    def test_step6_export_fields_mirror_structural_evidence(self):
        evidence = calculate_structural_evidence("C1CCOC1")
        exported = build_step6_export_fields(evidence, "C1CCOC1")

        self.assertEqual(exported["Structural Evidence Status"], "Calculated")
        self.assertEqual(exported["Structural Evidence Version"], "1.0.0")
        self.assertEqual(exported["Toxicophore Profile"], "None")
        self.assertEqual(exported["All Structural Alerts"], "None")
        self.assertEqual(exported["Alerts"], "None")
        self.assertEqual(exported["Scaffold"], "C1CCOC1")

    def test_step6_export_fields_use_safe_fallbacks(self):
        exported = build_step6_export_fields(None, "")

        self.assertEqual(exported["Structural Evidence Status"], "")
        self.assertEqual(exported["Structural Evidence Version"], "")
        self.assertEqual(exported["Toxicophore Profile"], "None")
        self.assertEqual(exported["All Structural Alerts"], "None")
        self.assertEqual(exported["Alerts"], "None")
        self.assertEqual(exported["Scaffold"], "Not available")


if __name__ == "__main__":
    unittest.main()
