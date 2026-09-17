import unittest

from descriptor_export import build_descriptor_export_fields
from structural_evidence import calculate_structural_evidence


class DescriptorExportTests(unittest.TestCase):
    def test_descriptor_export_fields_mirror_structural_evidence(self):
        evidence = calculate_structural_evidence("C1CCOC1")
        exported = build_descriptor_export_fields(evidence, {"molecular_weight": 72.06})

        self.assertEqual(exported["Structural Evidence Status"], "Calculated")
        self.assertEqual(exported["Structural Evidence Version"], "1.0.0")
        self.assertEqual(exported["MW"], 72.06)
        self.assertEqual(exported["LogP"], 0.8)
        self.assertEqual(exported["TPSA"], 9.23)
        self.assertEqual(exported["Min_EState"], evidence["Min_EState"])
        self.assertEqual(exported["Max_EState"], evidence["Max_EState"])
        self.assertEqual(exported["3D_Asphericity"], evidence["3D_Asphericity"])
        self.assertEqual(exported["3D_PMI1"], evidence["3D_PMI1"])
        self.assertEqual(exported["3D_PMI2"], evidence["3D_PMI2"])
        self.assertEqual(exported["3D_PMI3"], evidence["3D_PMI3"])
        self.assertEqual(exported["3D_RadiusOfGyration"], evidence["3D_RadiusOfGyration"])
        self.assertEqual(exported["Ring Count"], 1)
        self.assertEqual(exported["Aromatic Ring Count"], 0)
        self.assertEqual(exported["Fraction CSP3"], 1.0)
        self.assertEqual(exported["Ionisation Indicator"], "Neutral")
        self.assertEqual(exported["Toxicophore Profile"], "None")
        self.assertEqual(exported["All Structural Alerts"], "None")


if __name__ == "__main__":
    unittest.main()
