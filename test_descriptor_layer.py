import unittest

import pandas as pd

from clustering_descriptor_layer import (
    annotate_descriptor_layer_evidence,
    build_descriptor_distance_matrix,
    load_descriptor_layer_rules,
)


class DescriptorLayerTests(unittest.TestCase):
    def test_rules_are_review_only(self):
        rules = load_descriptor_layer_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertEqual(rules["rule_version"], "1.0.0-review-only")
        self.assertIn("Ionisation_Indicator", rules["descriptor_columns"])

    def test_distance_matrix_detects_descriptor_divergence(self):
        df = pd.DataFrame([
            {
                "HBD": 0,
                "HBA": 1,
                "Rotatable_Bonds": 0,
                "Formal_Charge": 0,
                "Ring_Count": 1,
                "Aromatic_Ring_Count": 0,
                "Fraction_CSP3": 1.0,
                "Heavy_Atom_Count": 5,
                "Molecular_Refractivity": 20.053,
                "Ionisation_Indicator": "Neutral",
            },
            {
                "HBD": 2,
                "HBA": 3,
                "Rotatable_Bonds": 5,
                "Formal_Charge": 1,
                "Ring_Count": 3,
                "Aromatic_Ring_Count": 1,
                "Fraction_CSP3": 0.15,
                "Heavy_Atom_Count": 12,
                "Molecular_Refractivity": 41.200,
                "Ionisation_Indicator": "Likely basic",
            },
        ])
        matrix, metadata = build_descriptor_distance_matrix(df)
        self.assertEqual(matrix.shape, (2, 2))
        self.assertEqual(matrix[0, 0], 0.0)
        self.assertGreater(matrix[0, 1], 0.0)
        self.assertEqual(metadata["status"], "Calculated")
        self.assertGreater(metadata["coverage_ratio"], 0.0)

    def test_cluster_annotation_marks_mixed_descriptor_cohesion(self):
        df = pd.DataFrame([
            {
                "Cluster ID": "Cluster_1",
                "HBD": 0,
                "HBA": 1,
                "Rotatable_Bonds": 0,
                "Formal_Charge": 0,
                "Ring_Count": 1,
                "Aromatic_Ring_Count": 0,
                "Fraction_CSP3": 1.0,
                "Heavy_Atom_Count": 5,
                "Molecular_Refractivity": 20.053,
                "Ionisation_Indicator": "Neutral",
            },
            {
                "Cluster ID": "Cluster_1",
                "HBD": 2,
                "HBA": 3,
                "Rotatable_Bonds": 5,
                "Formal_Charge": 1,
                "Ring_Count": 3,
                "Aromatic_Ring_Count": 1,
                "Fraction_CSP3": 0.15,
                "Heavy_Atom_Count": 12,
                "Molecular_Refractivity": 41.200,
                "Ionisation_Indicator": "Likely basic",
            },
        ])
        result = annotate_descriptor_layer_evidence(df)
        self.assertEqual(result["Descriptor_Layer_Status"].iloc[0], "Calculated")
        self.assertIn(result["Descriptor_Layer_Cohesion"].iloc[0], {"Chemically aligned", "Chemically mixed"})
        self.assertIn("Mean pairwise descriptor distance", result["Descriptor_Layer_Reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
