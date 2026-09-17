import unittest

import pandas as pd

from clustering_ionisation_layer import (
    annotate_ionisation_layer_evidence,
    build_ionisation_distance_matrix,
    extract_ionisation_signature,
    load_ionisation_layer_rules,
)


class IonisationLayerTests(unittest.TestCase):
    def test_rules_are_review_only(self):
        rules = load_ionisation_layer_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertEqual(rules["rule_version"], "1.0.0-review-only")
        self.assertIn("Ionisation_Indicator", rules["ionisation_columns"])

    def test_signature_uses_indicator_and_formal_charge(self):
        row = pd.Series({
            "Ionisation_Indicator": "Likely acidic",
            "Formal_Charge": -1,
        })
        signature = extract_ionisation_signature(row)
        self.assertEqual(signature["indicator"], "Likely acidic")
        self.assertEqual(signature["formal_charge"], -1)
        self.assertEqual(signature["charge_band"], "Negative")

    def test_distance_matrix_penalises_contrasting_ionisation_states(self):
        df = pd.DataFrame([
            {"Ionisation_Indicator": "Neutral", "Formal_Charge": 0},
            {"Ionisation_Indicator": "Neutral", "Formal_Charge": 0},
            {"Ionisation_Indicator": "Likely acidic", "Formal_Charge": -1},
        ])
        matrix, signatures = build_ionisation_distance_matrix(df)
        self.assertEqual(len(signatures), 3)
        self.assertEqual(matrix[0, 1], 0.0)
        self.assertGreater(matrix[0, 2], matrix[0, 1])
        self.assertGreaterEqual(matrix[0, 2], 0.5)

    def test_cluster_annotation_reports_ionisation_cohesion(self):
        df = pd.DataFrame([
            {"Cluster ID": "Cluster_1", "Ionisation_Indicator": "Neutral", "Formal_Charge": 0},
            {"Cluster ID": "Cluster_1", "Ionisation_Indicator": "Neutral", "Formal_Charge": 0},
        ])
        result = annotate_ionisation_layer_evidence(df)
        self.assertEqual(result["Ionisation_Layer_Status"].iloc[0], "Calculated")
        self.assertIn(result["Ionisation_Layer_Cohesion"].iloc[0], {"Ionisation coherent", "Ionisation aligned"})
        self.assertIn("Mean pairwise ionisation distance", result["Ionisation_Layer_Reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
