import unittest

import pandas as pd

from clustering_toxicophore_layer import (
    annotate_toxicophore_layer_evidence,
    build_toxicophore_distance_matrix,
    extract_toxicophore_signature,
    load_toxicophore_layer_rules,
)


class ToxicophoreLayerTests(unittest.TestCase):
    def test_rules_are_review_only(self):
        rules = load_toxicophore_layer_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertEqual(rules["rule_version"], "1.0.0-review-only")
        self.assertIn("Toxicophore_Profile", rules["toxicophore_columns"])

    def test_signature_uses_profile_and_alert_sets(self):
        row = pd.Series({
            "Toxicophore_Profile": "Aldehyde; Epoxide",
            "All_Structural_Alerts": "Aldehyde alert; Epoxide alert",
            "Alerts": "Aldehyde alert",
        })
        signature = extract_toxicophore_signature(row)
        self.assertEqual(signature["profile"], {"Aldehyde", "Epoxide"})
        self.assertEqual(signature["alerts"], {"Aldehyde alert", "Epoxide alert"})
        self.assertEqual(signature["first_alert"], "Aldehyde alert")

    def test_distance_matrix_prefers_shared_toxicophore_sets(self):
        df = pd.DataFrame([
            {"Toxicophore_Profile": "None", "All_Structural_Alerts": "None", "Alerts": "None"},
            {"Toxicophore_Profile": "None", "All_Structural_Alerts": "None", "Alerts": "None"},
            {"Toxicophore_Profile": "Aldehyde", "All_Structural_Alerts": "Aldehyde alert", "Alerts": "Aldehyde alert"},
        ])
        matrix, signatures = build_toxicophore_distance_matrix(df)
        self.assertEqual(len(signatures), 3)
        self.assertEqual(matrix[0, 1], 0.0)
        self.assertGreater(matrix[0, 2], matrix[0, 1])

    def test_cluster_annotation_reports_toxicophore_cohesion(self):
        df = pd.DataFrame([
            {"Cluster ID": "Cluster_1", "Toxicophore_Profile": "Aldehyde", "All_Structural_Alerts": "Aldehyde alert", "Alerts": "Aldehyde alert"},
            {"Cluster ID": "Cluster_1", "Toxicophore_Profile": "Aldehyde", "All_Structural_Alerts": "Aldehyde alert", "Alerts": "Aldehyde alert"},
        ])
        result = annotate_toxicophore_layer_evidence(df)
        self.assertEqual(result["Toxicophore_Layer_Status"].iloc[0], "Calculated")
        self.assertIn(result["Toxicophore_Layer_Cohesion"].iloc[0], {"Toxicophore coherent", "Toxicophore aligned"})
        self.assertIn("Mean pairwise toxicophore distance", result["Toxicophore_Layer_Reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
