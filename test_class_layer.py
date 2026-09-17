import unittest

import pandas as pd

from clustering_class_layer import (
    annotate_cluster_class_evidence,
    build_class_distance_matrix,
    extract_class_signature,
    load_class_layer_rules,
)


class ChemicalClassLayerTests(unittest.TestCase):
    def test_rules_are_review_only(self):
        rules = load_class_layer_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertEqual(rules["rule_version"], "1.0.0-review-only")
        self.assertGreater(rules["soft_layer_weight"], 0)

    def test_extracts_signature_from_taxonomy_hierarchy(self):
        row = pd.Series({
            "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Alcohol | Subclass: Ring-substituted alcohol | Structural Type: Primary Alcohol",
            "Corrected Chemical Class": "Alcohol",
        })
        signature = extract_class_signature(row)
        self.assertEqual(signature["parent_class"], "Oxygen-containing Compound")
        self.assertEqual(signature["functional_group"], "Alcohol")
        self.assertEqual(signature["subclass"], "Ring-substituted alcohol")
        self.assertEqual(signature["structural_type"], "Primary Alcohol")

    def test_class_distance_prefers_shared_hierarchy(self):
        df = pd.DataFrame([
            {
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Alcohol | Subclass: Ring-substituted alcohol | Structural Type: Primary Alcohol",
                "Corrected Chemical Class": "Alcohol",
            },
            {
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Alcohol | Subclass: Ring-substituted alcohol | Structural Type: Primary Alcohol",
                "Corrected Chemical Class": "Alcohol",
            },
            {
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Ester | Subclass: Carboxylic ester | Structural Type: Cyclic Ester",
                "Corrected Chemical Class": "Ester",
            },
        ])
        rules = load_class_layer_rules()
        matrix, signatures = build_class_distance_matrix(df, rules)
        self.assertEqual(len(signatures), 3)
        self.assertEqual(matrix[0, 1], 0.0)
        self.assertGreater(matrix[0, 2], matrix[0, 1])

    def test_cluster_class_annotation_labels_cohesive_clusters(self):
        df = pd.DataFrame([
            {
                "Cluster ID": "Cluster_1",
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Alcohol | Subclass: Ring-substituted alcohol | Structural Type: Primary Alcohol",
                "Corrected Chemical Class": "Alcohol",
            },
            {
                "Cluster ID": "Cluster_1",
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Alcohol | Subclass: Ring-substituted alcohol | Structural Type: Primary Alcohol",
                "Corrected Chemical Class": "Alcohol",
            },
        ])
        result = annotate_cluster_class_evidence(df)
        self.assertEqual(result["Cluster_Class_Cohesion"].iloc[0], "Chemically coherent")
        self.assertIn("All members share class signature", result["Cluster_Class_Reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
