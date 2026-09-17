import unittest

import pandas as pd

from clustering_chemont_layer import (
    annotate_chemont_layer_evidence,
    build_chemont_distance_matrix,
    extract_chemont_signature,
    load_chemont_layer_rules,
)


class ChemOntLayerTests(unittest.TestCase):
    def test_rules_are_review_only(self):
        rules = load_chemont_layer_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertEqual(rules["rule_version"], "1.0.0-review-only")
        self.assertIn("ChemOnt_Kingdom", rules["hierarchy_columns"])
        self.assertIn("Taxonomy Path", rules["taxonomy_fallback_columns"])

    def test_signature_prefers_chemont_then_taxonomy_then_class(self):
        chemont_row = pd.Series({
            "ChemOnt_Kingdom": "Organic compounds",
            "ChemOnt_Superclass": "Organic oxygen compounds",
            "ChemOnt_Class": "Organooxygen compounds",
            "ChemOnt_Subclass": "Ethers",
            "ChemOnt_Direct_Parent": "Dialkyl ethers",
            "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
            "ChemOnt_Substituents": "Ether; Alkyl",
            "Corrected Chemical Class": "Cyclic ether",
        })
        taxonomy_row = pd.Series({
            "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Ether | Subclass: Ether",
            "Taxonomy Path": "Oxygen-containing Compound → Ether → Ether",
            "Corrected Chemical Class": "Cyclic ether",
        })
        class_row = pd.Series({
            "Corrected Chemical Class": "Cyclic ether",
        })

        chemont_signature = extract_chemont_signature(chemont_row)
        taxonomy_signature = extract_chemont_signature(taxonomy_row)
        class_signature = extract_chemont_signature(class_row)

        self.assertEqual(chemont_signature["source"], "chemont")
        self.assertEqual(chemont_signature["chemont"]["ChemOnt_Substituents"], {"Alkyl", "Ether"})
        self.assertEqual(taxonomy_signature["source"], "taxonomy")
        self.assertEqual(taxonomy_signature["taxonomy_path"], "Oxygen-containing Compound → Ether → Ether")
        self.assertEqual(class_signature["source"], "class")

    def test_distance_matrix_prefers_shared_taxonomy(self):
        df = pd.DataFrame([
            {
                "ChemOnt_Kingdom": "Organic compounds",
                "ChemOnt_Superclass": "Organic oxygen compounds",
                "ChemOnt_Class": "Organooxygen compounds",
                "ChemOnt_Subclass": "Ethers",
                "ChemOnt_Direct_Parent": "Dialkyl ethers",
                "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
                "ChemOnt_Substituents": "Ether; Alkyl",
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Ether | Subclass: Ether",
                "Taxonomy Path": "Oxygen-containing Compound → Ether → Ether",
                "Corrected Chemical Class": "Cyclic ether",
            },
            {
                "ChemOnt_Kingdom": "Organic compounds",
                "ChemOnt_Superclass": "Organic oxygen compounds",
                "ChemOnt_Class": "Organooxygen compounds",
                "ChemOnt_Subclass": "Ethers",
                "ChemOnt_Direct_Parent": "Dialkyl ethers",
                "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
                "ChemOnt_Substituents": "Ether; Alkyl",
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Ether | Subclass: Ether",
                "Taxonomy Path": "Oxygen-containing Compound → Ether → Ether",
                "Corrected Chemical Class": "Cyclic ether",
            },
            {
                "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Alcohol | Subclass: Ring-substituted Alcohol",
                "Taxonomy Path": "Oxygen-containing Compound → Alcohol → Ring-substituted Alcohol",
                "Corrected Chemical Class": "Alcohol",
            },
        ])

        matrix, signatures = build_chemont_distance_matrix(df)
        self.assertEqual(len(signatures), 3)
        self.assertEqual(matrix[0, 1], 0.0)
        self.assertGreater(matrix[0, 2], matrix[0, 1])

    def test_cluster_annotation_reports_chemont_cohesion(self):
        df = pd.DataFrame([
            {
                "Cluster ID": "Cluster_1",
                "ChemOnt_Kingdom": "Organic compounds",
                "ChemOnt_Superclass": "Organic oxygen compounds",
                "ChemOnt_Class": "Organooxygen compounds",
                "ChemOnt_Subclass": "Ethers",
                "ChemOnt_Direct_Parent": "Dialkyl ethers",
                "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
                "ChemOnt_Substituents": "Ether; Alkyl",
                "Corrected Chemical Class": "Cyclic ether",
            },
            {
                "Cluster ID": "Cluster_1",
                "ChemOnt_Kingdom": "Organic compounds",
                "ChemOnt_Superclass": "Organic oxygen compounds",
                "ChemOnt_Class": "Organooxygen compounds",
                "ChemOnt_Subclass": "Ethers",
                "ChemOnt_Direct_Parent": "Dialkyl ethers",
                "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
                "ChemOnt_Substituents": "Ether; Alkyl",
                "Corrected Chemical Class": "Cyclic ether",
            },
        ])

        result = annotate_chemont_layer_evidence(df)
        self.assertEqual(result["ChemOnt_Layer_Status"].iloc[0], "Calculated")
        self.assertIn(result["ChemOnt_Layer_Cohesion"].iloc[0], {"Taxonomically coherent", "Taxonomically aligned"})
        self.assertIn("Mean pairwise taxonomy distance", result["ChemOnt_Layer_Reasons"].iloc[0])

    def test_singleton_cluster_is_not_assessable(self):
        df = pd.DataFrame([
            {
                "Cluster ID": "Cluster_1",
                "ChemOnt_Kingdom": "Organic compounds",
                "ChemOnt_Superclass": "Organic oxygen compounds",
                "ChemOnt_Class": "Organooxygen compounds",
                "ChemOnt_Subclass": "Ethers",
                "ChemOnt_Direct_Parent": "Dialkyl ethers",
                "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
                "ChemOnt_Substituents": "Ether; Alkyl",
                "Corrected Chemical Class": "Cyclic ether",
            }
        ])

        result = annotate_chemont_layer_evidence(df)
        self.assertEqual(result["ChemOnt_Layer_Status"].iloc[0], "Not assessable")
        self.assertIn("Singleton cluster", result["ChemOnt_Layer_Reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
