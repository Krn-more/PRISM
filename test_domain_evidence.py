import unittest
import json
import os
import tempfile

import pandas as pd

from domain_evidence import add_domain_evidence, load_domain_rules


class DomainEvidenceTests(unittest.TestCase):
    def test_rules_are_versioned_and_review_only(self):
        rules = load_domain_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertIn("rule_version", rules)
        self.assertEqual(rules["minimum_nearest_neighbour_similarity"], 0.6)

    def test_unsupported_configured_outlier_method_is_rejected(self):
        rules = load_domain_rules()
        rules["descriptor_outlier_method"] = "unsupported"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", encoding="utf-8", delete=False) as handle:
            json.dump(rules, handle)
            rules_path = handle.name
        try:
            with self.assertRaisesRegex(ValueError, "iqr_1.5"):
                load_domain_rules(rules_path)
        finally:
            os.unlink(rules_path)

    def test_coherent_members_are_inside_without_changing_legacy_fields(self):
        source = pd.DataFrame({
            "Cluster ID": ["Cluster_1", "Cluster_1"],
            "Standardized SMILES": ["CCOc1ccccc1", "CCCOc1ccccc1"],
            "Scaffold": ["c1ccccc1", "c1ccccc1"],
            "MW": [122.17, 136.19], "LogP": [2.0, 2.4], "TPSA": [9.23, 9.23],
            "HBD": [0, 0], "HBA": [1, 1], "Rotatable_Bonds": [2, 3],
            "Formal_Charge": [0, 0], "Fraction_CSP3": [0.29, 0.38],
            "Ionisation_Indicator": ["Neutral", "Neutral"],
            "Toxicophore_Profile": ["None", "None"],
            "Domain": ["INSIDE DOMAIN", "INSIDE DOMAIN"], "Decision": ["PASS", "PASS"],
        })
        result = add_domain_evidence(source)
        self.assertEqual(result["Domain_Status"].tolist(), ["Inside", "Inside"])
        self.assertEqual(result["Domain"].tolist(), source["Domain"].tolist())
        self.assertEqual(result["Decision"].tolist(), source["Decision"].tolist())
        self.assertTrue((result["Nearest_Neighbour_Tanimoto"] >= 0.6).all())
        self.assertIn("Not assessed", result["Property_Outlier_Flags"].iloc[0])

    def test_divergent_member_is_borderline_with_reasons(self):
        source = pd.DataFrame({
            "Cluster ID": ["Cluster_1", "Cluster_1"],
            "Standardized SMILES": ["CCO", "O=Cc1ccc(C=C)cc1"],
            "Scaffold": ["No_Scaffold", "c1ccccc1"],
            "MW": [46.07, 132.16], "LogP": [-0.3, 1.4], "TPSA": [20.23, 17.07],
            "HBD": [1, 0], "HBA": [1, 1], "Rotatable_Bonds": [0, 2],
            "Formal_Charge": [0, 0], "Fraction_CSP3": [1.0, 0.11],
            "Ionisation_Indicator": ["Neutral", "Neutral"],
            "Toxicophore_Profile": ["None", "Aldehyde"],
        })
        result = add_domain_evidence(source)
        self.assertTrue((result["Domain_Status"] == "Borderline").all())
        self.assertIn("below review threshold", result["Domain_Reasons"].iloc[1])
        self.assertIn("Toxicophore profile has no unique", result["Domain_Reasons"].iloc[1])

    def test_missing_structure_is_not_assessable(self):
        source = pd.DataFrame({"Cluster ID": ["Cluster_1"], "Standardized SMILES": [""], "Scaffold": ["Invalid"]})
        result = add_domain_evidence(source)
        self.assertEqual(result["Domain_Status"].iloc[0], "Not assessable")
        self.assertIn("missing or invalid", result["Domain_Reasons"].iloc[0])

    def test_tied_member_evidence_is_order_independent(self):
        source = pd.DataFrame({
            "Cluster ID": ["Cluster_1", "Cluster_1"],
            "Standardized SMILES": ["CCOc1ccccc1", "CCCOc1ccccc1"],
            "Scaffold": ["c1ccccc1", "c1ccccc1"],
            "Ionisation_Indicator": ["Neutral", "Likely basic"],
            "Toxicophore_Profile": ["None", "None"],
        })
        first = add_domain_evidence(source).set_index("Standardized SMILES")
        reversed_result = add_domain_evidence(source.iloc[::-1].reset_index(drop=True)).set_index("Standardized SMILES")
        columns = ["Domain_Status", "Domain_Reasons", "Ionisation_Consistency", "Toxicophore_Consistency"]
        self.assertTrue(first[columns].sort_index().equals(reversed_result[columns].sort_index()))
        self.assertTrue((first["Ionisation_Consistency"] == "No unique cluster majority").all())
        self.assertTrue(first["Domain_Reasons"].str.contains("no unique cluster majority").all())

    def test_surrogate_identity_is_borderline_and_unresolved_is_not_assessable(self):
        source = pd.DataFrame({
            "Cluster ID": ["Cluster_1", "Cluster_1"],
            "Standardized SMILES": ["CCOc1ccccc1", "CCCOc1ccccc1"],
            "Scaffold": ["c1ccccc1", "c1ccccc1"],
            "Ionisation_Indicator": ["Neutral", "Neutral"],
            "Toxicophore_Profile": ["None", "None"],
            "Identity Status": ["Surrogate SMILES (Homologue)", "Tier 3: Class Fallback Required (No Structure)"],
        })
        result = add_domain_evidence(source)
        self.assertEqual(result["Domain_Status"].tolist(), ["Borderline", "Not assessable"])
        self.assertIn("Identity status requires review", result["Domain_Reasons"].iloc[0])
        self.assertIn("Identity status prevents assessment", result["Domain_Reasons"].iloc[1])


if __name__ == "__main__":
    unittest.main()
