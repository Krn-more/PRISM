"""Regression tests for the one-authority Single Assessment classifier view."""

import inspect
import unittest
from unittest.mock import patch

import main
from structure_classification import classify_structure


class SingleAssessmentPrismContractTests(unittest.TestCase):
    class _DatabaseStub:
        def add(self, _record):
            pass

        def commit(self):
            pass

        def refresh(self, _record):
            pass

    def test_primary_functional_group_comes_from_prism_hierarchy(self):
        result = classify_structure("CCCCP(=O)(CCCC)CCCC")
        contract = main._prism_single_assessment_contract(result)

        self.assertEqual(contract["chemical_class"], "Trialkylphosphine oxide")
        self.assertEqual(contract["primary_functional_group"], "Phosphine Oxide")
        self.assertIn("STC-PHO-003", contract["classification_rationale"])

    def test_peroxide_and_siloxane_contracts_match_prism(self):
        controls = [
            ("CC(C)(C)OOC(C)(C)C", "Dialkyl peroxide", "Peroxide"),
            ("C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1", "Macrocyclic siloxane", "Siloxane"),
        ]
        for smiles, expected_class, expected_group in controls:
            with self.subTest(smiles=smiles):
                result = classify_structure(smiles)
                contract = main._prism_single_assessment_contract(result)
                self.assertEqual(contract["chemical_class"], expected_class)
                self.assertEqual(contract["primary_functional_group"], expected_group)

    def test_secondary_groups_are_classifier_derived(self):
        result = classify_structure("CCOC(=O)c1ccc(C(=O)O)cc1")
        contract = main._prism_single_assessment_contract(result)

        self.assertEqual(contract["primary_functional_group"], "Carboxylic Acid")
        self.assertTrue(all(group["SMARTS_pattern"] == "PRISM classifier-derived feature" for group in contract["secondary_functional_groups"]))

    def test_single_assessment_makes_no_chemont_runtime_call(self):
        source = inspect.getsource(main.assess_chemical)
        self.assertNotIn("chemont_adapter", source)
        self.assertNotIn("_resolve_chemont_bundle", source)
        self.assertIn("_prism_single_assessment_contract", source)

    def test_assessment_response_uses_the_prism_contract(self):
        request = main.schemas.AssessmentRequest(
            standardized_smiles="CCCCP(=O)(CCCC)CCCC",
            original_smiles="CCCCP(=O)(CCCC)CCCC",
        )
        with patch.object(main, "_persist_structure_classification"), \
             patch.object(main.services, "generate_fingerprint", return_value=object()), \
             patch.object(main.services, "analog_search_and_score", return_value=[]), \
             patch.object(main.services, "determine_confidence", return_value="High"):
            response = main.assess_chemical(request, self._DatabaseStub())

        self.assertEqual(response.chemical_class, response.corrected_chemical_class)
        self.assertEqual(response.corrected_chemical_class, "Trialkylphosphine oxide")
        self.assertEqual(response.primary_functional_group, "Phosphine Oxide")
        self.assertEqual(response.chemont_classification, {})
        self.assertIn("STC-PHO-003", response.classification_rationale)

    def test_generic_name_returns_manual_review_instead_of_an_http_error(self):
        request = main.schemas.AssessmentRequest(
            compound_name="Branched alkane (C=09) related compound 01",
        )
        with patch.object(main, "_persist_structure_classification"), \
             patch.object(main.services, "generate_fingerprint", return_value=None), \
             patch.object(main.services, "analog_search_and_score", return_value=[]), \
             patch.object(main.services, "determine_confidence", return_value="Low"):
            response = main.assess_chemical(request, self._DatabaseStub())

        self.assertEqual(response.corrected_chemical_class, "Unclassified")
        self.assertTrue(response.manual_review_flag)
        self.assertEqual(response.taxonomy_hierarchy["Parent Class"], "Unclassified")


if __name__ == "__main__":
    unittest.main()
