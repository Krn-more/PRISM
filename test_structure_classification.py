import unittest

import pandas as pd

from structure_classification import (
    apply_classification_contract,
    add_structure_classification,
    classification_export_fields,
    classification_secondary_functional_groups,
    classify_structure,
    load_ontology_compatibility_mappings,
    load_rules,
)


class StructureClassificationTests(unittest.TestCase):
    def test_rules_are_versioned_and_layered(self):
        rules = load_rules()
        self.assertEqual(rules["rule_version"], "1.5.1")
        self.assertIn("primary_class_vocabulary", rules)
        self.assertIn("subtype_vocabulary", rules)
        self.assertIn("Ketone", rules["primary_class_vocabulary"])
        self.assertIn("Cyclic ketal", rules["primary_class_vocabulary"])
        self.assertIn("Dialkyl peroxide", rules["primary_class_vocabulary"])
        self.assertIn("Macrocyclic siloxane", rules["primary_class_vocabulary"])
        self.assertIn("1,2-Oxasilinane", rules["primary_class_vocabulary"])

    def test_approved_ontology_compatibility_mappings(self):
        mappings = load_ontology_compatibility_mappings()
        self.assertEqual(mappings["ontology_mapping_version"], "1.0.0")

        phosphine_oxide = classify_structure("CCCCP(=O)(CCCC)CCCC")["ontology_compatibility"]
        self.assertEqual(phosphine_oxide["mapping_status"], "Mapped")
        self.assertEqual(phosphine_oxide["superclass"], "Organophosphorus compounds")
        self.assertEqual(phosphine_oxide["direct_parent"], "Organophosphine oxides")

        hydroperoxide = classify_structure("CC(C)(C)OO")["ontology_compatibility"]
        self.assertEqual(hydroperoxide["class"], "Organic hydroperoxides")

        peroxide = classify_structure("CC(C)(C)OOC(C)(C)C")
        self.assertEqual(peroxide["corrected_chemical_class"], "Dialkyl peroxide")
        self.assertEqual(peroxide["classification_rule_id"], "STC-PER-001")
        self.assertEqual(peroxide["taxonomy_hierarchy"]["Functional Group"], "Peroxide")
        self.assertEqual(peroxide["ontology_compatibility"]["class"], "Organic oxides")
        self.assertEqual(peroxide["ontology_compatibility"]["direct_parent"], "Dialkyl peroxides")

    def test_workbook_classification_flattens_ontology_mapping_columns(self):
        classified = add_structure_classification(pd.DataFrame({
            "Standardized SMILES": ["CC(C)(C)OOC(C)(C)C"],
            "Original SMILES": ["CC(C)(C)OOC(C)(C)C"],
        }))

        self.assertNotIn("Ontology Compatibility", classified.columns)
        self.assertIn(
            classified.loc[0, "Ontology Mapping Status"],
            {"Mapped", "Exact local taxonomy match"},
        )
        self.assertEqual(classified.loc[0, "Ontology Class"], "Organic oxides")
        self.assertEqual(classified.loc[0, "Ontology Direct Parent"], "Dialkyl peroxides")

    def test_core_primary_classifications(self):
        self.assertEqual(classify_structure("CC(=O)C")["corrected_chemical_class"], "Ketone")
        self.assertEqual(classify_structure("O=C1CCCCC1")["corrected_chemical_class"], "Cyclic ketone")
        self.assertEqual(classify_structure("CCOP(=O)(OCC)OCC")["corrected_chemical_class"], "Trialkyl phosphate")
        self.assertEqual(classify_structure("CCCCP(=O)(CCCC)CCCC")["corrected_chemical_class"], "Trialkylphosphine oxide")
        self.assertEqual(classify_structure("C1CCOC1")["corrected_chemical_class"], "Cyclic ether")
        self.assertEqual(classify_structure("O=C1CCCCCN1")["corrected_chemical_class"], "Lactam")
        self.assertEqual(classify_structure("O=C1OCCCC1")["corrected_chemical_class"], "Lactone")
        self.assertEqual(classify_structure("COCCOCCOCCOCC")["corrected_chemical_class"], "Polyether")
        self.assertEqual(classify_structure("CCOC1(CCCCC1)OCC")["corrected_chemical_class"], "Cyclic ketal")

    def test_fallback_and_manual_review_behavior(self):
        unresolved = classify_structure("", "")
        original = classify_structure("", "CCO")

        self.assertEqual(unresolved["corrected_chemical_class"], "Unclassified")
        self.assertTrue(unresolved["manual_review_flag"])
        self.assertEqual(unresolved["classification_input_source"], "Unavailable")
        self.assertEqual(unresolved["taxonomy_hierarchy"]["Parent Class"], "Unclassified")
        self.assertEqual(unresolved["taxonomy_hierarchy"]["Structural Type"], "Not assessable")
        self.assertEqual(original["corrected_chemical_class"], "Alcohol")
        self.assertEqual(original["classification_input_source"], "Original SMILES")
        self.assertFalse(original["manual_review_flag"])

    def test_feature_profile_retains_supporting_information(self):
        result = classify_structure("CCOC1(CCCCC1)OCC")
        self.assertIn("Ketal", result["detected_feature_profile"])
        self.assertIn("Cyclic ketal", result["detected_feature_profile"])
        self.assertIn("Ether", result["detected_feature_profile"])

    def test_taxonomy_path_is_built_from_broad_to_specific(self):
        result = classify_structure("CCOC1(CCCCC1)OCC")
        self.assertEqual(
            result["taxonomy_path"],
            "Oxygen-containing Compound → Acetal → Ketal → Cyclic Ketal",
        )
        self.assertEqual(
            result["taxonomy_path_steps"],
            ["Oxygen-containing Compound", "Acetal", "Ketal", "Cyclic Ketal"],
        )
        self.assertIn("Taxonomy path:", result["final_classification_record"])

    def test_alcohol_taxonomy_uses_topology_and_substitution(self):
        cyclic = classify_structure("OC1CCCCC1")
        ring_substituted = classify_structure("C1CCCCC1CO")
        bicyclic = classify_structure("OC1CCCCC1CC2CCCCC2")

        self.assertEqual(
            cyclic["taxonomy_path"],
            "Oxygen-containing Compound → Alcohol → Monocyclic Alcohol → Secondary Alcohol",
        )
        self.assertEqual(
            cyclic["taxonomy_path_steps"],
            ["Oxygen-containing Compound", "Alcohol", "Monocyclic Alcohol", "Secondary Alcohol"],
        )
        self.assertEqual(cyclic["taxonomy_hierarchy"]["Subclass"], "Monocyclic Alcohol")
        self.assertEqual(cyclic["taxonomy_hierarchy"]["Structural Type"], "Secondary Alcohol")

        self.assertEqual(
            ring_substituted["taxonomy_path"],
            "Oxygen-containing Compound → Alcohol → Monocyclic Alcohol → Primary Alcohol",
        )
        self.assertEqual(
            ring_substituted["taxonomy_path_steps"],
            ["Oxygen-containing Compound", "Alcohol", "Monocyclic Alcohol", "Primary Alcohol"],
        )
        self.assertEqual(
            ring_substituted["taxonomy_hierarchy"]["Subclass"],
            "Monocyclic Alcohol",
        )
        self.assertEqual(
            ring_substituted["taxonomy_hierarchy"]["Structural Type"],
            "Primary Alcohol",
        )

        self.assertEqual(
            bicyclic["taxonomy_hierarchy"]["Subclass"],
            "Bicyclic Alcohol",
        )
        self.assertEqual(
            bicyclic["taxonomy_hierarchy"]["Structural Type"],
            "Secondary Alcohol",
        )
        self.assertIn("Bicyclic Alcohol", bicyclic["taxonomy_path"])

    def test_topology_profile_is_structured_and_family_agnostic(self):
        acyclic = classify_structure("CCO")
        aromatic = classify_structure("c1ccccc1")

        self.assertEqual(acyclic["topology_profile"]["Topology Class"], "Acyclic")
        self.assertEqual(acyclic["topology_profile"]["Ring System"], "Open-chain")
        self.assertEqual(acyclic["topology_profile"]["Ring Profile"], "0 ring(s), 0 aromatic ring(s)")

        self.assertEqual(aromatic["topology_profile"]["Topology Class"], "Monocyclic")
        self.assertIn("Aromatic", aromatic["topology_profile"]["Topology Modifiers"])
        self.assertEqual(aromatic["topology_profile"]["Ring Profile"], "1 ring(s), 1 aromatic ring(s)")

    def test_monoethyl_terephthalate_prioritizes_acid_and_retains_secondary_groups(self):
        result = classify_structure("CCOC(=O)c1ccc(C(=O)O)cc1")

        self.assertEqual(result["corrected_chemical_class"], "Carboxylic acid")
        self.assertEqual(
            result["taxonomy_path"],
            "Oxygen-containing Compound → Carboxylic Acid → Aromatic Carboxylic Acid → Monocyclic Aromatic Carboxylic Acid",
        )
        self.assertEqual(
            result["taxonomy_hierarchy"]["Structural Type"],
            "Monocyclic Aromatic Carboxylic Acid",
        )

    def test_structural_evidence_is_available_for_standard_structures(self):
        result = classify_structure("CCO")
        self.assertEqual(result["topology_profile"]["Topology Class"], "Acyclic")

    def test_trialkylphosphine_oxide_has_phosphorus_taxonomy(self):
        result = classify_structure("CCCCP(=O)(CCCC)CCCC")

        self.assertEqual(result["classification_rule_id"], "STC-PHO-003")
        self.assertEqual(
            result["taxonomy_path"],
            "Phosphorus-containing Compound → Phosphine Oxide → Trialkylphosphine Oxide",
        )
        self.assertEqual(result["taxonomy_hierarchy"]["Functional Group"], "Phosphine Oxide")
        self.assertEqual(result["taxonomy_hierarchy"]["Structural Type"], "Acyclic Phosphine Oxide")
        self.assertIn("Trialkylphosphine oxide", result["detected_feature_profile"])

    def test_cyclic_and_macrocyclic_siloxanes_are_refined_locally(self):
        cyclic = classify_structure("C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O1")
        macrocyclic = classify_structure("C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1")

        self.assertEqual(cyclic["corrected_chemical_class"], "Cyclic siloxane")
        self.assertEqual(cyclic["classification_rule_id"], "STC-SIL-003")
        self.assertEqual(cyclic["taxonomy_hierarchy"]["Subclass"], "Cyclic Siloxane")
        self.assertEqual(macrocyclic["corrected_chemical_class"], "Macrocyclic siloxane")
        self.assertEqual(macrocyclic["classification_rule_id"], "STC-SIL-004")
        self.assertEqual(macrocyclic["taxonomy_hierarchy"]["Structural Type"], "Macrocyclic Siloxane")
        self.assertEqual(macrocyclic["ontology_compatibility"]["mapping_status"], "Not mapped")

    def test_acetal_and_trialkyl_phosphate_rules_exclude_known_lookalikes(self):
        geminal_diol = classify_structure("CCCCCC(O)(O)CCCO")
        triaryl_phosphate = classify_structure("CC(C)(C)c1ccc(OP(=O)(Oc2ccc(C(C)(C)C)cc2C(C)(C)C)Oc2ccc(C(C)(C)C)cc2C(C)(C)C)c(C(C)(C)C)c1")
        trialkyl_phosphate = classify_structure("CCOP(=O)(OCC)OCC")

        self.assertNotEqual(geminal_diol["corrected_chemical_class"], "Acetal")
        self.assertEqual(geminal_diol["corrected_chemical_class"], "Alcohol")
        self.assertEqual(triaryl_phosphate["corrected_chemical_class"], "Phosphate ester")
        self.assertEqual(trialkyl_phosphate["corrected_chemical_class"], "Trialkyl phosphate")

    def test_distinct_organosilicon_fallbacks_have_strict_local_rules(self):
        oxasilolane = classify_structure("C1CC[SiH2]OC1")
        silyl_alkyl_peroxide = classify_structure("CC(C)c1ccc(C(C)(C)OO[Si](C)(C)C)cc1")
        dialkyl_peroxide = classify_structure("CC(C)(C)OOC(C)(C)C")

        self.assertEqual(oxasilolane["corrected_chemical_class"], "1,2-Oxasilinane")
        self.assertEqual(oxasilolane["classification_rule_id"], "STC-SIL-005")
        self.assertEqual(oxasilolane["taxonomy_hierarchy"]["Functional Group"], "Cyclic Organosilicon Ether")
        self.assertEqual(silyl_alkyl_peroxide["corrected_chemical_class"], "Silyl alkyl peroxide")
        self.assertEqual(silyl_alkyl_peroxide["classification_rule_id"], "STC-SIL-006")
        self.assertEqual(silyl_alkyl_peroxide["taxonomy_hierarchy"]["Functional Group"], "Organosilicon Peroxide")
        self.assertEqual(dialkyl_peroxide["corrected_chemical_class"], "Dialkyl peroxide")

    def test_hydrocarbon_guard_and_elemental_fallbacks(self):
        benzene = classify_structure("c1ccccc1")
        cyclohexane = classify_structure("C1CCCCC1")
        hydroperoxide = classify_structure("CC(C)(C)OO")
        sulfide = classify_structure("CCSCC")
        silane = classify_structure("C[Si](C)C")

        self.assertEqual(benzene["corrected_chemical_class"], "Hydrocarbon")
        self.assertEqual(benzene["taxonomy_hierarchy"]["Subclass"], "Aromatic Hydrocarbon")
        self.assertEqual(cyclohexane["corrected_chemical_class"], "Hydrocarbon")
        self.assertEqual(cyclohexane["taxonomy_hierarchy"]["Subclass"], "Monocyclic Hydrocarbon")
        self.assertEqual(hydroperoxide["corrected_chemical_class"], "Organic hydroperoxide")
        self.assertEqual(hydroperoxide["classification_rule_id"], "STC-HYP-001")
        self.assertEqual(hydroperoxide["taxonomy_hierarchy"]["Functional Group"], "Hydroperoxide")
        self.assertEqual(sulfide["corrected_chemical_class"], "Sulfur-containing compound — not further classified")
        self.assertEqual(silane["corrected_chemical_class"], "Silicon-containing compound — not further classified")
        self.assertEqual(silane["classification_scope"], "Broad provisional class")
        self.assertTrue(silane["classification_review_recommended"])
        self.assertFalse(silane["manual_review_flag"])
        self.assertIn("review recommended", silane["classification_status"])
        self.assertIn("curate a specific structural rule", silane["classification_suggested_action"])
        self.assertNotEqual(sulfide["corrected_chemical_class"], "Hydrocarbon")
        self.assertNotEqual(silane["corrected_chemical_class"], "Hydrocarbon")

    def test_nitro_group_is_not_reported_as_an_amine(self):
        nitrophenol = classify_structure("CC(C)(C)c1cc([N+](=O)[O-])cc(C(C)(C)C)c1O")

        self.assertEqual(nitrophenol["corrected_chemical_class"], "Phenol")
        self.assertNotIn("Amine", nitrophenol["detected_feature_profile"])

    def test_export_fields_preserve_the_deterministic_contract(self):
        classification = classify_structure("CC(C)(C)OOC(C)(C)C")
        fields = classification_export_fields(classification)

        self.assertEqual(fields["Corrected Chemical Class"], "Dialkyl peroxide")
        self.assertEqual(fields["Classification Standardization Version"], "RDKit canonical-isomeric-v1")
        self.assertEqual(fields["Classification Structure SMILES"], "CC(C)(C)OOC(C)(C)C")
        self.assertEqual(fields["Classification Scope"], "Specific rule")
        self.assertEqual(fields["Classification Suggested Action"], "No classification action required")
        self.assertEqual(fields["Primary Functional Group"], "Peroxide")
        self.assertEqual(fields["Direct Parent"], "Organic Peroxide")
        self.assertEqual(fields["Matched Categories"], "Oxygen-containing Compound → Peroxide → Organic Peroxide → Dialkyl Peroxide")
        self.assertEqual(fields["Taxonomy Dictionary Version"], "1.0.0")
        self.assertEqual(fields["Ontology Mapping Status"], "Mapped")

    def test_secondary_features_require_non_overlapping_evidence(self):
        aminophenol = classify_structure("Nc1ccc(O)cc1")
        secondary, basis = classification_secondary_functional_groups(aminophenol)

        self.assertEqual(aminophenol["corrected_chemical_class"], "Phenol")
        self.assertEqual(secondary, "Amine")
        self.assertEqual(basis, "Non-overlapping detected feature atoms")

    def test_batch_contract_refresh_replaces_stale_classification_values(self):
        source = pd.DataFrame([{
            "Standardized SMILES": "CC(C)(C)OOC(C)(C)C",
            "Original SMILES": "CC(C)(C)OOC(C)(C)C",
            "Corrected Chemical Class": "Stale value",
        }])
        result = apply_classification_contract(source)

        self.assertEqual(result.loc[0, "Corrected Chemical Class"], "Dialkyl peroxide")
        self.assertEqual(result.loc[0, "Classification Scope"], "Specific rule")
        self.assertEqual(result.loc[0, "Classification Rule Version"], "1.5.1")


if __name__ == "__main__":
    unittest.main()
