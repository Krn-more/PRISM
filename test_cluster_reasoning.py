import json
import unittest
from unittest.mock import patch

import pandas as pd

from cluster_reasoning import add_cluster_membership_evidence
import xai_engine


class FakeCompletion:
    def __init__(self, content):
        self.choices = [type("Choice", (), {"message": type("Message", (), {"content": content})()})()]


class FakeClient:
    def __init__(self, content):
        self.content = content
        self.chat = type("Chat", (), {"completions": self})()

    def create(self, **kwargs):
        return FakeCompletion(self.content)


class ClusterReasoningTests(unittest.TestCase):
    def setUp(self):
        frame = pd.DataFrame([
            {
                "Compound Name": "ethanol", "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "Standardized SMILES": "CCO", "Cluster ID": "Cluster_1", "Cluster Size": 2,
                "Scaffold": "", "Cluster_Representative_Scaffold": "Not supplied",
                "Cluster_Scaffold_Coverage": 0.0, "Domain_Status": "Inside",
                "Ionisation_Consistency": "Consistent", "Toxicophore_Consistency": "Consistent",
                "Property_Outlier_Flags": "None", "Identity Status": "Resolved",
                "Corrected Chemical Class": "Alcohol", "Primary Functional Group": "Alcohol",
                "Taxonomy Path": "Oxygen-containing Compound → Alcohol", "Domain_Reasons": "Within domain",
                "Uncertainty_Summary": "Low uncertainty across Stage 3 categories",
            },
            {
                "Compound Name": "propanol", "InChIKey": "BDERNNFJNOPAEC-UHFFFAOYSA-N",
                "Standardized SMILES": "CCCO", "Cluster ID": "Cluster_1", "Cluster Size": 2,
                "Scaffold": "", "Cluster_Representative_Scaffold": "Not supplied",
                "Cluster_Scaffold_Coverage": 0.0, "Domain_Status": "Inside",
                "Ionisation_Consistency": "Consistent", "Toxicophore_Consistency": "Consistent",
                "Property_Outlier_Flags": "None", "Identity Status": "Resolved",
                "Corrected Chemical Class": "Alcohol", "Primary Functional Group": "Alcohol",
                "Taxonomy Path": "Oxygen-containing Compound → Alcohol", "Domain_Reasons": "Within domain",
                "Uncertainty_Summary": "Low uncertainty across Stage 3 categories",
            },
        ])
        self.result = add_cluster_membership_evidence(frame)

    def test_deterministic_membership_evidence_is_row_specific_and_requires_sme_review(self):
        row = self.result.iloc[0]
        self.assertEqual(row["Cluster_Reasoning_Status"], "Calculated deterministic evidence")
        self.assertIn("propanol", row["Cluster_Nearest_Members"])
        self.assertIn("Compatibility-gated Morgan/ECFP4", row["Cluster_Membership_Rationale"])
        self.assertIn("SME review is required", row["SME_Review_Boundary"])
        self.assertIn("nearest_member_tanimoto", row["Cluster_Membership_Parameters"])

    def test_grounded_ai_explanation_requires_selected_evidence_fields_and_sme_boundary(self):
        row = self.result.iloc[0].to_dict()
        content = json.dumps({
            "cluster_membership_rationale": "The selected compound has a documented in-cluster structural neighbour and was retained by deterministic consensus.",
            "supporting_parameters": "Nearest-member Tanimoto and the recorded domain-consistency fields support the technical grouping rationale.",
            "assessment_boundary": "SME review is required before this grouping is used for analogue selection, read-across, or any endpoint-specific assessment.",
            "evidence_fields": ["Cluster_Membership_Rationale", "Cluster_Nearest_Members", "SME_Review_Boundary"],
        })
        with patch.object(xai_engine, "_build_client", return_value=(FakeClient(content), "test-deployment", "Test", "")):
            output = xai_engine.generate_compound_cluster_reasoning(row)
        self.assertEqual(output["XAI_Status"], "Generated grounded (Test)")
        self.assertIn("SME review is required", output["AI_Assessment_Boundary"])
        self.assertIn("Cluster_Nearest_Members", output["AI_Evidence_Fields"])

    def test_unassessable_structure_never_calls_ai(self):
        row = self.result.iloc[0].to_dict()
        row["Cluster_Reasoning_Status"] = "Not assessable"
        with patch.object(xai_engine, "_build_client") as build_client:
            output = xai_engine.generate_compound_cluster_reasoning(row)
        build_client.assert_not_called()
        self.assertEqual(output["XAI_Status"], "Not applicable")

    def test_batch_reasoning_accepts_fenced_json_and_records_validation_audit(self):
        content = """```json
        {"cluster_interpretation":"The deterministic structural-consensus record supports this technical grouping.","supporting_parameters":"The supplied compatibility and consistency fields support the grouping.","assessment_boundary":"SME review is required before any analogue or endpoint-specific use.","evidence_fields":["Cluster ID","Cluster_Membership_Rationale"],"gateway_trace":"ignored"}
        ```"""
        with patch.object(xai_engine, "_cache_get", return_value=None), \
             patch.object(xai_engine, "_cache_put"), \
             patch.object(xai_engine, "_build_client", return_value=(FakeClient(content), "test-deployment", "Test", "")):
            output = xai_engine.generate_batch_cluster_reasoning(self.result, enabled=True)
        self.assertTrue(output["XAI_Status"].str.startswith("Generated grounded (Test)").all())
        self.assertEqual(set(output["AI_Validation_Status"]), {"validated_json_object"})
        self.assertTrue(output["AI_Response_Hash"].str.match(r"^[0-9a-f]{64}$").all())


if __name__ == "__main__":
    unittest.main()
