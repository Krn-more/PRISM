import os
import unittest
from io import BytesIO
from unittest.mock import patch

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import models
from chemont_adapter import enrich_with_chemont, fetch_chemont_by_inchikey, load_chemont_rules, normalise_classification
from cluster_evidence_matrix import build_cluster_evidence_matrix
from excel_formatter import format_excel_output
from schemas import AssessmentResponse


class FakeResponse:
    status_code = 200
    def __init__(self, payload): self.payload = payload
    def json(self): return self.payload


class ChemontAdapterTests(unittest.TestCase):
    def setUp(self):
        engine = create_engine("sqlite:///:memory:")
        models.Base.metadata.create_all(engine)
        self.Session = sessionmaker(bind=engine)
        self.session = self.Session()
        self.engine = engine
        self.rules = load_chemont_rules()
        self.inchikey = "QVGXLLKOCUKJST-UHFFFAOYSA-N"
        self.payload = {
            "classification_version": "2.1", "kingdom": {"name": "Organic compounds"},
            "superclass": {"name": "Organic oxygen compounds"}, "class": {"name": "Organooxygen compounds"},
            "subclass": {"name": "Ethers"}, "direct_parent": {"name": "Dialkyl ethers"},
            "molecular_framework": "Aliphatic acyclic compounds", "substituents": [{"name": "Ether"}, {"name": "Alkyl"}],
        }

    def tearDown(self):
        self.session.close(); self.engine.dispose()

    def test_rules_limit_external_identifier_to_inchikey(self):
        self.assertEqual(self.rules["allowed_identifier_columns"], ["InChIKey"])

    def test_normalisation_retains_hierarchy_and_version(self):
        value = normalise_classification(self.payload, self.rules)
        self.assertEqual(value["ChemOnt_Direct_Parent"], "Dialkyl ethers")
        self.assertEqual(value["ChemOnt_Substituents"], "Alkyl; Ether")
        self.assertEqual(value["ChemOnt_Classification_Version"], "2.1")

    def test_cache_prevents_repeat_external_request(self):
        calls = []
        def get(*args, **kwargs): calls.append(args); return FakeResponse(self.payload)
        first, first_status = fetch_chemont_by_inchikey(self.inchikey, self.session, self.rules, get)
        second, second_status = fetch_chemont_by_inchikey(self.inchikey, self.session, self.rules, get)
        self.assertEqual(first, second)
        self.assertEqual((first_status, second_status), ("Retrieved", "Cached"))
        self.assertEqual(len(calls), 1)

    def test_disabled_adapter_never_requests_or_changes_decision(self):
        source = pd.DataFrame([{"InChIKey": self.inchikey, "Cluster ID": 1, "Decision": "PASS"}])
        with patch.dict(os.environ, {"CHEMONT_CLASSYFIRE_ENABLED": "false"}, clear=False):
            output = enrich_with_chemont(source, self.session, self.rules, lambda *args, **kwargs: self.fail("request made"))
        self.assertEqual(output.loc[0, "Decision"], "PASS")
        self.assertTrue(output.loc[0, "ChemOnt_Retrieval_Status"].startswith("Unavailable:"))

    def test_live_failure_is_visible_and_does_not_change_cluster(self):
        source = pd.DataFrame([{"InChIKey": self.inchikey, "Cluster ID": 1, "Decision": "PASS"}])
        failure = type("Failure", (), {"status_code": 503, "json": lambda self: {}})()
        with patch.dict(os.environ, {"CHEMONT_CLASSYFIRE_ENABLED": "true"}, clear=False):
            output = enrich_with_chemont(source, self.session, self.rules, lambda *args, **kwargs: failure)
        self.assertEqual(output.loc[0, "Cluster ID"], 1)
        self.assertEqual(output.loc[0, "Decision"], "PASS")
        self.assertEqual(output.loc[0, "ChemOnt_Retrieval_Status"], "Unavailable: ChemontError")

    def test_enabled_adapter_appends_hierarchy_without_changing_decision(self):
        source = pd.DataFrame([{"InChIKey": self.inchikey, "Cluster ID": 1, "Decision": "PASS"}])
        with patch.dict(os.environ, {"CHEMONT_CLASSYFIRE_ENABLED": "true"}, clear=False):
            output = enrich_with_chemont(source, self.session, self.rules, lambda *args, **kwargs: FakeResponse(self.payload))
        self.assertEqual(output.loc[0, "ChemOnt_Class"], "Organooxygen compounds")
        self.assertEqual(output.loc[0, "Decision"], "PASS")

    def test_hierarchy_is_in_matrix_and_separate_excel_audit_sheet(self):
        source = pd.DataFrame([{"InChIKey": self.inchikey, "Standardized SMILES": "CCOCC", "Cluster ID": 1, "Decision": "PASS", "Cluster_Status": "Stable"}])
        with patch.dict(os.environ, {"CHEMONT_CLASSYFIRE_ENABLED": "true"}, clear=False):
            output = enrich_with_chemont(source, self.session, self.rules, lambda *args, **kwargs: FakeResponse(self.payload))
        matrix = build_cluster_evidence_matrix(output)
        self.assertTrue((matrix["Evidence_Type"] == "ChemOnt hierarchy").any())
        output.attrs["chemont_classifications"] = output[["InChIKey", "Standardized SMILES", "ChemOnt_Kingdom", "ChemOnt_Superclass", "ChemOnt_Class", "ChemOnt_Subclass", "ChemOnt_Direct_Parent", "ChemOnt_Molecular_Framework", "ChemOnt_Substituents", "ChemOnt_Classification_Version", "ChemOnt_Source", "ChemOnt_Retrieval_Status"]]
        buffer = BytesIO()
        format_excel_output(output, buffer)
        with pd.ExcelFile(buffer) as workbook:
            self.assertIn("ChemOnt Classification", workbook.sheet_names)

    def test_response_schema_accepts_chemont_classification_bundle(self):
        response = AssessmentResponse(
            identity={},
            structure={},
            functional_groups=[],
            parent_moiety="",
            chemical_class="",
            analogs=[],
            confidence="Low",
            rationale="",
            chemont_classification={
                "InChIKey": self.inchikey,
                "Standardized SMILES": "CCOCC",
                "ChemOnt_Retrieval_Status": "Unavailable: CHEMONT_CLASSYFIRE_ENABLED is false",
            },
        )
        self.assertIn("ChemOnt_Retrieval_Status", response.chemont_classification)


if __name__ == "__main__":
    unittest.main()
