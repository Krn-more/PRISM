import os
import unittest
from unittest.mock import patch

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import models
from epa_ctx_adapter import enrich_with_epa_ctx, fetch_hazard_by_dtxsid, load_epa_ctx_rules, map_hazard_payload_to_ledger
from evidence_ledger import attach_evidence_ledger, export_evidence_dataframe


class FakeResponse:
    status_code = 200

    def __init__(self, payload):
        self.payload = payload

    def json(self):
        return self.payload


class EpaCtxAdapterTests(unittest.TestCase):
    def setUp(self):
        self.engine = create_engine("sqlite:///:memory:")
        models.Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.rules = load_epa_ctx_rules()
        self.dtxsid = "DTXSID7020182"
        self.payload = [{
            "id": "EPA-TEST-001", "dtxsid": self.dtxsid,
            "effect": "Genotoxicity", "toxvalNumeric": 1.0,
            "toxvalUnits": "mg/kg", "studyType": "In vivo",
            "species": "Rat", "exposureRoute": "Oral",
        }]

    def tearDown(self):
        self.session.close()
        self.engine.dispose()

    def test_rules_are_complete_and_identifier_scope_is_public_only(self):
        self.assertEqual(self.rules["rule_version"], "1.0.0-read-only")
        self.assertEqual(self.rules["allowed_identifier_columns"], ["DTXSID", "InChIKey", "CAS", "CASRN"])
        self.assertGreater(self.rules["min_request_interval_seconds"], 0)

    def test_hazard_response_is_cached_without_a_second_request(self):
        calls = []
        def post(*args, **kwargs):
            calls.append((args, kwargs))
            return FakeResponse(self.payload)
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            first, first_status = fetch_hazard_by_dtxsid(self.dtxsid, self.session, self.rules, post)
            second, second_status = fetch_hazard_by_dtxsid(self.dtxsid, self.session, self.rules, post)
        self.assertEqual(first, self.payload)
        self.assertEqual(second, self.payload)
        self.assertEqual((first_status, second_status), ("Retrieved", "Cached"))
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][1]["headers"]["x-api-key"], "test-key")

    def test_mapping_retains_only_unambiguous_endpoint_evidence(self):
        mapped = map_hazard_payload_to_ledger(self.dtxsid, self.payload + [{"id": "EPA-UNKNOWN", "effect": "Body weight"}], self.rules)
        self.assertEqual(len(mapped), 1)
        self.assertEqual(mapped[0]["endpoint"], "Genotoxicity")
        self.assertEqual(mapped[0]["result"], "Not reported")
        self.assertEqual(mapped[0]["compound_mapping_method"], "DTXSID")

    def test_no_record_response_is_retained_as_status_not_a_negative_conclusion(self):
        source = pd.DataFrame([{"DTXSID": self.dtxsid, "Cluster": "Cluster_1", "Decision": "PASS"}])
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            output = enrich_with_epa_ctx(source, self.session, self.rules, lambda *args, **kwargs: FakeResponse([]))
        self.assertEqual(output.loc[0, "EPA_CTX_Status"], "Retrieved: 0 Stage 4-mappable record(s)")
        self.assertTrue(export_evidence_dataframe(self.session).empty)

    def test_service_failure_is_a_visible_status_and_not_a_pipeline_failure(self):
        source = pd.DataFrame([{"DTXSID": self.dtxsid, "Cluster": "Cluster_1", "Decision": "PASS"}])
        class FailureResponse:
            status_code = 503
            def json(self):
                return {"message": "temporarily unavailable"}
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            output = enrich_with_epa_ctx(source, self.session, self.rules, lambda *args, **kwargs: FailureResponse())
        self.assertEqual(output.loc[0, "Cluster"], "Cluster_1")
        self.assertEqual(output.loc[0, "Decision"], "PASS")
        self.assertEqual(output.loc[0, "EPA_CTX_Status"], "Unavailable: EpaCtxError")

    def test_auth_rate_limit_timeout_and_malformed_payloads_fail_safely(self):
        source = pd.DataFrame([{"DTXSID": self.dtxsid, "Cluster": "Cluster_1", "Decision": "PASS"}])
        class BadJsonResponse:
            status_code = 200
            def json(self):
                raise ValueError("not json")
        cases = [
            lambda *args, **kwargs: type("Auth", (), {"status_code": 401, "json": lambda self: {}})(),
            lambda *args, **kwargs: type("RateLimit", (), {"status_code": 429, "json": lambda self: {}})(),
            lambda *args, **kwargs: (_ for _ in ()).throw(TimeoutError("timeout")),
            lambda *args, **kwargs: BadJsonResponse(),
        ]
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            for post in cases:
                output = enrich_with_epa_ctx(source, self.session, self.rules, post)
                self.assertEqual(output.loc[0, "Cluster"], "Cluster_1")
                self.assertEqual(output.loc[0, "Decision"], "PASS")
                self.assertTrue(output.loc[0, "EPA_CTX_Status"].startswith("Unavailable:"))

    def test_disabled_adapter_never_requests_and_preserves_existing_decision(self):
        source = pd.DataFrame([{"DTXSID": self.dtxsid, "Cluster": "Cluster_1", "Decision": "PASS"}])
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "false", "EPA_CTX_API_KEY": ""}, clear=False):
            output = enrich_with_epa_ctx(source, self.session, self.rules, lambda *args, **kwargs: self.fail("request made"))
        self.assertEqual(output.loc[0, "Cluster"], "Cluster_1")
        self.assertEqual(output.loc[0, "Decision"], "PASS")
        self.assertTrue(output.loc[0, "EPA_CTX_Status"].startswith("Unavailable:"))

    def test_enabled_adapter_stores_provenance_without_changing_grouping(self):
        source = pd.DataFrame([{"DTXSID": self.dtxsid, "Cluster": "Cluster_1", "Decision": "PASS"}])
        def post(url, **kwargs):
            self.assertIn("hazard/toxval/search/by-dtxsid", url)
            return FakeResponse(self.payload)
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            output = enrich_with_epa_ctx(source, self.session, self.rules, post)
        evidence = export_evidence_dataframe(self.session, [self.dtxsid])
        self.assertEqual(output.loc[0, "Cluster"], "Cluster_1")
        self.assertEqual(output.loc[0, "Decision"], "PASS")
        self.assertEqual(output.loc[0, "EPA_CTX_DTXSID"], self.dtxsid)
        self.assertTrue(output.loc[0, "EPA_CTX_Status"].startswith("Retrieved:"))
        self.assertEqual(list(evidence["endpoint"]), ["Genotoxicity"])
        self.assertEqual(evidence.loc[0, "reviewer_status"], "Not reviewed")

    def test_resolved_dtxsid_is_used_when_attaching_the_report_ledger(self):
        source = pd.DataFrame([{"CAS": "100-00-0", "Cluster": "Cluster_1", "Decision": "PASS"}])
        def post(url, **kwargs):
            if "chemical/search/equal" in url:
                return FakeResponse([{"dtxsid": self.dtxsid}])
            return FakeResponse(self.payload)
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            output = enrich_with_epa_ctx(source, self.session, self.rules, post)
        attached = attach_evidence_ledger(output, session_factory=self.Session)
        ledger = attached.attrs["endpoint_evidence_ledger"]
        self.assertEqual(output.loc[0, "EPA_CTX_DTXSID"], self.dtxsid)
        self.assertEqual(list(ledger["identity_key"]), [self.dtxsid])

    def test_changed_epa_source_record_creates_a_superseding_ledger_row(self):
        source = pd.DataFrame([{"DTXSID": self.dtxsid, "Cluster": "Cluster_1", "Decision": "PASS"}])
        fresh_rules = dict(self.rules)
        fresh_rules["cache_ttl_days"] = -1
        payloads = [self.payload, [{**self.payload[0], "toxvalNumeric": 2.0}]]
        def post(url, **kwargs):
            return FakeResponse(payloads.pop(0))
        with patch.dict(os.environ, {"EPA_CTX_ENABLED": "true", "EPA_CTX_API_KEY": "test-key"}, clear=False):
            enrich_with_epa_ctx(source, self.session, fresh_rules, post)
            enrich_with_epa_ctx(source, self.session, fresh_rules, post)
        evidence = export_evidence_dataframe(self.session, [self.dtxsid])
        self.assertEqual(len(evidence), 2)
        self.assertEqual(evidence.iloc[1]["supersedes_evidence_id"], evidence.iloc[0]["evidence_id"])


if __name__ == "__main__":
    unittest.main()
