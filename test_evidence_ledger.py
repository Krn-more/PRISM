import csv
import unittest
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import models
from evidence_ledger import EvidenceValidationError, attach_evidence_ledger, export_evidence_csv, export_evidence_dataframe, ingest_evidence_records, validate_evidence_record
from excel_formatter import format_excel_output


class EvidenceLedgerTests(unittest.TestCase):
    def setUp(self):
        self.engine = create_engine("sqlite:///:memory:")
        models.Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        with open(Path("fixtures") / "stage4_endpoint_evidence_valid.csv", newline="", encoding="utf-8") as handle:
            self.fixture = list(csv.DictReader(handle))

    def tearDown(self):
        self.session.close()
        self.engine.dispose()

    def test_valid_fixture_imports_every_endpoint_family_with_provenance(self):
        outcomes = ingest_evidence_records(self.session, self.fixture, "test.reviewer")
        self.assertEqual([outcome["status"] for outcome in outcomes], ["Imported"] * 5)
        exported = export_evidence_dataframe(self.session, ["TEST-INCHIKEY"])
        self.assertEqual(set(exported["endpoint"]), {"Genotoxicity", "Sensitisation", "Systemic toxicity", "Reproductive/developmental toxicity", "Carcinogenicity"})
        self.assertTrue(exported["raw_record_hash"].str.fullmatch(r"[0-9a-f]{64}").all())
        self.assertTrue((exported["vocabulary_version"] == "1.0.0").all())

    def test_duplicate_import_is_idempotent(self):
        ingest_evidence_records(self.session, [self.fixture[0]], "test.reviewer")
        outcome = ingest_evidence_records(self.session, [self.fixture[0]], "test.reviewer")[0]
        self.assertEqual(outcome["status"], "Duplicate")
        self.assertEqual(len(export_evidence_dataframe(self.session)), 1)

    def test_missing_source_url_is_retained_as_incomplete(self):
        record = dict(self.fixture[0])
        record["source_url"] = ""
        outcome = ingest_evidence_records(self.session, [record], "test.reviewer")[0]
        exported = export_evidence_dataframe(self.session)
        self.assertEqual(outcome["status"], "Imported")
        self.assertEqual(exported.iloc[0]["evidence_status"], "Incomplete")
        self.assertIn("provenance is incomplete", exported.iloc[0]["notes"])

    def test_structured_dose_representations_are_retained_and_invalid_numeric_rejected(self):
        for dose, unit in (("Not reported", ""), ("10-50", "mg/kg"), ("0.1 mg/L", ""), ("", "")):
            record = dict(self.fixture[0])
            record["dose"], record["dose_unit"] = dose, unit
            validate_evidence_record(record)
        invalid = dict(self.fixture[0])
        invalid["dose"], invalid["dose_unit"] = "-1", "mg/kg"
        with self.assertRaisesRegex(EvidenceValidationError, "non-negative"):
            validate_evidence_record(invalid)

    def test_no_record_is_not_negative_and_corrections_supersede(self):
        no_record = dict(self.fixture[3])
        self.assertEqual(validate_evidence_record(no_record)["evidence_status"], "No record retrieved")
        invalid = dict(no_record)
        invalid["result"] = "Negative"
        with self.assertRaisesRegex(EvidenceValidationError, "not a negative"):
            validate_evidence_record(invalid)
        first = ingest_evidence_records(self.session, [self.fixture[0]], "test.reviewer")[0]
        corrected = dict(self.fixture[0])
        corrected["evidence_id"] = "CORRECTED-001"
        corrected["notes"] = "Corrected record"
        corrected["supersedes_evidence_id"] = first["evidence_id"]
        ingest_evidence_records(self.session, [corrected], "test.reviewer")
        self.assertEqual(len(export_evidence_dataframe(self.session)), 2)

    def test_changed_source_record_requires_valid_supersession(self):
        original = ingest_evidence_records(self.session, [self.fixture[0]], "test.reviewer")[0]
        changed = dict(self.fixture[0])
        changed["notes"] = "Corrected interpretation"
        with self.assertRaisesRegex(EvidenceValidationError, "requires supersedes_evidence_id"):
            ingest_evidence_records(self.session, [changed], "test.reviewer")
        changed["evidence_id"] = "CORRECTED-002"
        changed["supersedes_evidence_id"] = original["evidence_id"]
        self.assertEqual(ingest_evidence_records(self.session, [changed], "test.reviewer")[0]["status"], "Imported")

    def test_failed_batch_is_atomic(self):
        invalid = dict(self.fixture[1])
        invalid["endpoint"] = "Unsupported endpoint"
        with self.assertRaises(EvidenceValidationError):
            ingest_evidence_records(self.session, [self.fixture[0], invalid], "test.reviewer")
        self.assertTrue(export_evidence_dataframe(self.session).empty)

    def test_ledger_import_does_not_change_clustering_fields(self):
        source = pd.DataFrame({"Cluster ID": ["Cluster_1"], "Decision": ["PASS"], "InChIKey": ["TEST-INCHIKEY"]})
        baseline = source.copy(deep=True)
        self.assertTrue(source[["Cluster ID", "Decision"]].equals(baseline[["Cluster ID", "Decision"]]))

    def test_excel_export_contains_endpoint_evidence_worksheet(self):
        ingest_evidence_records(self.session, [self.fixture[0]], "test.reviewer")
        source = pd.DataFrame({"Cluster ID": ["Cluster_1"], "Decision": ["PASS"]})
        source.attrs["endpoint_evidence_ledger"] = export_evidence_dataframe(self.session)
        output = BytesIO()
        format_excel_output(source, output)
        with ZipFile(BytesIO(output.getvalue())) as workbook:
            self.assertIn(b"Endpoint Evidence Ledger", workbook.read("xl/workbook.xml"))

    def test_empty_evidence_export_still_contains_visible_ledger_sheet(self):
        source = pd.DataFrame({"Cluster ID": ["Cluster_1"], "Decision": ["PASS"]})
        source.attrs["endpoint_evidence_ledger"] = export_evidence_dataframe(self.session)
        output = BytesIO()
        format_excel_output(source, output)
        with ZipFile(BytesIO(output.getvalue())) as workbook:
            self.assertIn(b"Endpoint Evidence Ledger", workbook.read("xl/workbook.xml"))
            shared_strings = workbook.read("xl/sharedStrings.xml")
            self.assertIn(b"No matching local evidence records", shared_strings)

    def test_audit_csv_export_contains_all_provenance_columns(self):
        ingest_evidence_records(self.session, [self.fixture[0]], "test.reviewer")
        exported = export_evidence_csv(self.session)
        self.assertIn("raw_record_hash", exported.splitlines()[0])
        self.assertIn("source_url", exported.splitlines()[0])
        self.assertIn("Public fixture", exported)

    def test_dtxsid_is_matched_into_report_evidence(self):
        record = dict(self.fixture[0])
        record["identity_key"] = "DTXSID123"
        record["compound_mapping_method"] = "DTXSID"
        ingest_evidence_records(self.session, [record], "test.reviewer")
        report = pd.DataFrame({"Compound Name": ["test"], "DTXSID": ["DTXSID123"]})
        attached = attach_evidence_ledger(report, session_factory=self.Session)
        self.assertEqual(len(attached.attrs["endpoint_evidence_ledger"]), 1)


if __name__ == "__main__":
    unittest.main()
