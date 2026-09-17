import asyncio
import csv
import io
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import main
import models
from evidence_ledger import export_evidence_dataframe


class FakeRequest:
    def __init__(self, payload):
        self.payload = payload

    async def json(self):
        return self.payload


class FakeUploadFile:
    def __init__(self, content: bytes):
        self.content = content

    async def read(self):
        return self.content


async def response_body(response) -> bytes:
    body = b""
    async for chunk in response.body_iterator:
        body += chunk.encode("utf-8") if isinstance(chunk, str) else chunk
    return body


class EvidenceApiWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.engine = create_engine("sqlite:///:memory:")
        models.Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)
        with open(Path("fixtures") / "stage4_endpoint_evidence_valid.csv", newline="", encoding="utf-8") as handle:
            self.records = list(csv.DictReader(handle))

    def tearDown(self):
        self.engine.dispose()

    def test_json_csv_duplicate_invalid_and_export_workflow(self):
        with patch("database.SessionLocal", self.Session):
            first = asyncio.run(main.import_endpoint_evidence(FakeRequest({"records": [self.records[0]], "imported_by": "workflow.test"})))
            self.assertEqual(first["outcomes"][0]["status"], "Imported")

            csv_content = io.StringIO()
            writer = csv.DictWriter(csv_content, fieldnames=self.records[1].keys())
            writer.writeheader()
            writer.writerows(self.records[1:])
            second = asyncio.run(main.import_endpoint_evidence_csv(FakeUploadFile(csv_content.getvalue().encode("utf-8")), "workflow.test"))
            self.assertEqual([outcome["status"] for outcome in second["outcomes"]], ["Imported"] * 4)

            duplicate = asyncio.run(main.import_endpoint_evidence(FakeRequest({"records": [self.records[0]], "imported_by": "workflow.test"})))
            self.assertEqual(duplicate["outcomes"][0]["status"], "Duplicate")

            invalid = dict(self.records[0])
            invalid["source_record"] = "INVALID-DOSE"
            invalid["dose"] = "-1"
            invalid["dose_unit"] = "mg/kg"
            with self.assertRaises(HTTPException) as raised:
                asyncio.run(main.import_endpoint_evidence(FakeRequest({"records": [invalid], "imported_by": "workflow.test"})))
            self.assertEqual(raised.exception.status_code, 400)

            exported_response = asyncio.run(main.export_endpoint_evidence_csv("TEST-INCHIKEY"))
            exported = asyncio.run(response_body(exported_response)).decode("utf-8")
            self.assertIn("raw_record_hash", exported)
            self.assertIn("supersedes_evidence_id", exported)
            self.assertIn("GEN-001", exported)
            session = self.Session()
            try:
                self.assertEqual(len(export_evidence_dataframe(session, ["TEST-INCHIKEY"])), 5)
            finally:
                session.close()


if __name__ == "__main__":
    unittest.main()
