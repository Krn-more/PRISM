"""Tests for the guarded OPSIN identity-resolution fallback."""

import asyncio
import os
import unittest
from unittest.mock import patch

from opsin_identity import opsin_identity_from_payload, opsin_name_is_eligible
from schemas import AssessmentRequest
import services
from services_async import resolve_identity_async


OPSIN_PAYLOAD = {
    "status": "SUCCESS",
    "smiles": "CC(C)(C)O",
    "inchi": "InChI=1/C4H10O/c1-4(2,3)5/h5H,1-3H3",
    "stdInChIKey": "DKGAVHZHDRPRBM-UHFFFAOYSA-N",
}


class _SyncResponse:
    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload or {}
        self.text = text

    def json(self):
        return self._payload


class _AsyncResponse:
    def __init__(self, status, payload=None, text=""):
        self.status = status
        self._payload = payload or {}
        self._text = text

    async def __aenter__(self):
        return self

    async def __aexit__(self, *_args):
        return False

    async def json(self):
        return self._payload

    async def text(self):
        return self._text


class _AsyncSession:
    def __init__(self):
        self.urls = []

    def get(self, url, **_kwargs):
        self.urls.append(url)
        if "ebi.ac.uk/opsin/ws" in url:
            return _AsyncResponse(200, OPSIN_PAYLOAD)
        return _AsyncResponse(404)


class OpsinIdentityTests(unittest.TestCase):
    def test_eligibility_rejects_generic_analytical_names(self):
        for name in (
            "Branched alkane (C=09) related compound 01",
            "Cyanox 425 (loss of tert-Butyl)",
            "Pentafluoro aminio polydimethylsiloxane (n=03) related compound 01",
            "Pentafluoro aminio polydimethylsiloxane",
            "Cyanox 2246 tert-Butylcyclopentadiene",
        ):
            with self.subTest(name=name):
                eligible, _reason = opsin_name_is_eligible(name)
                self.assertFalse(eligible)

    def test_payload_is_rdkit_validated_and_canonicalized(self):
        result = opsin_identity_from_payload("2-methylpropan-2-ol", OPSIN_PAYLOAD)
        self.assertEqual(result["smiles"], "CC(C)(C)O")
        self.assertEqual(result["inchikey"], "DKGAVHZHDRPRBM-UHFFFAOYSA-N")
        self.assertEqual(result["status"], "OPSIN parsed — identity review recommended")

    def test_sync_resolution_uses_opsin_only_after_other_sources_fail(self):
        calls = []

        def fake_get(url, **_kwargs):
            calls.append(url)
            if "ebi.ac.uk/opsin/ws" in url:
                return _SyncResponse(200, OPSIN_PAYLOAD)
            return _SyncResponse(404)

        with patch.object(services, "_http_get", side_effect=fake_get), patch.dict(os.environ, {"EPA_API_KEY": ""}):
            result = services.resolve_identity(AssessmentRequest(compound_name="2-methylpropan-2-ol"))

        self.assertEqual(result["source"], "OPSIN systematic-name parser (unverified identity)")
        self.assertEqual(result["smiles"], "CC(C)(C)O")
        self.assertTrue(any("pubchem.ncbi.nlm.nih.gov" in url for url in calls))
        self.assertTrue(any("ebi.ac.uk/opsin/ws" in url for url in calls))

    def test_sync_resolution_does_not_call_opsin_for_generic_name(self):
        calls = []

        def fake_get(url, **_kwargs):
            calls.append(url)
            return _SyncResponse(404)

        with patch.object(services, "_http_get", side_effect=fake_get), patch.dict(os.environ, {"EPA_API_KEY": ""}):
            result = services.resolve_identity(AssessmentRequest(compound_name="Branched alkane (C=09) related compound 01"))

        self.assertEqual(result["status"], "Unresolved")
        self.assertFalse(any("ebi.ac.uk/opsin/ws" in url for url in calls))

    def test_async_resolution_uses_the_same_opsin_guard(self):
        async def run():
            session = _AsyncSession()
            with patch.dict(os.environ, {"EPA_API_KEY": ""}):
                result = await resolve_identity_async(
                    AssessmentRequest(compound_name="2-methylpropan-2-ol"), session
                )
            return result, session.urls

        result, urls = asyncio.run(run())
        self.assertEqual(result["source"], "OPSIN systematic-name parser (unverified identity)")
        self.assertTrue(any("ebi.ac.uk/opsin/ws" in url for url in urls))


if __name__ == "__main__":
    unittest.main()
