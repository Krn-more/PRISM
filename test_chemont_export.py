import unittest

import pandas as pd

from chemont_export import attach_chemont_classifications, build_chemont_bundle, chemont_display_frame, chemont_fields_present


class ChemOntExportTests(unittest.TestCase):
    def test_attach_chemont_classifications_populates_attrs_from_present_columns(self):
        df = pd.DataFrame([{
            "InChIKey": "ABCDE-FGHIJ",
            "Standardized SMILES": "CCO",
            "ChemOnt_Kingdom": "Organic compounds",
            "ChemOnt_Class": "Organooxygen compounds",
            "ChemOnt_Retrieval_Status": "Retrieved",
        }])

        result = attach_chemont_classifications(df)

        self.assertIn("chemont_classifications", result.attrs)
        bundle = result.attrs["chemont_classifications"]
        self.assertListEqual(list(bundle.columns), [
            "InChIKey",
            "Standardized SMILES",
            "ChemOnt_Kingdom",
            "ChemOnt_Class",
            "ChemOnt_Retrieval_Status",
        ])

    def test_attach_chemont_classifications_skips_when_no_columns_present(self):
        df = pd.DataFrame([{"Compound Name": "Acetone"}])

        result = attach_chemont_classifications(df)

        self.assertNotIn("chemont_classifications", result.attrs)
        self.assertEqual(chemont_fields_present(result), [])

    def test_build_chemont_bundle_normalizes_required_keys(self):
        bundle = build_chemont_bundle(
            "ABCDE-FGHIJ",
            "CCO",
            "ClassyFire ChemOnt API",
            "Retrieved",
            ChemOnt_Class="Organooxygen compounds",
            ChemOnt_Subclass="Ethers",
        )

        self.assertEqual(bundle["InChIKey"], "ABCDE-FGHIJ")
        self.assertEqual(bundle["Standardized SMILES"], "CCO")
        self.assertEqual(bundle["ChemOnt_Source"], "ClassyFire ChemOnt API")
        self.assertEqual(bundle["ChemOnt_Retrieval_Status"], "Retrieved")
        self.assertEqual(bundle["ChemOnt_Class"], "Organooxygen compounds")
        self.assertEqual(bundle["ChemOnt_Subclass"], "Ethers")

    def test_display_frame_renames_visible_headers(self):
        df = pd.DataFrame([{
            "ChemOnt_Kingdom": "Organic compounds",
            "ChemOnt_Class": "Organooxygen compounds",
            "ChemOnt_Retrieval_Status": "Retrieved",
        }])

        displayed = chemont_display_frame(df)

        self.assertIn("ChemOnt Kingdom", displayed.columns)
        self.assertIn("ChemOnt Class", displayed.columns)
        self.assertIn("ChemOnt Retrieval Status", displayed.columns)


if __name__ == "__main__":
    unittest.main()
