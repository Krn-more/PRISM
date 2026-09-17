import unittest
from io import BytesIO
from zipfile import ZipFile

import pandas as pd

from excel_formatter import SUMMARY_COLUMNS, format_excel_output
from workbook_contract import apply_workbook_contract, build_workbook_contract_frame


class WorkbookAlignmentTests(unittest.TestCase):
    def test_export_uses_sme_review_sort_order(self):
        frame = pd.DataFrame([
            {
                'Compound Name': 'Zeta family two',
                'Compatibility Group': 'Family_002',
                'Cluster ID': 'Cluster_001',
                'Cluster Size': 1,
            },
            {
                'Compound Name': 'Beta cluster two',
                'Compatibility Group': 'Family_001',
                'Cluster ID': 'Cluster_002',
                'Cluster Size': 2,
            },
            {
                'Compound Name': 'Alpha cluster two',
                'Compatibility Group': 'Family_001',
                'Cluster ID': 'Cluster_002',
                'Cluster Size': 2,
            },
            {
                'Compound Name': 'Unresolved compound',
                'Compatibility Group': 'Not assessable',
                'Cluster ID': 'Not assessable',
                'Cluster Size': 'Not assessable',
            },
        ])
        out_buffer = BytesIO()
        format_excel_output(frame, out_buffer)

        expected = [
            'Alpha cluster two', 'Beta cluster two', 'Zeta family two',
            'Unresolved compound',
        ]
        workbook = BytesIO(out_buffer.getvalue())
        detailed = pd.read_excel(workbook, sheet_name='Detailed Analysis')
        self.assertEqual(detailed['Compound Name'].tolist(), expected)
        summary = pd.read_excel(BytesIO(out_buffer.getvalue()), sheet_name='Summary')
        self.assertEqual(summary['Compound Name'].tolist(), expected)

    def test_contract_starts_with_status_definitions(self):
        contract = build_workbook_contract_frame()
        definitions = contract.loc[contract["Section"] == "Status definitions"]
        self.assertEqual(definitions["Column Headers - Detailed Analysis sheet"].tolist(), ["Computed", "Preview only", "Retired", "Cluster identifier convention"])
        self.assertTrue(definitions["Notes"].str.len().gt(40).all())

    def test_review_reason_reconciles_surrogate_and_unresolved_rows(self):
        df = pd.DataFrame([
            {
                "Identity Status": "Surrogate SMILES (Typo Corrected)",
                "Standardized SMILES": "C1CCOC1",
                "Original SMILES": "C1CCOC1",
                "Manual Review Flag": 1,
                "Manual Review Reason": "No valid standardized or original SMILES could be parsed.",
                "Cluster ID": "Cluster_1",
            },
            {
                "Identity Status": "Tier 3: Class Fallback Required (No Structure)",
                "Standardized SMILES": "",
                "Original SMILES": "",
                "Manual Review Flag": 0,
                "Manual Review Reason": "",
                "Cluster ID": "Cluster_2",
            },
        ])
        result = apply_workbook_contract(df)
        self.assertEqual(result.loc[0, "Manual Review Reason"], "Surrogate structure used; identity uncertainty requires toxicologist review.")
        self.assertEqual(result.loc[0, "Manual Review Flag"], 1)
        self.assertEqual(result.loc[1, "Manual Review Reason"], "No valid standardized or original SMILES could be parsed.")
        self.assertEqual(result.loc[1, "Cluster ID"], "Not assessable")

    def test_export_includes_frozen_contract_and_review_sheets(self):
        df = pd.DataFrame([{
            "Compound Name": "tetrahydrofuran",
            "CASRN": "109-99-9",
            "Concentration": 47,
            "Source": "fixture",
            "Cleaned Compound Name": "tetrahydrofuran",
            "Conc_Float": 47,
            "Resolved Name": "tetrahydrofuran",
            "Original SMILES": "C1CCOC1",
            "Standardized SMILES": "C1CCOC1",
            "InChIKey": "WYURNTSHIVDZCO-UHFFFAOYSA-N",
            "InChI": "InChI=1S/C4H8O/c1-2-4-5-3-1/h1-4H2",
            "Chemical Class": "Aliphatic Compound (Miscellaneous)",
            "Corrected Chemical Class": "Cyclic ether",
            "Final Classification Record": "Taxonomy path: Oxygen-containing Compound → Ether → Cyclic Ether",
            "Classification Input Source": "Standardized SMILES",
            "Classification Status": "Calculated",
            "Classification Rule ID": "STC-ETH-003",
            "Classification Rule Version": "1.0.0",
            "Manual Review Flag": "No",
            "Manual Review Reason": "None",
            "Detected Feature Profile": "Ether; Cyclic ether",
            "Confidence": "High",
            "Tanimoto Analogs Found": 9,
            "Primary Functional Group": "Ether",
            "Secondary Functional Groups": "",
            "Taxonomy Path": "Oxygen-containing Compound → Ether → Cyclic Ether",
            "Taxonomy Path Steps": "Oxygen-containing Compound → Ether → Cyclic Ether",
            "Taxonomy Hierarchy": "Parent Class: Oxygen-containing Compound | Functional Group: Ether | Subclass: Ether | Structural Type: Cyclic Ether",
            "Topology Profile": "Topology Class: Monocyclic | Topology Modifiers: Heterocyclic | Ring System: Monocyclic | Ring Profile: 1 ring(s), 0 aromatic ring(s)",
            "MW": 72.06,
            "LogP": 0.8,
            "TPSA": 9.23,
            "Min_EState": 1,
            "Max_EState": 4.94,
            "3D_Asphericity": 0.153,
            "3D_PMI1": 68.167,
            "3D_PMI2": 70.589,
            "3D_PMI3": 119.994,
            "3D_RadiusOfGyration": 1.339,
            "Structural_Evidence_Status": "Calculated",
            "Structural_Evidence_Version": "1.0.0",
            "HBD": 0,
            "HBA": 1,
            "Rotatable_Bonds": 0,
            "Formal_Charge": 0,
            "Ring_Count": 1,
            "Aromatic_Ring_Count": 0,
            "Fraction_CSP3": 1.0,
            "Heavy_Atom_Count": 5,
            "Molecular_Refractivity": 20.053,
            "Ionisation_Indicator": "Neutral",
            "Toxicophore_Profile": "None",
            "All_Structural_Alerts": "None",
            "Alerts": "None",
            "Scaffold": "C1CCOC1",
            "Cluster ID": "Cluster_4",
            "Cluster Size": 19,
            "Cluster_Status": "Valid Group",
            "Domain": "INSIDE DOMAIN",
            "Decision": "PASS",
            "Domain_Status": "Inside",
            "Domain_Reasons": "Representative scaffold C1CCOC1; coverage 100%",
            "Domain_Rule_Version": "1.0.0-review-only",
            "Nearest_Neighbour_Tanimoto": 1.0,
            "Cluster_Min_Pairwise_Tanimoto": 0.16,
            "Cluster_Median_Pairwise_Tanimoto": 1.0,
            "Cluster_Scaffold_Coverage": 1.0,
            "Cluster_Representative_Scaffold": "C1CCOC1",
            "Property_Outlier_Flags": "None",
            "Ionisation_Consistency": "Consistent",
            "Toxicophore_Consistency": "Consistent",
            "Uncertainty_Summary": "High: Biological, Exposure",
            "Uncertainty_Rule_Version": "1.1.0-review-only",
            "EPA_CTX_Status": "Unavailable: EPA_CTX_ENABLED is false or EPA_CTX_API_KEY is not configured",
            "EPA_CTX_DTXSID": "",
            "EPA_CTX_Source_Version": "",
            "ChemOnt_Kingdom": "Organic compounds",
            "ChemOnt_Superclass": "Organic oxygen compounds",
            "ChemOnt_Class": "Organooxygen compounds",
            "ChemOnt_Subclass": "Ethers",
            "ChemOnt_Direct_Parent": "Dialkyl ethers",
            "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
            "ChemOnt_Substituents": "Ether; Alkyl",
            "ChemOnt_Classification_Version": "2.1",
            "ChemOnt_Source": "ClassyFire ChemOnt API",
            "ChemOnt_Retrieval_Status": "Retrieved",
            "ChEMBL_Max_Phase": "",
            "ChEMBL_Targets": "Not supplied",
            "XAI_Status": "Unavailable",
            "Read_Across_Justification": "AI read-across summary unavailable.",
        }])
        df.attrs["workbook_contract"] = pd.DataFrame([{
            "Workbook Order": 1,
            "Section": "Input capture and identity resolution",
            "Column": "Compound Name",
            "Status": "computed",
            "Notes": "Smoke test contract",
        }])
        df.attrs["uncertainty_register"] = pd.DataFrame([{"Scope": "Cluster", "Category": "Domain"}])
        df.attrs["cluster_evidence_matrix"] = pd.DataFrame([{"Cluster_ID": "Cluster_4", "Evidence_ID": "E-1", "Matrix_Version": "1.0.0"}])
        df.attrs["chemont_classifications"] = pd.DataFrame([{
            "InChIKey": "WYURNTSHIVDZCO-UHFFFAOYSA-N",
            "Standardized SMILES": "C1CCOC1",
            "ChemOnt_Kingdom": "Organic compounds",
            "ChemOnt_Superclass": "Organic oxygen compounds",
            "ChemOnt_Class": "Organooxygen compounds",
            "ChemOnt_Subclass": "Ethers",
            "ChemOnt_Direct_Parent": "Dialkyl ethers",
            "ChemOnt_Molecular_Framework": "Aliphatic acyclic compounds",
            "ChemOnt_Substituents": "Ether; Alkyl",
            "ChemOnt_Classification_Version": "2.1",
            "ChemOnt_Source": "ClassyFire ChemOnt API",
            "ChemOnt_Retrieval_Status": "Retrieved",
        }])

        buffer = BytesIO()
        format_excel_output(df, buffer, raw_tables=[])

        with ZipFile(BytesIO(buffer.getvalue())) as archive:
            workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
            sheet_xml = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")

        self.assertIn("Workbook Contract", workbook_xml)
        self.assertIn("Cluster Evidence Matrix", workbook_xml)
        self.assertIn("Summary", workbook_xml)
        self.assertIn("Detailed Analysis", workbook_xml)
        self.assertNotIn('name="Compounds"', workbook_xml)
        self.assertNotIn("ChemOnt Classification", workbook_xml)
        self.assertIn("Uncertainty Register", workbook_xml)
        self.assertNotIn("Cramer Class", sheet_xml)
        self.assertNotIn("TTC Limit", sheet_xml)

        summary = pd.read_excel(BytesIO(buffer.getvalue()), sheet_name="Summary")
        self.assertEqual(list(summary.columns), SUMMARY_COLUMNS)
        self.assertEqual(summary.loc[0, "Compound Name"], "tetrahydrofuran")
        self.assertTrue(str(summary.loc[0, "Cluster_Class_Reasons"]).strip())


if __name__ == "__main__":
    unittest.main()
