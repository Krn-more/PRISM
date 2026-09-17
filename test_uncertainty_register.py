import unittest
from io import BytesIO
from zipfile import ZipFile

import pandas as pd

from excel_formatter import format_excel_output
from uncertainty_register import add_uncertainty_register, load_uncertainty_rules


class UncertaintyRegisterTests(unittest.TestCase):
    def _coherent_frame(self):
        return pd.DataFrame({
            "Compound Name": ["member_a", "member_b"],
            "Cluster ID": ["Cluster_1", "Cluster_1"],
            "Identity Status": ["Resolved", "Exact Match"],
            "Standardized SMILES": ["CCO", "CCCO"],
            "Structural_Evidence_Status": ["Calculated", "Calculated"],
            "All_Structural_Alerts": ["None", "None"],
            "Toxicophore_Profile": ["None", "None"],
            "Property_Outlier_Flags": ["None", "None"],
            "Domain_Status": ["Inside", "Inside"],
            "Domain_Reasons": ["Coherent", "Coherent"],
        })

    def test_every_cluster_has_seven_categories_and_actions(self):
        result, register = add_uncertainty_register(self._coherent_frame())
        cluster_rows = register[register["Scope"] == "Cluster"]
        self.assertEqual(set(cluster_rows["Category"]), set(load_uncertainty_rules()["categories"]))
        self.assertTrue(result["Uncertainty_Summary"].str.contains("High: Biological, Exposure").all())
        action_required = register[register["Level"].isin(["Medium", "High"])]
        self.assertTrue(action_required["Required_Action"].str.strip().ne("").all())

    def test_coherent_member_is_low_for_local_evidence_categories(self):
        _, register = add_uncertainty_register(self._coherent_frame())
        rows = register[(register["Scope"] == "Compound") & register["Scope_ID"].str.contains("member_a", regex=False)].set_index("Category")
        for category in ("Identity", "Structure", "Property", "Domain"):
            self.assertEqual(rows.at[category, "Level"], "Low")

    def test_surrogate_and_outside_domain_are_escalated(self):
        frame = self._coherent_frame().iloc[[0]].copy()
        frame.loc[frame.index[0], "Identity Status"] = "Surrogate SMILES (Homologue)"
        frame.loc[frame.index[0], "Domain_Status"] = "Outside"
        _, register = add_uncertainty_register(frame)
        rows = register[register["Scope"] == "Compound"].set_index("Category")
        self.assertEqual(rows.at["Identity", "Level"], "Medium")
        self.assertEqual(rows.at["Domain", "Level"], "High")
        self.assertIn("Do not use direct read-across", rows.at["Domain", "Required_Action"])

    def test_missing_exposure_is_high_and_excel_contains_register_sheet(self):
        result, register = add_uncertainty_register(self._coherent_frame())
        rows = register[(register["Scope"] == "Compound") & register["Scope_ID"].str.contains("member_a", regex=False)].set_index("Category")
        self.assertEqual(rows.at["Exposure", "Level"], "High")
        self.assertIn("Request concentration", rows.at["Exposure", "Required_Action"])
        output = BytesIO()
        format_excel_output(result, output)
        with ZipFile(BytesIO(output.getvalue())) as workbook:
            self.assertIn(b"Uncertainty Register", workbook.read("xl/workbook.xml"))

    def test_missing_identity_and_domain_evidence_are_high(self):
        frame = self._coherent_frame().iloc[[0]].drop(columns=["Identity Status", "Domain_Status"])
        _, register = add_uncertainty_register(frame)
        rows = register[register["Scope"] == "Compound"].set_index("Category")
        self.assertEqual(rows.at["Identity", "Level"], "High")
        self.assertIn("Resolve and document identity", rows.at["Identity", "Required_Action"])
        self.assertEqual(rows.at["Domain", "Level"], "High")
        self.assertIn("Run or resolve Stage 2", rows.at["Domain", "Required_Action"])

    def test_cluster_aggregation_is_stable_under_row_reordering(self):
        frame = self._coherent_frame()
        frame.loc[frame.index[0], "Domain_Status"] = "Outside"
        frame.loc[frame.index[1], "Domain_Status"] = "Not assessable"
        _, first = add_uncertainty_register(frame)
        _, reversed_result = add_uncertainty_register(frame.iloc[::-1].reset_index(drop=True))
        first_domain = first[(first["Scope"] == "Cluster") & (first["Category"] == "Domain")].iloc[0]
        reversed_domain = reversed_result[(reversed_result["Scope"] == "Cluster") & (reversed_result["Category"] == "Domain")].iloc[0]
        self.assertEqual(first_domain["Reason"], reversed_domain["Reason"])
        self.assertEqual(first_domain["Required_Action"], reversed_domain["Required_Action"])

    def test_compound_scope_ids_are_unique_for_duplicate_names(self):
        frame = self._coherent_frame()
        frame["Compound Name"] = "duplicate"
        _, register = add_uncertainty_register(frame)
        scope_ids = register[register["Scope"] == "Compound"]["Scope_ID"]
        self.assertEqual(scope_ids.nunique(), len(scope_ids) // 7)


if __name__ == "__main__":
    unittest.main()
