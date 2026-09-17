import unittest

import pandas as pd

from workbook_contract import (
    RETIRED_WORKBOOK_COLUMNS,
    WORKBOOK_COLUMN_ORDER,
    apply_workbook_contract,
    build_workbook_contract_frame,
)


class WorkbookContractTests(unittest.TestCase):
    def test_contract_has_unique_column_order(self):
        self.assertEqual(len(WORKBOOK_COLUMN_ORDER), len(set(WORKBOOK_COLUMN_ORDER)))

    def test_retired_columns_are_not_in_frozen_contract(self):
        for column in RETIRED_WORKBOOK_COLUMNS:
            self.assertNotIn(column, WORKBOOK_COLUMN_ORDER)

    def test_apply_workbook_contract_orders_known_columns_and_keeps_extras(self):
        df = pd.DataFrame([{
            "Extra Field": "keep-me",
            "Chemical Class": "Alcohol",
            "Compound Name": "Tetrahydrofuran",
            "TTC Limit": "1.5 µg/kg bw/day",
        }])

        result = apply_workbook_contract(df)

        self.assertEqual(result.columns[0], "Compound Name")
        self.assertEqual(result.columns[1], "CASRN")
        self.assertIn("Extra Field", result.columns)
        self.assertNotIn("TTC Limit", result.columns)

    def test_contract_frame_includes_retired_fields_for_audit(self):
        contract_frame = build_workbook_contract_frame()

        self.assertIn("Workbook Order", contract_frame.columns)
        self.assertIn("Section", contract_frame.columns)
        self.assertIn("Column Headers - Detailed Analysis sheet", contract_frame.columns)
        retired_rows = contract_frame[contract_frame["Status"] == "retired"]
        self.assertGreaterEqual(len(retired_rows), len(RETIRED_WORKBOOK_COLUMNS))


if __name__ == "__main__":
    unittest.main()
