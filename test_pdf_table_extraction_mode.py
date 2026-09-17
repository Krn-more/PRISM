import unittest
from pathlib import Path


class PdfTableExtractionModeTests(unittest.TestCase):
    def test_pdf_table_extraction_does_not_force_windows_process_pool(self):
        script = (Path(__file__).with_name("static") / "app.js").read_text(encoding="utf-8")
        start = script.index("document.getElementById('pdfContinue1')")
        end = script.index("document.getElementById('pdfContinue2')", start)
        extraction_handler = script[start:end]
        self.assertIn("const useParallel = false;", extraction_handler)
        self.assertNotIn("const useParallel = true;", extraction_handler)


if __name__ == "__main__":
    unittest.main()
