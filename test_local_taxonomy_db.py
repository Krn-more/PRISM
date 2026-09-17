import os
import unittest
from pathlib import Path
from unittest.mock import patch

import local_taxonomy_db


class LocalTaxonomyDatabaseTests(unittest.TestCase):
    def test_frozen_app_prefers_sidecar_beside_executable(self):
        executable = Path(r"C:\PRISM-Beta\PRISM_Beta.exe")
        sidecar = executable.parent / local_taxonomy_db.SIDECAR_FILENAME
        with (
            patch.dict(os.environ, {"PRISM_TAXONOMY_DB": ""}),
            patch.object(local_taxonomy_db.sys, "frozen", True, create=True),
            patch.object(local_taxonomy_db.sys, "executable", str(executable)),
            patch.object(Path, "is_file", return_value=True),
        ):
            self.assertEqual(local_taxonomy_db.configured_database_path(), sidecar)

    def test_explicit_database_path_overrides_portable_sidecar(self):
        explicit = Path(r"C:\PRISM-Beta\reference.sqlite")
        with patch.dict(os.environ, {"PRISM_TAXONOMY_DB": str(explicit)}):
            self.assertEqual(local_taxonomy_db.configured_database_path(), explicit)


if __name__ == "__main__":
    unittest.main()
