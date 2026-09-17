"""Emit the current local PRISM classification payload for workbook rendering."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from build_classification_assessment_data import build_rows


if __name__ == "__main__":
    payload = build_rows(Path("List of Inchikey and SMILES.xlsx"))
    rendered = json.dumps(payload, ensure_ascii=False, default=str)
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(rendered, encoding="utf-8")
    else:
        print(rendered)
