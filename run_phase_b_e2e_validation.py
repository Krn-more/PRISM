"""Run the live PRISM B1a -> B1b -> B2 sequence for a Phase A workbook.

Uses localhost deliberately so the test exercises the same persisted run state
and B2 validation gate as the browser workflow.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pandas as pd


SOURCE = Path("phase_a_D61148-1 C33511 Draft 1 (2).xlsx")
OUTPUT_DIR = Path("outputs/e2e_phase_b_validation_20260907")
URL = "http://localhost:8055/assess/pdf/run_deterministic_enrichment"


def _records(frame: pd.DataFrame) -> list[dict[str, object]]:
    safe = frame.where(pd.notna(frame), "")
    return safe.to_dict(orient="records")


def _post(payload: dict[str, object]) -> dict[str, object]:
    request = Request(
        URL,
        data=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urlopen(request, timeout=300) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Phase B HTTP {exc.code}: {detail}") from exc


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # Phase A exports intentionally split the primary records by confidence.
    # Recombine only those row sheets; do not accidentally feed evidence or
    # legacy external-enrichment tabs back into Phase B.
    phase_a_sheets = ["High Confidence", "Medium-Low Analogs", "Unresolved (Review)"]
    input_frame = pd.concat(
        [pd.read_excel(SOURCE, sheet_name=sheet_name) for sheet_name in phase_a_sheets],
        ignore_index=True,
    )
    payload: dict[str, object] = {
        "filename": SOURCE.name,
        "compounds": _records(input_frame),
        "phase_b_stage": "b1a",
    }
    b1a = _post(payload)
    b1b = _post({
        "filename": SOURCE.name,
        "compounds": b1a["compounds"],
        "run_id": b1a["run_id"],
        "phase_b_stage": "b1b4",
    })
    b2 = _post({
        "filename": SOURCE.name,
        "compounds": b1b["compounds"],
        "run_id": b1b["run_id"],
        "phase_b_stage": "b2",
    })
    if b2.get("phase_b_stage") != "B2" or not b2.get("excel_base64"):
        raise RuntimeError(f"B2 did not return an export: {b2}")

    target = OUTPUT_DIR / str(b2["filename"])
    target.write_bytes(base64.b64decode(str(b2["excel_base64"])))
    trace = {
        "input_rows": len(input_frame),
        "b1a_stage": b1a.get("phase_b_stage"),
        "b1b_stage": b1b.get("phase_b_stage"),
        "b2_stage": b2.get("phase_b_stage"),
        "run_id": b2.get("run_id"),
        "output": str(target),
        "b1a_trace": b1a.get("phase_b_trace", []),
        "b1b_trace": b1b.get("phase_b_trace", []),
        "b2_trace": b2.get("phase_b_trace", []),
    }
    (OUTPUT_DIR / "phase_b_e2e_trace.json").write_text(json.dumps(trace, indent=2), encoding="utf-8")
    print(json.dumps({key: trace[key] for key in ("input_rows", "b1a_stage", "b1b_stage", "b2_stage", "run_id", "output")}, indent=2))


if __name__ == "__main__":
    main()
