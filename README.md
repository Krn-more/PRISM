# PRISM — Extractables & Leachables Assessment

PRISM is a local Windows application for extracting chemical records from analytical PDF reports, resolving chemical identity, applying deterministic structure-based classification, grouping compatible chemicals, and exporting an auditable Excel workbook for SME review.

> **Beta status.** PRISM supports scientific review and traceability. It does not establish toxicological endpoint equivalence, safety, or a release decision. SME review remains required for unresolved identities, surrogate structures, borderline clusters, and read-across use.

## What PRISM does

1. Extracts tables, captions, compound names, CASRN, concentrations, and source references from analytical PDFs.
2. Resolves usable chemical structures from names and CASRN where evidence is available.
3. Standardizes structures and calculates deterministic classifications, functional groups, descriptors, structural alerts, and taxonomy evidence.
4. Groups compatible structures with Morgan/ECFP4 fingerprints, Tanimoto similarity, Butina clustering, and complete-linkage refinement.
5. Exports a workbook containing raw extraction evidence, detailed analysis, cluster evidence, uncertainty, and audit sheets.

## Technology stack

| Layer | Technology | Use |
|---|---|---|
| Backend API | Python, FastAPI, Uvicorn | Local web server and assessment endpoints |
| Front end | HTML, CSS, vanilla JavaScript, Lucide icons | Browser user interface served by the local API |
| Chemistry | RDKit | Structure parsing, normalization, descriptors, SMARTS rules, Morgan fingerprints, and Tanimoto similarity |
| Data processing | pandas, NumPy, scikit-learn | Tabular processing, clustering refinement, and analysis |
| PDF extraction | pdfplumber, PyMuPDF | Extraction of analytical-report tables and text |
| Excel generation | XlsxWriter, openpyxl | Final `.xlsx` workbook creation and verification |
| Local persistence | SQLite, SQLAlchemy | Application state and offline taxonomy lookup |
| Local taxonomy | `chemical_taxonomies_v1_1.sqlite` | Exact InChIKey reference lookup for ontology evidence |
| Optional AI | Azure OpenAI, Azure Identity | Grounded cluster-reasoning draft generation when explicitly configured and enabled |
| Packaging | PyInstaller | Windows portable executable build |

## Quick start: developer installation

### Prerequisites

- Windows 10 or 11, 64-bit
- Python 3.11, 64-bit, available on `PATH`
- Internet access only where public identity-resolution services or optional Azure AI are permitted

```powershell
git clone <approved-repository-url> PRISM
Set-Location PRISM

py -3.11 -m venv .build-venv
.\.build-venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

Copy-Item .env.example .env
python main.py
```

Open `http://localhost:8055` in a browser.

For a concise build and deployment guide, see [docs/BUILD_AND_INSTALL.md](docs/BUILD_AND_INSTALL.md). For implementation architecture and network/security considerations, see [docs/TECHNICAL_ARCHITECTURE.md](docs/TECHNICAL_ARCHITECTURE.md).

## Portable Beta installation

The portable package must contain these two files in the **same folder**:

```text
PRISM_Beta.exe
chemical_taxonomies_v1_1.sqlite
```

Unzip both files to a user-writable folder, for example `C:\PRISM-Beta`, then launch `PRISM_Beta.exe`. The application starts a local server and opens at `http://localhost:8055`.

Do not rename or separate the SQLite sidecar. Without it, PRISM still runs but cannot provide its local exact-match taxonomy evidence.

## Repository structure

```text
static/                         Browser UI assets
main.py                         FastAPI application and workflow orchestration
extractor.py                    PDF/table extraction
services.py                     Identity resolution using approved public sources
structure_classification.py     Deterministic structure classification rules
clustering_engine.py            Compatibility gates and structural grouping
excel_formatter.py              Final workbook schema and formatting
workbook_contract.py            Workbook field definitions and output gate
local_taxonomy_db.py            Offline SQLite taxonomy lookup
build_local_taxonomy_db.py      SQLite sidecar builder
config/                         Versioned local rule and policy files
fixtures/                       Regression and workflow fixtures
test_*.py                       Automated tests
prism_beta.spec                 PyInstaller specification for PRISM_Beta.exe
```

## Chemical grouping summary

PRISM groups only records with a valid standardized structure. It first applies structural compatibility gates, then clusters unique structures using Morgan/ECFP4 fingerprints (radius 2, 1,024 bits), a Tanimoto threshold of 0.60, Butina clustering, and complete-linkage refinement. It maps assignments back to source rows after deduplication.

Cluster IDs are membership labels for review and sorting; they are not similarity scores, priority rankings, or toxicological measurements.

## Configuration and secrets

The committed [.env.example](.env.example) is an empty, credential-free template. It contains only blank Azure AI variable names and disabled AI flags. Never commit a real `.env` file.

The existing Portable Beta EXE does **not** require `.env` to run the local PDF extraction, identity resolution, deterministic structural classification, clustering, and Excel export workflow. `.env` is not embedded in `PRISM_Beta.exe`.

Create `.env` beside `PRISM_Beta.exe` only when Azure AI cluster reasoning is approved and configured:

```powershell
Copy-Item .env.example .env
```

Populate the approved Azure values and set both `XAI_ENABLED=true` and `XAI_STAGE6_ENABLED=true`. Without these settings, PRISM remains fully usable in deterministic local mode and records AI reasoning as unavailable rather than blocking the workbook workflow.

The default mode is deterministic local processing. Public identity-resolution lookups are used only when permitted by the deployment. See [docs/TECHNICAL_ARCHITECTURE.md](docs/TECHNICAL_ARCHITECTURE.md) for the current data flow and network allow-list information.

## Testing

Run the automated suite from the activated virtual environment:

```powershell
python -m unittest discover -p "test_*.py"
```

For the main end-to-end regression fixture:

```powershell
python -m unittest test_end_to_end_smoke.EndToEndSmokeTest.test_frozen_fixture_runs_through_batch_pipeline_and_workbook_export
```

## GitHub transfer rules

- Commit source code, `static/`, `config/`, fixtures, tests, `.env.example`, and documentation.
- Do **not** commit `.env`, virtual environments, PyInstaller `build/` and `dist/` folders, runtime databases, temporary PDFs, logs, or generated workbooks.
- Do **not** commit the 445 MB source taxonomy text file or 514 MB SQLite sidecar to ordinary Git. Store them as approved release assets, an internal package artifact, or Git LFS if IT explicitly approves it.
- Build and sign the portable ZIP in the approved CI or controlled Windows build environment.

## License and intended use

Add the organisation-approved license, data-governance statement, and support contact before publishing outside the internal organisation.
