# Build and installation guide

This guide is for IT teams building, validating, or distributing PRISM Beta on Windows.

## 1. Supported platform

| Item | Requirement |
|---|---|
| Operating system | Windows 10 or Windows 11, 64-bit |
| Build Python | Python 3.11, 64-bit |
| Browser | Modern Chromium-based browser, Edge, or Chrome |
| Local port | TCP port 8055 on loopback or local host |
| Disk space | At least 4 GB free during build; approximately 700 MB for the portable release and taxonomy sidecar |

PRISM is a local web application. The executable starts a FastAPI/Uvicorn server on port 8055 and serves the UI locally. It is not an internet-hosted web site.

### VM and database clarification

A standalone Windows VM does **not** need SQL Server, Oracle, PostgreSQL, MySQL, or another database server. The portable release needs only:

```text
PRISM_Beta.exe
chemical_taxonomies_v1_1.sqlite
```

Keep both files in the same writable folder. PRISM reads the taxonomy SQLite file locally and manages any runtime SQLite state, temporary files, uploaded PDFs, and exported workbooks on the VM. A central database is optional and is needed only for future multi-user history, central audit retention, or enterprise-scale workflow management.

## 2. Source installation

Open PowerShell in the repository root.

```powershell
py -3.11 -m venv .build-venv
.\.build-venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
Copy-Item .env.example .env
```

Start the application:

```powershell
python main.py
```

Browse to `http://localhost:8055`.

If PowerShell prevents activation, use Command Prompt instead:

```cmd
.build-venv\Scripts\activate.bat
```

## 3. Configuration

`.env` is local deployment configuration and must not be committed to Git.

| Setting | Required? | Purpose |
|---|---:|---|
| `XAI_ENABLED` and `XAI_STAGE6_ENABLED` | No | Enables Azure OpenAI cluster-reasoning drafts when explicitly approved. |
| `AZURE_OPENAI_ENDPOINT`, `AZURE_OPENAI_DEPLOYMENT_NAME` | Only for AI | Azure OpenAI deployment configuration. |
| `AZURE_OPENAI_API_KEY` or Azure service-principal values | Only for AI | Azure OpenAI authentication. |
| `PRISM_TAXONOMY_DB` | No | Overrides the location of the local taxonomy SQLite file. |

`DETERMINISTIC_LOCAL_MODE` defaults to `true` in the application. Core identity handling, structural classification, grouping, and workbook export remain deterministic and local.

## 4. Build the local taxonomy SQLite sidecar

The portable application expects `chemical_taxonomies_v1_1.sqlite` beside the executable. Build it from the approved ChemicalTaxonomies v1.1 tab-delimited source:

```powershell
python build_local_taxonomy_db.py `
  --source .\ChemicalTaxonomies\ChemicalTaxonomies_v1.1.txt `
  --output .\dist\chemical_taxonomies_v1_1.sqlite
```

This indexed file supports exact local InChIKey matches. It is evidence enrichment only; PRISM’s local structural rules remain authoritative for classification and review decisions.

## 5. Build the Windows executable

From the activated virtual environment and repository root:

```powershell
python -m PyInstaller --clean --noconfirm prism_beta.spec
```

Expected artifact:

```text
dist\PRISM_Beta.exe
```

The PyInstaller specification packages:

- `main.py` and required backend modules
- `static/` UI assets
- `config/` rule and policy files
- RDKit runtime dependencies
- FastAPI/Uvicorn and workbook-generation dependencies

The taxonomy SQLite database is deliberately not embedded in the executable because it is large and is independently versioned. It must be created or copied into `dist\` after the executable build.

## 6. Create the portable release ZIP

Verify that the two files are present:

```powershell
Get-Item .\dist\PRISM_Beta.exe
Get-Item .\dist\chemical_taxonomies_v1_1.sqlite
```

Create the distributable package:

```powershell
Compress-Archive `
  -Path .\dist\PRISM_Beta.exe, .\dist\chemical_taxonomies_v1_1.sqlite `
  -DestinationPath .\dist\PRISM_Beta_Portable.zip `
  -Force
```

The ZIP must contain exactly these two top-level files:

```text
PRISM_Beta.exe
chemical_taxonomies_v1_1.sqlite
```

## 7. Install the portable release

1. Unzip the package to a user-writable folder, for example `C:\PRISM-Beta`.
2. Keep the EXE and SQLite file together in that folder.
3. Run `PRISM_Beta.exe`.
4. If a browser does not open automatically, browse to `http://localhost:8055`.
5. Upload an analytical PDF and follow the guided workflow.

Do not install directly from a ZIP file, a read-only directory, or an untrusted email attachment. Do not move the SQLite sidecar elsewhere unless `PRISM_TAXONOMY_DB` is configured.

## 8. Build verification

Before distribution:

```powershell
python -m unittest discover -p "test_*.py"
python -m unittest test_end_to_end_smoke.EndToEndSmokeTest.test_frozen_fixture_runs_through_batch_pipeline_and_workbook_export
```

Then perform a clean-machine smoke test:

1. Extract the portable ZIP to a new folder.
2. Launch `PRISM_Beta.exe`.
3. Confirm `http://localhost:8055` loads.
4. Run one representative PDF through export.
5. Confirm the final workbook contains `Raw Extraction (Pre-Dedup)`, `Summary`, `Cluster Summary`, `Detailed Analysis`, and audit sheets.
6. Confirm the taxonomy SQLite file was found by checking populated local-taxonomy evidence for a known exact-match structure.

## 9. Troubleshooting

| Symptom | Likely cause | Action |
|---|---|---|
| Browser cannot open `localhost:8055` | Application did not start or port is occupied | Review the console output; release port 8055 or restart the application. |
| Taxonomy evidence is unavailable | SQLite sidecar is missing, renamed, or separated from the EXE | Place `chemical_taxonomies_v1_1.sqlite` beside `PRISM_Beta.exe`, or configure `PRISM_TAXONOMY_DB`. |
| AI reasoning is unavailable | AI is disabled or Azure configuration is absent | This is expected in deterministic local mode. Configure Azure values only after approval. |
| Some compounds are `Not assessable` | No defensible standardized structure was resolved | Review identity evidence; obtain an exact CASRN, name, or validated structure. |
| Security software blocks the EXE | Unsigned executable or organisational policy | Code-sign the release and distribute through the approved internal software channel. |

## 10. Release controls

- Build from a tagged commit in a controlled Windows environment.
- Record Python version, dependency lock/version set, Git commit SHA, and taxonomy sidecar version.
- Scan and code-sign the EXE according to organisational policy.
- Publish source and release artifacts separately.
- Store secrets only in the approved secret-management system.
