# Technical architecture

## Overview

PRISM is a local single-page web application hosted by a Python FastAPI backend. A Windows portable build packages the backend and browser assets into `PRISM_Beta.exe`; the browser connects only to the local PRISM service at `http://localhost:8055`.

```text
Browser UI
  └─ static/index.html + static/app.js
       └─ FastAPI endpoints in main.py
            ├─ PDF extraction and row normalization
            ├─ Identity resolution and structure standardization
            ├─ Deterministic structural classification
            ├─ Compatibility-gated structural clustering
            ├─ Workbook export and audit generation
            └─ Approved AI cluster-reasoning integration
```

## Processing flow

```text
Analytical PDF
  → raw table extraction
  → identity resolution from compound name/CASRN
  → standardized SMILES and InChIKey
  → deterministic PRISM structural classification
  → descriptors, alerts, scaffold, and taxonomy evidence
  → compatibility gate
  → Morgan/ECFP4 fingerprint and Tanimoto clustering
  → cluster evidence and review rationale
  → Excel workbook export
```

## Deterministic chemistry and grouping

PRISM uses RDKit to parse structures, generate Morgan/ECFP4 fingerprints (radius 2, 1,024 bits), and calculate Tanimoto similarity. Only structures with valid standardized SMILES are eligible for structural clustering.

Before clustering, PRISM applies compatibility gates based on structural class, functional-group profile, topology, charge state, and repeated-unit restrictions. Within a compatible space it applies:

1. Tanimoto similarity threshold of 0.60.
2. Butina clustering with reordering.
3. Complete-linkage refinement to prevent bridge-like clusters.
4. Molecular-weight, LogP, TPSA, ionisation, toxicophore, scaffold, and outlier evidence checks.

Cluster assignment is deterministic for a fixed input and configuration. A Cluster ID is a display and review identifier only; it is not a risk score or toxicity conclusion.

## Data stores

| Store | Purpose | Deployment treatment |
|---|---|---|
| `chem_engine.db` | Local application state and development data | Runtime data; do not commit to Git. |
| `chemical_taxonomies_v1_1.sqlite` | Offline exact-match local taxonomy evidence | Required beside the portable EXE; distribute as a versioned release asset. |
| `config/*.json` | Structural classification, clustering, uncertainty, ontology, and AI policy rules | Version-controlled application configuration. |

## VM deployment and database requirements

PRISM does **not** require SQL Server, Oracle, PostgreSQL, MySQL, or any other separate database server for a standalone VM deployment.

The portable application operates as follows:

```text
Windows VM
  └─ PRISM_Beta.exe
       ├─ starts a local FastAPI/Uvicorn service on port 8055
       ├─ serves the browser interface at http://localhost:8055
       ├─ reads chemical_taxonomies_v1_1.sqlite beside the EXE
       └─ creates/uses local runtime data, temporary files, uploads, and exports
```

### Required files and storage

| Requirement | How it works |
|---|---|
| `PRISM_Beta.exe` | Starts the local application on the VM. |
| `chemical_taxonomies_v1_1.sqlite` | Must remain in the same folder as the EXE, unless the `PRISM_TAXONOMY_DB` environment variable points to another approved path. It is read locally for exact InChIKey taxonomy matches. |
| Writable working location | The Windows account running PRISM needs write permission for application runtime data, temporary files, uploaded PDFs, generated workbooks, and any local SQLite state. A user-writable folder such as `C:\PRISM-Beta` is suitable for Beta deployment. |
| Local port 8055 | Used only by the VM's local browser to reach PRISM. It does not require an inbound database connection. |

`chem_engine.db` is a local SQLite runtime database. It is created or used locally when required by the workflow; it is not a remote database dependency and should not be shared through Git.

### When a central database is optional

A separate central database is needed only if Philips chooses to add capabilities outside the current standalone Beta scope, for example:

- shared multi-user run history;
- central audit retention and reporting;
- cross-VM job queues or workload management;
- central role-based access control; or
- centrally managed configuration and reference-data updates.

Until those capabilities are introduced, each VM runs independently with its own local application state and its own copy of the versioned taxonomy SQLite sidecar.

## Network services used by the portable application

The core classification and clustering logic runs locally. The portable application may contact the following services only when the relevant workflow and deployment configuration permit it.

| Service | Purpose | Data sent | Default role |
|---|---|---|---|
| PubChem PUG REST | Identity and structure lookup | CASRN or compound name | Identity-resolution source |
| EPA CompTox | Identity lookup | CASRN or compound name | Identity-resolution source |
| NCI Cactus | SMILES lookup fallback | Compound name | Identity-resolution fallback |
| OPSIN | Structure generation for valid systematic names | Compound name | Restricted fallback; generic/polymer names are excluded |
| Azure OpenAI | Optional grounded cluster-reasoning draft | Approved cluster evidence payload | Disabled unless explicitly enabled and approved |

Network identity-resolution failures must not change PRISM’s deterministic structural class or cluster membership. A record without a defensible standardized structure is routed to manual review rather than structurally grouped.

## Security and deployment considerations

- Bind the portable application to the local host deployment context; do not expose it through an unapproved reverse proxy or public network interface.
- Treat uploaded analytical PDFs and generated workbooks as sensitive project data.
- Keep `.env`, Azure credentials, API keys, and client secrets outside source control.
- Allow-list only the required outbound identity-resolution destinations and the approved Azure OpenAI endpoint when AI reasoning is enabled.
- The current UI references Lucide icons from `https://unpkg.com`; host this asset internally or bundle it before a fully offline or restricted-network deployment.
- Code-sign the EXE and distribute through an approved internal channel.
- Retain workbook outputs and the associated source PDF according to the applicable data-retention policy.

### Recommended Beta access controls

For a Philips-only Beta release, apply the following controls in the approved build and deployment process:

| Control | Recommended implementation | Security boundary |
|---|---|---|
| Managed-device / Windows-domain check | At application launch, verify that the current Windows device is joined to an approved Philips Active Directory or Entra-managed domain. Maintain the allowed domain names in signed or centrally managed configuration. | A local launch gate. It deters accidental use on unmanaged devices but is not a substitute for enterprise access control. |
| Philips home-page or internal-service check | Before enabling the workflow, validate reachability of an approved Philips internal URL over TLS. Do not use a public web-page redirect as the sole control. | A network-presence signal only. It does not prove the identity of the user or device. |
| Beta expiry date | Store a signed release-expiry date in the build or centrally managed configuration. After expiry, show a clear message and block new assessments until an approved updated build is installed. | Limits use of an obsolete Beta build. Local system-clock changes can bypass an unsigned offline-only date check. |
| Code signing | Authenticode-sign the executable and validate the publisher through the normal enterprise software-distribution process. | Establishes publisher integrity and supports endpoint-security allow-listing. |
| Enterprise deployment | Distribute through an approved Philips software-management channel with device/user targeting, version control, and removal capability. | Stronger than controls embedded in a local executable. |

Use the Windows-domain, internal-service, and expiry checks as defence-in-depth controls. For enforceable access restriction, rely on enterprise identity, managed-device compliance, code signing, and controlled software distribution rather than a check implemented only inside the executable.

## Audit boundaries

The workbook records identity status, classification rule/version, structural evidence, cluster membership rationale, uncertainty, and applicable AI provenance. These provide traceability for review. They do not replace SME judgement, endpoint-specific evidence, or formal toxicological assessment.
