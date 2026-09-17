# Phase B Workbook Gap Remediation Plan

## Objective

Make the Excel workflow a complete, auditable continuation of Phase A and the
B1a/B2 local workflow, without changing the already-working identity,
classification, descriptor, or PDF extraction behavior.

The remediation is gated: the next stage cannot be signed off until the
acceptance criteria and regression tests for the current stage pass.

## Non-negotiable invariants

- Phase A identity and structure values remain unchanged.
- Valid standardized SMILES remain eligible for deterministic descriptors,
  alerts, topology, clustering, and domain assessment.
- Rows without a valid standardized structure never receive a similarity
  cluster, an inside-domain decision, or `PASS`.
- External services are supporting evidence only; an unavailable service must
  produce an explicit status and must not erase deterministic results.
- Optional B1b stages, when enabled, are independently rerunnable and preserve
  completed fields; they are not required for the local B1a → B2 path.
- B2 exports only the latest validated state and never silently exports a
  placeholder evidence sheet after a completed stage.
- Cramer/TTC fields remain retired unless separately approved later.

### Frozen status and null vocabulary

- `Cluster ID`: `Not assessable` when no valid structure exists.
- `Cluster Size`: `Not assessable` when no valid structure exists.
- `Cluster_Status`: `Review Required` for unresolved rows.
- `Domain`: `Not assessable` when no valid structure exists.
- `Domain_Status`: `Not assessable` when no valid structure exists.
- `Decision`: `REVIEW` when a structural prerequisite is missing.
- Valid structure with no finding: literal `None`.
- Missing/invalid structure: `Not assessable`.
- External failure: `Unavailable: <reason>`.
- Numeric zero remains numeric `0` or `0.0`, never blank.

## Consolidated implementation status (2026-09-03)

| Plan area | Status | Current assessment |
| --- | --- | --- |
| Stage 0 — Fixture/contract | Complete | Frozen 219-row baseline recorded: 190 structured, 29 unresolved; contract and sheet counts are reproducible. |
| Stage 1 — Unresolved-row eligibility | Complete | Invalid structures are excluded from clustering/domain assessment and cannot receive `PASS`. |
| Stage 2 — Missing/no-finding normalization | Complete | `None`, `Not assessable`, external `Unavailable`, and numeric-zero handling are implemented and tested. |
| Stage 3 — Persisted Phase B state | Complete | `run_id`, server-side evidence restoration, browser refresh recovery, and four-hour expiry cleanup are implemented. |
| Stage 3.5 — Review-reason reconciliation | Complete | Contradictory surrogate/review reasons are reconciled by the workbook contract. |
| Deterministic local mode | Implemented | ChemOnt/EPA external calls can be skipped without blocking local processing; statuses remain explicit. |
| Stage 4 — Bounded external enrichment | Deferred/out of scope | EPA CTX and ChemOnt are intentionally removed from the active workflow; deterministic local enrichment proceeds without them. No EPA 2a/2b or ChemOnt service work is pending for this release. |
| Stage 5 — Evidence-sheet population | Complete for in-scope layers | Local evidence sheets and explicit not-requested statuses are generated. EPA CTX/ChemOnt fields remain traceable as deferred, not as required blocking inputs. |
| Stage 6 — B2 validation gate | Complete | Required columns, prerequisites, unresolved safety, and evidence-bundle checks gate export. |
| Stage 7 — UI workflow | Complete | Ordered checkpoints, progress, retry-only failed stage, disabled completed buttons, and refresh recovery are implemented. |
| Stage 8 — Release/EXE validation | Substantially complete | Onedir EXE rebuilt and smoke-tested; targeted release tests pass. One stale end-to-end assertion still needs contract update; EPA CTX/ChemOnt are intentionally excluded. |
| Stage 9 — Final sign-off | Partial | Baseline and targeted checks are documented; bounded full post-Phase-B regeneration/render and final sign-off remain. |

## Stage 0 — Freeze the regression fixture and contract

### Work

1. Use `processed_D61148-1 C33511 Draft 1 (2).xlsx` as the initial regression
   fixture, retaining its original copy unchanged.
2. Record the approved population from the current fixture: 219 total rows,
   190 rows with valid standardized structures, and 29 unresolved rows.
3. Freeze the workbook contract and controlled status vocabulary:
   `Calculated`, `Not assessable`, `Review Required`, `PASS`, `WARNING`,
   `FAIL`, `Unavailable`, and `None`.
4. Add a machine-readable fixture summary for use in automated tests.
5. Record the expected sheet relationships: `Compounds` is the single compound
   population; the cluster matrix is one or more evidence rows per eligible
   cluster; external ChemOnt/EPA sheets are out of scope for this release.

### Acceptance criteria

- The fixture can be loaded without modifying it.
- Row counts and key sheet names are reproducible.
- The frozen counts and SHA-256 fingerprint are recorded in
  `phase_b_regression_fixture_summary.json`.
- Any older 43/145/31 count is treated as historical, not as the active
  acceptance baseline.
- Existing Phase A fields are byte/value-equivalent before and after a no-op
  export test.
- The baseline fixture is retained unchanged and all later comparisons are
  against a fresh generated copy.

## Stage 1 — Correct unresolved-row eligibility and domain decisions

### Implementation status (2026-09-03)

Validated against the frozen 219-row fixture: all 29 unresolved rows are
review-only, excluded from similarity assessment, and carry `Not assessable`
cluster/domain values. No unresolved row has `PASS` or `INSIDE DOMAIN`.

### Work

1. Define one shared `has_valid_structure` predicate based on a parseable
   standardized SMILES.
2. Exclude invalid/unresolved rows before distance-matrix construction and
   consensus clustering.
3. For excluded rows, set:
   - `Cluster ID`: `Not assessable`
   - `Cluster Size`: `Not assessable`
   - `Cluster_Status`: `Review Required`
   - `Domain`: `Not assessable`
   - `Decision`: `REVIEW`
   - `Domain_Status`: `Not assessable`
   - `Domain_Reasons`: `Standardized structure missing or invalid`
4. Ensure no downstream layer overwrites these values.

### Acceptance criteria

- All 29 unresolved rows are excluded from similarity clustering.
- No row has `Domain_Status = Not assessable` with `Decision = PASS`.
- No unresolved row receives a numeric Tanimoto or scaffold-coverage value.
- Structured rows retain their existing cluster/domain outcomes.
- Cluster sizes and representative scaffolds for structured rows are computed
  only from the eligible structured population.

### Tests

- Unit test with one valid and one invalid SMILES.
- Regression test against all 219 fixture rows.
- Assertion that unresolved rows cannot enter the distance matrix.

## Stage 2 — Normalize missing and no-finding values

### Implementation status (2026-09-03)

Validated against the frozen 219-row fixture: valid rows retain explicit
no-finding values, unresolved rows use `Not assessable`, required status fields
are populated, and numeric zero descriptors remain preserved as numeric values.

### Work

1. Apply one export normalization function to every compound sheet.
2. Use explicit values:
   - `None` for a valid structure with no alert/toxicophore finding.
   - `Not assessable` for a missing or invalid structure.
   - `Unavailable: <reason>` for a failed external dependency.
3. Preserve numeric zero values (`0`, `0.0`) as real numbers.
4. Do not convert valid identifiers, statuses, or concentrations to empty
   strings.

### Acceptance criteria

- No silent blanks remain in the required structural evidence fields.
- A zero descriptor is distinguishable from missing data.
- Toxicophore and alert no-findings are displayed as `None`.
- Unresolved rows show `Not assessable` where appropriate.
- Valid no-finding structural rows show literal `None`, not an empty cell.

### Tests

- Test values `0`, `None`, `NaN`, and an external-error string.
- Verify numeric descriptor types after Excel export and reload.

## Stage 3 — Persist B1b state and evidence bundles

### Implementation status (2026-09-01)

Implemented for the local progressive workflow: the API now issues a `run_id`,
stores compound snapshots and pandas evidence bundles server-side, merges
substage artifacts, restores them before B2 formatting, exposes a recovery
endpoint, and expires abandoned state after four hours. The browser persists
the run ID and restores the current Phase B state after refresh.

### Work

1. Introduce a versioned, server-side Phase B run state keyed by a generated
   `run_id`. The serializable state envelope contains:
   - `run_id`
   - `compounds`
   - `cluster_evidence_matrix`
   - `endpoint_evidence_ledger`
   - `chemont_classifications`
   - `uncertainty_register`
   - stage completion flags
   - phase trace
2. Return the `run_id`, stage summary, trace, and compact row status from every
   B1b substage; keep the full evidence bundles in the local run state.
3. Load the prior state by `run_id` for the next substage; do not rely on
   pandas `DataFrame.attrs` across HTTP requests.
4. Merge each substage's fields into the prior state rather than replacing
   the whole bundle.
5. Allow retry of a failed substage without rerunning completed substages.
6. Add state versioning and an expiry/cleanup policy for abandoned local runs.

### Acceptance criteria

- B1b-1 → B1b-2a → B1b-2b → B1b-3 → B1b-4 preserves all previously calculated fields.
- The state envelope survives JSON serialization and deserialization.
- Evidence bundles remain available to B2.
- A failed optional stage does not remove prior successful data.
- A browser refresh can recover the current run from its `run_id`.

### Tests

- Round-trip serialization test for a one-row state.
- Sequential four-stage test with one external stage forced to fail.
- Assert that the final state still contains the successful earlier bundle.

## Stage 3.5 — Reconcile identity, review reason, and clustering eligibility

### Objective

Remove stale manual-review messages while preserving the distinction between
computationally usable surrogate structures and toxicologically confirmed
identities. This stage must not change valid descriptor or cluster calculations.

### Work

1. Add a canonical identity category produced by identity resolution:
   `Resolved`, `Surrogate`, or `Unresolved`.
2. Preserve the identity provenance fields: resolution source, surrogate reason,
   structure used for computation, and identity confidence.
3. Derive `Manual Review Flag` and `Manual Review Reason` from the canonical
   category and classification outcome:
   - `Unresolved` → `No valid standardized or original SMILES could be parsed.`
   - `Surrogate` → `Surrogate structure used; identity uncertainty requires toxicologist review.`
   - valid structure with no priority class → `No priority class matched the detected feature profile.`
   - resolved and classified → blank review reason and flag false.
4. Keep clustering eligibility based on parseable structure, not identity
   confidence: valid surrogate structures may be clustered; unresolved rows may
   not.
5. Add an export-boundary consistency validator as a safety net. It may repair
   only stale combinations or emit a validation failure; it must not invent
   scientific classifications.
6. Apply the same reconciliation before the single-assessment response is
   rendered so UI and Excel use identical values.

### Acceptance criteria

- No parseable clustered row contains the no-structure review reason.
- Every surrogate row retains a review flag and uses the surrogate reason.
- Every unresolved row retains the no-structure reason, `REVIEW`, and
  `Not assessable` domain/cluster values.
- Resolved, successfully classified rows have no manual-review reason.
- Descriptor values, cluster IDs, cluster sizes, and domain outcomes are
  unchanged by the reconciliation.
- UI and regenerated Excel values match for all five identity/review fields.
- The 43/147/29 sheet distribution is compared with the approved baseline and
  any intentional change is recorded.

### Tests

- Unit fixtures for Resolved, Surrogate, and Unresolved rows.
- Regression assertion for the 36 currently contradictory surrogate rows.
- Mixed valid/invalid clustering test proving eligibility is unchanged.
- UI/API versus Excel parity test for identity and review fields.
- Fresh 219-row workbook audit confirming zero stale no-structure reasons on
  clustered rows.

## Temporary operating mode — Deterministic local workflow

### Purpose

Provide a stable validation mode while ChemOnt and EPA CTX connectivity is
intermittent. This mode must preserve the scientific grouping logic and must
not silently treat unavailable external evidence as a negative finding.

### Rules

1. Disable ChemOnt and EPA CTX explicitly for the run; do not leave their
   clustering/evidence layers partially configured.
2. Set the ChemOnt soft-layer weight to `0` when no usable ChemOnt signatures
   are available, so an empty layer cannot distort or rescale distances.
3. Continue local deterministic processing for structures, descriptors,
   functional groups, topology, toxicophores, alerts, clustering, and domain.
4. Mark external fields as `Not requested` with a source/version status rather
   than blank or inferred values.
5. Add an uncertainty note to the compound and cluster evidence stating that
   external taxonomy and EPA endpoint evidence were not included.
6. Keep the external stages re-enableable later as optional enrichment without
   changing the local classification contract.

### Acceptance criteria

- A complete local run finishes without ChemOnt/EPA network calls.
- Primary classification, descriptors, structural alerts, cluster IDs, cluster
  sizes, and local domain calculations remain populated.
- ChemOnt/EPA unavailability cannot produce a 500 or block B2 export.
- Re-running the same fixture in local mode is deterministic.
- A comparison report records any cluster changes when external layers are
  later enabled.

## Phase A → Phase B structure carry-forward remediation

### Audit findings (2026-09-03)

The reviewed workbook contains 219 rows and all required columns, but it
regresses from the frozen baseline: structured rows fall from 190 to 172 and
unresolved rows rise from 29 to 47. Eighteen rows have valid standardized
SMILES in the approved Phase A workbook but are emitted as `Tier 3: Class
Fallback Required (No Structure)`. The affected rows include valid structures
and surrogate/homologue structures. Their downstream descriptors, alerts,
topology, and clustering are therefore lost. This is an identity/structure
handoff regression, not an EPA CTX/ChemOnt or clustering-policy issue.

### Implementation plan

1. Treat nonblank, parseable Phase A structure fields as authoritative during
   Phase B. Do not re-query PubChem for a row that already has a valid
   `Standardized SMILES`.
2. Preserve the incoming `Original SMILES`, `Standardized SMILES`, `InChIKey`,
   `InChI`, `Resolved Name`, and `Identity Status` before any optional lookup.
3. If a lookup is attempted for a missing/invalid structure and fails, retain
   the preserved Phase A values rather than entering the exception path that
   blanks structure fields.
4. Keep parseable surrogate structures computationally eligible while retaining
   their toxicologist-review flag and surrogate explanation.
5. Restrict Tier 3 fallback to rows with no valid incoming structure and no
   successful resolution.
6. Add a source-precedence trace to the Phase B run state: `phase_a_workbook`,
   `user_supplied`, `pubchem_resolved`, `surrogate_mapping`, or `unresolved`.
7. Run the corrected handoff through B1a and confirm descriptors/alerts/topology
   are calculated before clustering; do not alter the clustering weights or
   unresolved-row safety gate.

### Acceptance criteria

- Structured count returns to 190 and unresolved count to 29 for the frozen
  219-row fixture.
- The 18 identified compounds retain byte/value-equivalent standardized SMILES
  and identity fields after B1a.
- No valid/surrogate row loses descriptors, alerts, topology, or cluster
  eligibility because a repeat lookup failed.
- True unresolved rows remain `Not assessable`, `Review Required`, and `REVIEW`.
- EPA CTX and ChemOnt remain deferred and do not affect the handoff decision.

### Regression tests

- Fixture comparison before/after Phase B handoff for all identity/structure
  columns.
- Explicit test covering the 18 affected compounds and at least one surrogate.
- Simulated PubChem failure with a valid incoming SMILES; assert preservation.
- Simulated missing incoming structure; assert Tier 3 fallback.
- B1a workbook export count/sheet reconciliation against the Stage 0 baseline.

## Stage 4 — Make external enrichment bounded and observable

### Implementation status (2026-09-03)

EPA CTX and ChemOnt integration is intentionally deferred/removed from the
active release after the decision to avoid unreliable external dependencies.
Deterministic local processing remains the in-scope path. The adapters retain
bounded/error-safe behavior for possible future re-enablement, but no EPA 2a/2b
or ChemOnt service work is a release blocker.

### Work

1. Keep the current separate substages and split EPA internally:
   - B1b-1 ChemOnt
   - B1b-2a EPA DTXSID resolution
   - B1b-2b EPA hazard retrieval
   - B1b-3 ChEMBL
2. Add per-request connect/read timeouts and a defined stage-level deadline.
3. Add per-row progress counters and the last processed identifier to the
   phase trace.
4. Continue after row-level errors; record the error on that row.
5. Use cached results first and avoid repeating identical identifiers within a
   run.
6. For EPA CTX, split internally into DTXSID resolution and hazard retrieval
   with separate completion flags and retry behavior.
7. Return `None` when a valid compound has no ChEMBL match, and
   `Unavailable: <reason>` only when the request itself failed.

### Acceptance criteria

- One slow or unavailable service cannot hold the request indefinitely.
- Every row receives a terminal status for the stage.
- A full 219-row run completes within the approved local time budget or ends
  with a controlled timeout status.
- No external failure causes a 500 response for the whole workbook.
- The time budget, timeout values, retry count, and post-deadline behavior are
  recorded in the Phase B trace and configuration.

### Tests

- Mock timeout, HTTP 403, HTTP 404, empty response, and valid response.
- Verify cached rows do not issue duplicate network requests.
- Verify progress trace reports start, row counts, completion, and errors.

## Stage 5 — Populate evidence and taxonomy sheets from persisted state

### Implementation status (2026-09-03)

Implemented for local and persisted Phase B runs: local evidence and taxonomy
outputs are restored before B2, while EPA CTX and ChemOnt fields retain explicit
deferred/not-requested statuses. Their external completeness is not a release
requirement under the approved scope decision.

### Work

1. Export the consolidated `Compounds` sheet with deterministic local fields.
2. Export the complete cluster evidence matrix generated by deterministic
   clustering when available.
3. Export the uncertainty register when present.
4. Use placeholder rows only when the stage was genuinely not run; label them
   with `Stage not run`, not `No data`, when distinguishable.

### Acceptance criteria

- `Compounds` contains every input row exactly once.
- Cluster Evidence Matrix contains actual evidence rows for every eligible
  cluster after deterministic clustering.
- No completed stage leaves a one-row placeholder sheet.

### Tests

- Count reconciliation between compound rows and ChemOnt rows.
- Verify EPA ledger counts against row-level EPA status values.
- Verify every exported matrix row has a cluster ID and evidence status.
- Verify placeholder sheets are permitted only when their stage is explicitly
  marked `not_run`.

## Stage 6 — Add the B2 export validation gate

### Implementation status (2026-09-03)

Implemented for the local workflow: B2 now validates required compound fields,
unresolved-row safety, retired Cramer/TTC absence, and Phase B prerequisites
before formatting Excel. Blocking structural failures return HTTP 409 with an
explicit validation error list; B2 is permitted after B1a, while optional
external stages remain non-blocking.

### Work

Before generating Excel, validate:

1. Required compound columns exist.
2. Row counts reconcile across compound sheets and audit sheets.
3. No unresolved row has `PASS` or `INSIDE DOMAIN`.
4. No required completed-stage bundle is missing.
5. External failures are explicit, not blank.
6. Retired Cramer/TTC columns are absent.
7. Phase trace records each completed substage.

If a structural logic invariant fails, block B2 and show the exact failing
rows. If only an external service is unavailable, allow export with explicit
statuses.

### Acceptance criteria

- Invalid workbook state cannot be downloaded as a final workbook.
- External unavailability does not prevent a controlled export.
- The validation result is visible in the UI and stored in the phase trace.
- B2 requires each declared prerequisite stage to be either `complete` or
  explicitly `skipped`; it cannot infer completion from partially populated
  columns.

### Tests

- Deliberately inject an unresolved `PASS` and confirm B2 is blocked.
- Deliberately remove ChemOnt state after B1b-1 and confirm B2 explains the
  missing bundle.
- Inject an EPA timeout and confirm B2 remains exportable with status text.

## Stage 7 — UI workflow and layout

### Implementation status (2026-09-03)

Implemented in the localhost UI. Phase B is presented as ordered checkpoints
with a visible status/progress area. Phase A remains independently downloadable
and a Phase A workbook can be loaded directly for Phase B testing. After a
successful checkpoint, its button remains disabled and only the next valid
checkpoint is enabled. On failure, the failed checkpoint is the only retryable
button and prior completed results are retained. The `run_id` is persisted in
browser storage and server-side state can restore the checkpoint state after a
refresh. B2 unlocks after B1a and the backend validation gate passes. Optional
B1b stages may add evidence but do not block local export.

Known limitation: the EPA 2a/2b split and full external-service progress
telemetry remain Stage 4 work; the UI currently exposes the implemented B1b
substage sequence and bounded failure/retry behavior.

### Work

1. Show B1a, B1b-1, B1b-2a, B1b-2b, B1b-3, B1b-4, and B2 as a single ordered
   workflow (responsive layout may wrap the buttons).
2. Keep completed buttons visibly inactive/disabled in the normal flow; enable
   only the failed stage for retry.
3. Unlock only the next valid stage.
4. Show progress, row count, current dependency, and terminal status.
5. Preserve the Phase A download button independently.
6. Display a concise failure card with retry for the failed substage.

### Acceptance criteria

- The user can identify exactly which stage is running.
- Completed stages cannot be accidentally rerun from the normal path.
- A failed stage can be retried without losing previous results.
- The page never shows a completed stage as pending after its response returns.
- A refresh restores the stage state from `run_id` and does not silently reset
  completed stages.

### Tests

- Browser test of the complete button sequence.
- Refresh/reload test during and after each stage.
- Error-and-retry test for EPA CTX.
- Responsive-layout test confirms all buttons remain readable when wrapped.

## Stage 8 — Performance, security, and release validation

### Implementation status (2026-09-03)

Release validation is substantially complete for the deterministic/local path.
The full offline regression suite was run in the bundled build environment; the
distance-engine compatibility regression and evidence-ledger placeholder check
were corrected and their targeted tests pass. The onedir executable was rebuilt
at `dist/xchem_engine/xchem_engine.exe` and a fresh EXE startup smoke test
returned HTTP 200 on port 8055. The packaged build includes the cache-busted UI
assets and the Phase B checkpoint workflow.

Security review confirms that API credentials are passed through runtime
configuration and are not part of workbook fields or UI traces. External
enrichment remains bounded by the configured stage deadline. One legacy
end-to-end assertion still expects every row—including unresolved review-only
rows—to receive a `Cluster_*` identifier; that expectation conflicts with the
Stage 1 safety rule and is retained as a known test-contract update for final
sign-off rather than weakening the production behavior.

### Work

1. Run the complete 219-row workflow with all external services available,
   cached, and unavailable.
2. Verify that API keys are never written to traces, workbook cells, or logs.
3. Verify local run-state expiry and cleanup.
4. Rebuild the executable and confirm it contains the updated backend and
   cache-busted UI assets.
5. Run the EXE smoke test through B1a, every B1b substage, and B2.

### Acceptance criteria

- The 219-row run completes within the configured budget or terminates with a
  controlled, resumable status.
- Secrets are masked everywhere outside the request header.
- The EXE exposes the same stage buttons and produces the same workbook schema
  as the source application.

### Tests

- 219-row performance test with timing per stage.
- Secret/log scan for API-key leakage.
- EXE endpoint and workbook smoke test.

## Stage 9 — Final workbook verification and sign-off

### Implementation status (2026-09-03)

The authoritative 219-row baseline was re-opened and reconciled, and the
required workbook sheets and safety statuses are documented in
`stage9_final_verification.md`. Localhost and the rebuilt EXE are reachable and
the targeted release regression suite passes. Final sign-off remains pending a
bounded full post-Phase-B regeneration/render and updating one stale test that
expects unresolved review-only rows to receive cluster IDs.

### Work

1. Generate a fresh workbook from the corrected application.
2. Reconcile all sheet names, row counts, headers, statuses, and key values.
3. Render each sheet for visual review.
4. Record the final verification results and known external limitations.
5. Only then mark the implementation complete.

### Final sign-off criteria

- 219 compound rows are accounted for.
- 190 structured rows retain full deterministic evidence.
- 29 unresolved rows are clearly review-only and not clustered/passed.
- ChemOnt and EPA CTX are absent from the active export; ChEMBL, when retained,
  is explicit and non-blocking, and cluster evidence remains traceable.
- No placeholder audit sheet remains after its corresponding stage completes.
- B2 export is repeatable and bounded.
- The executable and localhost builds pass the same final verification record.

## Export schema consolidation — single compound sheet and ordered evidence

### Status (2026-09-03)

Implemented. ChemOnt/EPA CTX fields and tabs are excluded from the active
export contract, all compound populations are emitted on one deterministic
`Compounds` sheet, and XAI fields are last. B2 can run after deterministic B1a;
optional B1b stages remain additive and non-blocking.

### Objective

Make the workbook easier to review by presenting all compounds in one sorted
sheet, while retaining confidence and identity status as explicit columns.

### Work

1. Replace the separate `High Confidence`, `Medium-Low Analogs`, and
   `Unresolved (Review)` output sheets with one `Compounds` sheet.
2. Sort the compound sheet by valid numeric `Cluster ID`, then cluster size,
   confidence, and finally unresolved/review-only rows.
3. Preserve `Confidence`, `Identity Status`, and `Manual Review Flag` so the
   former sheet categories remain filterable without splitting the dataset.
4. Keep `ChEMBL_Max_Phase` and `ChEMBL_Targets` as optional supplementary
   fields. Populate `None` or `Unavailable: <reason>` explicitly; never allow
   ChEMBL latency or failure to block clustering or export.
5. Move `XAI_Status` and `Read_Across_Justification` to the final columns of
   the `Compounds` sheet. In AI-isolated mode, use explicit not-run values.
6. Keep `Raw Data`, `Workbook Contract`, `Cluster Evidence Matrix`, and
   `Uncertainty Register` as the visible supporting tabs.
7. Omit confidence-specific tabs and external ChemOnt/EPA tabs from the final
   workbook. Retain retired fields only in the contract's audit record.
8. Keep policy tabs hidden only when they are required for reproducibility;
   otherwise consolidate their metadata into `Workbook Contract`.

### Acceptance criteria

- Exactly one compound output sheet contains all input, identity, structure,
  classification, descriptor, alert, clustering, domain, ChEMBL, and XAI fields.
- No compound row is lost or duplicated during the merge.
- Cluster sorting is deterministic; unresolved rows appear at the end and keep
  their review-only values.
- XAI fields are the final columns and contain explicit AI-isolated statuses.
- ChEMBL `None`/`Unavailable` values are distinguishable from missing columns.
- ChemOnt/EPA CTX columns and tabs are absent from the active export.
- Supporting audit tabs remain readable and do not contain duplicate compound
  populations.

### Tests

- Compare pre/post merge row counts and unique identity keys.
- Verify cluster ordering for numeric, non-numeric, and unresolved cluster IDs.
- Confirm confidence and review filters reproduce the former three-sheet
  populations.
- Simulate ChEMBL timeout and verify export still completes with explicit
  unavailable values.
- Inspect the final header order and assert XAI fields are last.
- Verify no `ChemOnt`, `EPA_CTX`, or confidence-specific sheet is present in
  the final workbook.

## Stage 10 — Embedded SME review guide and workbook dashboard

### Status

Planned. The reviewed source instruction is maintained in
`docs/SME_WORKBOOK_REVIEW_INSTRUCTIONS.md`. This stage embeds a concise,
read-only version in each final workbook and adds an auditable management and
SME overview dashboard. It does not change identity resolution,
classification, clustering, or any underlying compound result.

### Objective

Make every final workbook self-explanatory and reviewable without a separate
document, while providing a reliable high-level view of the source population,
unique chemicals, identity-resolution outcome, and structural grouping.

### Workbook order

The final workbook will place the following sheets at the beginning:

1. `SME Review Guide`
2. `Dashboard`
3. `Workbook Contract`
4. `Raw Extraction (Pre-Dedup)`
5. `Summary`
6. `Detailed Analysis`
7. `Cluster Summary`

Existing evidence/support sheets follow these in their current logical order.

### Work item 10.1 — SME Review Guide worksheet

1. Add a read-only `SME Review Guide` worksheet to every final export.
2. Base its content on `docs/SME_WORKBOOK_REVIEW_INSTRUCTIONS.md`, condensed
   for workbook reading and formatted with section headers, short decision
   tables, and usable column widths.
3. Include the required review sequence: `Workbook Contract` → `Raw Extraction
   (Pre-Dedup)` → `Summary` → `Detailed Analysis` → `Cluster Summary`.
4. Define `Computed`, `Preview only`, `Retired`, `Not assessable`,
   `Unclassified`, and `SME Review Required` exactly as stated in the workbook
   contract.
5. Include recommended Summary/Detailed Analysis sorting and filtering,
   identity/classification/cluster review checks, and accept / accept with
   limitations / do not use / follow-up outcomes.
6. State prominently that PRISM structural grouping does not establish
   endpoint equivalence, read-across acceptability, exposure equivalence,
   safety, or a final toxicological conclusion; SME review remains required.
7. Include a compact decision-record template: reviewer, date, intended use or
   endpoint, outcome, rationale, and required follow-up.

### Work item 10.2 — Dashboard worksheet

1. Add a read-only `Dashboard` worksheet after `SME Review Guide`.
2. Populate all dashboard values from the exported data frames during the same
   export operation. Do not use manually maintained workbook formulas or
   static counts.
3. Show clearly labelled KPI values for:
   - raw extracted source rows before deduplication;
   - final deduplicated chemical records;
   - unique reported CASRNs;
   - unique resolved standardised structures;
   - structure-resolved versus unresolved/manual-review records;
   - records eligible for grouping;
   - compatibility groups;
   - clusters, singleton clusters, and clustered compounds; and
   - records with alerts, property/compatibility limitations, or SME review
     requirements.
4. Add compact, data-driven visual summaries for:
   - compounds by compatibility group;
   - compounds by corrected chemical class or primary functional group;
   - resolved versus unresolved/manual-review status; and
   - the largest clusters and singleton-cluster count.
5. Include precise Dashboard definitions:
   - **Raw data:** every extracted source row before duplicate handling.
   - **Unique chemical records:** deduplicated output records, with CASRN and
     standardised-structure counts shown separately.
   - **Grouping eligible:** usable standardised structures that passed the
     compatibility-gating prerequisites.
   - **Clustered:** a grouping-eligible structure assigned to a structural
     cluster; singleton clusters are counted separately.
6. Display the generated workbook timestamp, PRISM version/rule version where
   available, and a clear statement that the dashboard is an overview; the
   detailed evidence sheets are authoritative.
7. Use accessible, print-friendly formatting with descriptive chart titles and
   no chart that implies toxicological risk or ranking.

### Work item 10.3 — Workbook Contract header clarification

1. In the `Workbook Contract` worksheet, rename the existing `Column` header
   to `Column Headers - Detailed Analysis sheet`.
2. Retain all field names, definitions, status values, ordering, and workbook
   logic exactly as they are; this is a label-only clarification so reviewers
   understand that the listed fields correspond to `Detailed Analysis`.

### Acceptance criteria

- Every final workbook contains `SME Review Guide` and `Dashboard` in the
  stated sequence.
- The guide is legible without an external document and correctly states the
  SME responsibility and review boundary.
- Dashboard counts reconcile exactly to `Raw Extraction (Pre-Dedup)`,
  `Summary`/`Detailed Analysis`, and `Cluster Summary`.
- Raw-row, final-record, unique-CASRN, and unique-standardised-structure counts
  are explicitly distinct and never conflated.
- Unresolved structures are visible in the Dashboard but are not counted as
  grouping eligible or structurally clustered.
- Charts are generated solely from the exported data and remain correct for a
  different input report.
- The `Workbook Contract` field-list header reads `Column Headers - Detailed
  Analysis sheet`, with no changes to the associated contract data.
- Existing output values, cluster IDs, classifications, and evidence sheets are
  unchanged except for the addition of the two new worksheets.

### Tests

- Generate a workbook from a known fixture and independently reconcile every
  Dashboard count against the relevant worksheet.
- Generate a workbook with duplicates, aliases, unresolved records, singleton
  clusters, and multiple compatibility groups; verify all dashboard definitions
  and counts.
- Render `SME Review Guide` and `Dashboard` to verify readability, page layout,
  frozen panes, widths, and chart labels.
- Regression-test that dashboard generation neither modifies Detailed Analysis
  data nor changes clustering/classification results.
