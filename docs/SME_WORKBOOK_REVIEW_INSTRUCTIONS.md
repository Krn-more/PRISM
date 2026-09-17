# PRISM SME Workbook Review Instructions

## Purpose and review boundary

PRISM provides a traceable, structure-based organisation of extracted chemical records to support subject-matter-expert (SME) review. It does **not** establish toxicological endpoint equivalence, biological equivalence, safety, or a final read-across decision. Those decisions remain the responsibility of the SME.

Use this instruction to review each final PRISM Excel workbook consistently. The intended workbook sequence is:

1. `SME Review Guide` (when included in the exported workbook)
2. `Workbook Contract`
3. `Raw Extraction (Pre-Dedup)`
4. `Summary`
5. `Detailed Analysis`
6. `Cluster Summary`

## 1. Read the Workbook Contract first

Before reviewing results, read `Workbook Contract` to confirm the meaning of the exported fields, status values, calculation boundaries, and traceability expectations.

Interpret status values as follows:

| Status | Meaning for review |
|---|---|
| `Computed` | PRISM calculated the value from available data; it remains subject to scientific review. |
| `Preview only` | Informational/intermediate content, not intended as a final decision field. |
| `Retired` | Historical support field retained only for compatibility or traceability; do not use it in a new decision. |
| `Not assessable` | PRISM could not make a defensible calculation from the available information. This is not a negative result. |
| `Unclassified` | No sufficiently supported structural classification was available. Review the identity and source evidence. |
| `SME Review Required` | A defined limitation, uncertainty, or boundary requires SME judgement before the record is used. |

## 2. Verify the raw extraction record

Use `Raw Extraction (Pre-Dedup)` as the source audit trail. For priority records, check:

- compound name and CASRN against the analytical report;
- source, page, and table reference;
- reported concentration and unit; and
- whether repeated records are genuine duplicates, aliases, or separate reported findings.

This worksheet confirms what was extracted before deduplication. Do not use it alone to make structural grouping decisions.

## 3. Triage the Summary worksheet

Use `Summary` for the initial review queue. Apply the following sort order:

1. `Manual Review Flag` — review flagged records first;
2. `Uncertainty Summary` — highest concern first;
3. `Cluster ID` — ascending;
4. `Cluster Size` — largest first within cluster; and
5. `Compound Name` — ascending.

Prioritise records with missing or questionable structures, identity-resolution limitations, high reported concentrations, singleton clusters, small/borderline clusters, inconsistent alerts or ionisation, and property-outlier flags.

## 4. Confirm identity and structure

For each prioritised compound, use `Detailed Analysis` to check:

- `Compound Name` and `CASRN`;
- `Original SMILES`, `Standardized SMILES`, `InChIKey`, and `InChI`;
- identity-resolution status and source evidence; and
- report source/page/table reference.

### SME decision

- If identity and structure are credible, continue to classification and cluster review.
- If structure is unavailable, invalid, or inconsistent with the reported identity, retain the record as unclassified/manual review. It must not be used as structural analogue evidence.
- If an identity correction is needed, correct the controlled input/source record and reprocess; do not edit downstream classification or cluster fields in isolation.

## 5. Review deterministic chemical classification

Assess whether the classification is chemically plausible from the standardised structure. Review:

- `Corrected Chemical Class`;
- `Primary Functional Group` and `Secondary Functional Groups`;
- `Taxonomy Path` and `Taxonomy Hierarchy`;
- `Detected Feature Profile`;
- `Toxicophore Profile`; and
- `All Structural Alerts` / `Alerts`.

Accept a classification when it is chemically reasonable and adequately specific for its intended use. A broad but scientifically defensible class may be retained with a documented limitation. Flag incorrect or insufficiently specific classifications for controlled rule refinement.

## 6. Review the Cluster Summary

Review clusters one at a time in `Cluster Summary`. Check:

- `Compatibility Group` (or legacy `Cluster Family`);
- `Cluster ID` and `Cluster Size`;
- representative compound and standardised SMILES;
- member chemical classes and primary-functional-group profile;
- minimum/median pairwise Tanimoto similarity and nearest-neighbour range;
- scaffold, ionisation, structural-alert, and property-compatibility profiles;
- compatibility-gate reason and membership rationale; and
- stability assessment.

### How to interpret the grouping

`Compatibility Group` is the broad structural comparison space. It prevents unsuitable comparisons, for example a peroxide with a siloxane, before clustering is attempted.

`Cluster ID` is a membership identifier only. `Cluster01` is not more important, hazardous, or similar than `Cluster02`; the number is not a score.

Within a compatible group, PRISM represents standardised structures using Morgan/ECFP4 fingerprints (radius 2, 1,024 bits), compares them with Tanimoto similarity, applies a 0.60 minimum similarity threshold using Butina clustering, and then performs complete-linkage refinement. Complete linkage reduces chain clustering: A and B, and B and C, cannot remain in one cluster where A and C are not sufficiently similar.

## 7. Make an SME cluster decision

Use one of the following outcomes for each cluster:

| SME outcome | When it is appropriate |
|---|---|
| Accept for analogue review | Members have defensible shared structural features, adequate cohesion, and no material identity or compatibility concern. |
| Accept with limitations | The cluster is structurally coherent but small, borderline, or has incomplete descriptor evidence. State the limitation. |
| Do not use for analogue review | Functional groups, scaffold, ionisation, alert profile, or relevant properties differ materially, or cohesion is weak. |
| SME review required / follow-up | Identity is uncertain, the class is provisional, the cluster is a singleton/borderline, or endpoint relevance is unresolved. |

An accepted structural cluster supports analogue consideration only. It does not by itself support read-across, endpoint equivalence, exposure equivalence, or a safety conclusion.

## 8. Use individual membership rationale as supporting evidence

In `Detailed Analysis`, use `Cluster_Membership_Rationale` to understand the assigned cluster for one compound. It may report nearest in-cluster neighbours and Tanimoto values, representative-scaffold coverage, ionisation and toxicophore consistency, property-outlier assessment, and domain status.

This is deterministic evidence for structural membership. It should be read together with the Cluster Summary and the intended toxicological endpoint.

## 9. Handle AI cluster reasoning appropriately

Where `AI_Cluster_Reasoning` is present, treat it as a readable interpretation of the deterministic evidence supplied to the model. Verify that it agrees with the corresponding cluster fields and does not introduce unsupported chemical, biological, exposure, or endpoint claims.

The deterministic cluster evidence remains authoritative. AI reasoning is supporting audit text only. An SME must review it before analogue selection, read-across, or endpoint-specific assessment.

## 10. Document the SME decision

For each reviewed compound or cluster, record:

- reviewer name and review date;
- scope of review and intended endpoint/use case;
- outcome: accepted, accepted with limitations, not suitable, or follow-up required;
- scientific rationale and key supporting evidence;
- identity, classification, or clustering issue requiring correction; and
- any reprocessing or controlled rule-refinement request.

## Final review checkpoint

Before using PRISM output in an assessment, confirm that identities are reliable, classifications are chemically plausible, compatibility groups and clusters are defensible, material outliers are addressed, and all endpoint-specific scientific judgement is documented by the SME.
