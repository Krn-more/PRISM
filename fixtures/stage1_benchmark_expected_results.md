# Stage 1 Benchmark Expected Results

This fixture validates calculation behavior only. It does not validate or change
clustering, domain, read-across, endpoint, or safety conclusions.

| Compound | Structural evidence status | Expected ionisation indicator | Expected key evidence |
|---|---|---|---|
| Ethanol | Calculated | Neutral | HBD 1; HBA 1; ring count 0 |
| Acetic acid | Calculated | Likely acidic | carboxylic-acid structural motif |
| Triethylamine | Calculated | Likely basic | non-amide amine structural motif |
| Glycine | Calculated | Amphoteric | acidic and basic structural motifs |
| Tetramethylammonium | Calculated | Permanently charged | formal charge +1 |
| Benzene | Calculated | Neutral | aromatic ring count 1; HBD/HBA 0 |
| Phenethyl alcohol | Calculated | Neutral | HBD 1; HBA 1 |
| 4-vinylbenzaldehyde | Calculated | Neutral | aldehyde toxicophore; aldehyde alert expected |
| Unresolved | Missing structure | Unknown | descriptor fields null/not assessable |
| Invalid structure | Invalid structure | Unknown | descriptor fields null/not assessable |

The executable test suite is `test_structural_evidence.py`. Expected values are
asserted there and are versioned with the Stage 1 rules.
