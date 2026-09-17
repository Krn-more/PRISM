# Stage 1–3 representative benchmark expected results

Fixture: `stage1_to_3_representative_benchmark.csv` (13 non-confidential rows).
The CSV is the controlled input; this document describes review-relevant expected
outcomes. The automated test asserts these outcomes and a SHA-256 manifest.

| Fixture IDs | Expected Stage 1 | Expected Stage 2 | Expected Stage 3 review signal |
|---|---|---|---|
| B01–B03 | Calculated; neutral aromatic ether structures | `Inside` or review-only property limitation; legacy decision unchanged | Exposure remains incomplete at cluster level; biological evidence remains high pending endpoint review |
| B04 | Calculated; neutral same-scaffold member | Property-range/outlier review is evaluated against Cluster_A | Property differences require review before shared assessment |
| B05 | Calculated; neutral | Cluster_B is compared with B06 | Reactivity/domain evidence remains review-only |
| B06 | Calculated; aldehyde toxicophore/alert | `Borderline` due to divergent toxicophore/structural evidence | Reactivity and domain actions required |
| B07 | Calculated; likely acidic | Singleton is `Borderline` | Domain and missing exposure require review |
| B08 | Calculated; likely basic | Singleton is `Borderline` | Domain and missing exposure require review |
| B09 | Calculated; amphoteric | Singleton is `Borderline` | Domain and missing exposure require review |
| B10 | Calculated; permanently charged | Singleton is `Borderline` | Domain and missing exposure require review |
| B11 | Calculated; neutral | `Borderline` due to surrogate identity | Identity review action is required |
| B12 | Missing structure; unknown ionisation | `Not assessable` | High identity, structure, reactivity, and domain uncertainty |
| B13 | Invalid structure; unknown ionisation | `Not assessable` | High structure and domain uncertainty |

This benchmark verifies review evidence only. It must not be used to calibrate
or change cluster membership, legacy `Domain`, or `Decision`.
