# Animal neighborhood repeat load

Research code for repeat-load correspondence around homologous animal genes. The original three-hemichordate comparison and conditional 54-species extension are represented by the retained implementations below. Code is separated from versioned data and evidence archives.

## Implementation coverage

| Component | Code |
| --- | --- |
| Original discovery comparison, structural adjustment, cyclic backgrounds and chromosome influence | analysis/discovery/run_comparison.py |
| HOG-to-member bridge, annotation indexing, complete gene/all-exon evidence and genome callability | analysis/coordinates/ |
| Frozen membership/weights and A50/K50/A25 coefficients/deletions | analysis/effects/backend.py, original_caller.py, json_safe.py, effect_inputs.py |
| Own/common mapping and crossed definition/scale sensitivity | analysis/effects/cross_scale_inputs.py, compute_effect_packet_v2.py, contract_checks.py |
| Selection comparison, Unknown weighting, nested coverage, gene composition and summary/cache adapters | analysis/summaries/ |
| Dependence scores, interval information gates and layer omission | analysis/effects/joint_scores_v2.py |
| True-bp physical backgrounds | analysis/spatial/ |
| Controlled-distance supplement | analysis/location/ |
| Integrated evolutionary composite and other displays | figures/ |
| Offline evidence checking | scripts/verify_evidence.py |
| Synthetic diagnostic | Root synthetic_zero.py and calibration-frozen-spec.md |

CODE_PROVENANCE.json records original hashes and path-only substitutions. All 70 exported research implementation files preserve source function/class ASTs. Historical dispatchers, credentials and unrelated research are excluded.

## Start here

Read [the execution guide](docs/REPRODUCTION.md), [input/output map](docs/INPUT_OUTPUT.tsv), [figure map](docs/FIGURES.tsv) and [upstream scope](docs/UPSTREAM_AND_SCOPE.md). The [archive catalog](config/data-archives.json) identifies versioned data. Zenodo is currently a restricted, unpublished draft; no public DOI or anonymous access is claimed.

Root requirements belong to the synthetic diagnostic. Original expanded effects use a different retained runtime identity; see environment/. This upload checks code identity, syntax, import closure and documentation. It does not rerun scientific analyses or certify complete cross-environment numerical reproduction.

## Scope and availability

All NA states and literal reasons remain. OG0000725 remains unresolved, with 154 genes in 41 species. Source-unavailable species stay in their original cohorts. The historical remote-job UNKNOWN remains recorded outside this repository and is not resolved by uploading code.

Author/affiliation details are in the manuscript. A project-code license has not yet been selected; this release does not assert a new redistribution license. External tools retain their upstream licenses.


## Additional reproducibility deposit

The supplementary deposition audit adds the original 1,896-group main input and saved outputs, complete accepted native HOG tables, original/resolved native gene trees, 54 aggregate orthology files, a pinned repository snapshot and the current manuscript. See docs/DEPOSITION_COMPLETENESS.md. Existing data archives and exceptions remain retained.
