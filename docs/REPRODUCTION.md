# Execution and input guide

## Release checks and environments

Existing implementations are exported with source hashes and path bindings. No accepted fit, coordinate QC, RepeatMasker, alignment, tree, STRIDE, HOG or orthology stage was run during curation.

Expanded-effect adapters retain the original identity Py3.12.14/NumPy2.3.5/SciPy1.16.3/one_BLAS_thread. environment/effects.txt records package pins. Root requirements are the separately accepted synthetic environment. Discovery also uses pandas; plotting uses matplotlib; optional saved-PDF QA uses pypdf; the original Figure 2 renderer additionally uses reportlab and pypdfium2. Exact historical builds for these additional packages are not asserted here.

Set OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1. Use fresh output directories, separate from original results and once ledgers.

## Offline inventory checking

Put the three archives in config/catalog-flat.json under one chosen directory:

```sh
python scripts/verify_evidence.py --catalog config/catalog-flat.json --archive-root /path/to/archives --mode manifest --report fresh-evidence-report.json
```

Manifest mode checks root manifests, member names/sizes and dependencies. It follows no historical source paths and executes no analysis. Content verification is a separate explicit mode.

## Original three-species analysis

analysis/discovery/run_comparison.py expects the accepted Stage-A table all-1896-triplets-two-scales-input-statistics.tsv.gz and execution-scope.json with frozen identities. OUT is the script directory and STAGE its parent. For independent reproduction copy the script into a NEW Stage-B-main-comparison-v1 directory in a separate evidence view; supply the original relative inputs and execution scope, then inspect its main entrypoint before executing. The consumer writes exclusive outputs and a once claim. Source genomes, candidate selection and unretained historical parameters are not reconstructed by this export.

## Expanded effects and sensitivity

v12 contains frozen panels, accepted results and prepared input packets. Not every logical state is numerical. A prepared sensitivity packet contains status PREPARED_NUMERIC_INPUT, main_panel, observations, variant, narrow_rows, numeric_fingerprint and effective input identities. Typed-NA records are reported directly.

The retained consumer uses this exact numerical call after binding validation:

```python
import sys
sys.path.insert(0, "analysis/effects")
import backend

# inp is a decoded, hash-verified prepared packet from the evidence archive.
if inp["status"] != "PREPARED_NUMERIC_INPUT":
    result = inp
else:
    kw = {"narrow_rows": inp["narrow_rows"]} if inp["variant"] == "A25" else {}
    result = backend.spatial_jackknife(
        inp["main_panel"], inp["observations"], variant=inp["variant"], **kw)
```

This is a future reproduction recipe, not executed during upload. Use original_caller.decode_result for retained bytes/float tags and original_caller.encode_result or json_safe.json_safe for outputs. compute_effect_packet_v2.build(cfg, out, check, events) retains the complete binding/claim contract. Its code and source descriptors must point to a separately constructed evidence view. Historical receipt paths provide provenance, not current file authority.

cross_scale_inputs supplies off-diagonal and fixed-member scale inputs. joint_scores_v2 consumes saved accepted coefficients, global blocks and fixed policy. Full-S35 unavailable summaries remain unavailable; dropping source-NA species changes the estimand and does not repair the original result.

## Spatial and controlled-distance wrappers

Put the repository root and analysis/spatial, analysis/location, analysis/effects on PYTHONPATH as needed. Pure interval and geometry modules have no source dependencies. Historical execution wrappers need explicit external roots, listed in config/paths.example.json:

- ANRL_PROJECT_ROOT: separate original project evidence view with named source/legacy inputs.
- ANRL_ORIGINAL_STATS_ROOT: statistics view containing Stage-B-main-comparison-v1 accepted summaries and saved backgrounds.
- ANRL_R3_ROOT: round3 evidence view with geometry, offsets, tracks and accepted outputs under their original logical names.
- ANRL_R4_ROOT: round4 view, with round3 as the parent expected by historical wrappers.

Cloning this repository alone does not create those views. Archive manifests retain exact member identities and original directory names; rebind them explicitly. Source-preparation scripts additionally need the named original sources. Reuse accepted tracks where available. Freeze geometry and offsets before evaluating response tracks, preserve the fixed panel and all allowed zero/duplicate shifts, and do not replace inadmissible shifts by nearby values.

## Displays

v12 figure functions consume accepted compact dictionaries and return SVG. figure3_svg renders the frozen tree order without inferring a topology. draw_v12_figures_once additionally expects source descriptors and explicit code bindings. Figure6/7 consume saved accepted statistics. No figures were regenerated during upload.

