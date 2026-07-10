# gnomad_chets — Project Instructions

## What goes in this file vs `CLAUDE.local.md`

This file (`CLAUDE.md`) is **committed** and visible to everyone touching the repo. It should hold:

- Project structure, pipeline architecture, and how the pieces fit together.
- Code conventions and style choices that apply repo-wide.
- Hail / pipeline gotchas that any contributor would hit.
- Generic invocation patterns (commands without machine-specific paths or cluster names).
- Coupling notes (e.g. browser-side interfaces that need to stay in sync).

`CLAUDE.local.md` is **gitignored** and per-developer. It should hold:

- Hardcoded paths (conda env activation, gcloud SDK, adjacent repos on this machine).
- Personal Dataproc cluster names, GCP project / region, and one-liner `gcloud` recipes.
- Personal scratch-bucket / output-postfix conventions.
- Workflow preferences (confirm-before-commit, etc.) — overrides of `~/.claude/CLAUDE.md` scoped to this repo.

**Claude: keep these files current.** When you learn something genuinely useful that future you would want to know — a recurring gotcha, a new pipeline convention, a non-obvious file location, a constraint that bit us — add it to whichever file fits. Be concise: prefer one tight bullet over a paragraph, and remove or compress entries that have become stale. Don't bloat either file with conversation-specific or one-off details (those belong in commit messages or `analysis/`).

## What this repo does

`v2/` contains the production v2 phasing / variant-cooccurrence pipeline ([Guo, Francioli et al., *Nat Genet* 2024](https://www.nature.com/articles/s41588-023-01608-3)). The phasing logic in `phasing.py` is the historical reference; `compute_gnomad_phase.py` is the orchestration entry point.

`v4/` contains the in-progress port to gnomAD v4 plus the in-trans observed-vs-expected (in-trans-OE) depletion feature. The main scripts:

- `v4/create_vp_list.py` — builds the sites HT, variant-filter HT, filtered VMT, and variant-pair list from gnomAD v4 exomes (the upstream steps).
- `v4/compute_vp_counts.py` — consumes the variant-pair list and produces the per-pair genotype-counts table (multi-step; see below).
- `v4/run_in_trans_oe.py` — consumes the gt-counts table, computes per-candidate in-trans depletion against ClinVar P/LP partner sets.

`analysis/` holds one-off characterization scripts and reports.

## Conda environment

Use the `gnomad_chets` conda env (cloned from `hail`). This **overrides** any global preference for the `hail` env in `~/.claude/CLAUDE.md`.

(Activation path is per-machine — see `CLAUDE.local.md`.)

## v4 pipeline architecture

The pipeline is a chain of flags in dependency order. The upstream group runs in `create_vp_list.py`; the counts group in `compute_vp_counts.py`:

| script | flag | input | output |
|---|---|---|---|
| create_vp_list | `--preprocess-sites-ht` | gnomAD VEP / freq / final-filter | `sites.{postfix}.ht` |
| create_vp_list | `--create-variant-filter-ht` | sites HT | `variant_filter.{postfix}.ht` |
| create_vp_list | `--filter-vmt` | gnomAD VDS, filter HT | `filtered_vmt.{postfix}.mt` |
| create_vp_list | `--create-variant-pair-list-ht` | filtered VMT, filter HT | `variant_pairs.{postfix}.ht` |
| compute_vp_counts | `--encode-genotypes` | gnomAD VDS, pair list | `genotype_count_intermediates.{postfix}/...` (var_idx + encoded gt sets) |
| compute_vp_counts | `--build-variant-size-info` | encoded intermediates, pair list | `variant_size_info.{postfix}.ht` |
| compute_vp_counts | `--compute-counts-light` / `--compute-counts-heavy` | encoded intermediates, pair list, size-info | `.../counts_{light,heavy}.ht` |
| compute_vp_counts | `--combine-counts` | light + heavy counts | `variant_pairs.genotype_counts.{postfix}.ht` |

`--encode-genotypes` densifies only the pair-list variants out of the VDS in-memory (transient scratch checkpoint, no persisted dense MT) before encoding — there is no separate dense-MT step.

`run_in_trans_oe.py` then runs `annotate_pair_oe_terms` → `aggregate_oe_per_candidate` → output HT.

`phase_gnomad.py` runs `hl.experimental.haplotype_freq_em` on the gt-counts HT. By default it repartitions/checkpoints the input to balance the EM shuffle, but if the input's native partitioning is already good (e.g. a freshly-written HT from a Tier 3 post-hoc fix or from `--compute-counts-per-sample`) you can pass `--em-partitions 0` to skip the shuffle entirely. Bounding partition count also caps worst-case per-partition wall time if EM stalls (see the `haplotype_freq_em` gotcha below).

### Trio phasing (`trio_phasing.py`)

`v4/trio_phasing.py` consolidates the whole v2 trio/PBT analysis (`phase_by_transmission.py` + the `create_pbt_*` steps in v2 `create_vp_matrix.py` + `chet_utils` helpers + the gnomAD-comparison export) into one script (resources via `get_trio_phasing_resources`). It **reuses** the v4 counts machinery and EM phasing rather than duplicating: imports `create_variant_pair_ht` from `create_vp_list.py`, `create_variant_pair_filter_ht` / `filter_pairs_by_an_pct` / `encode_genotypes` / `count_all_pairs_via_index` from `compute_vp_counts.py`, and `get_phased_gnomad_ht` from `phase_gnomad.py`.

Trio side (high-quality samples, incl. unreleasable):

| flag | input | output |
|---|---|---|
| `--create-pbt-trio-matrix` | gnomAD VDS, `gnomad_qc.v4` finalized pedigree | `pbt_trio_matrix.{postfix}.mt` |
| `--explode-pbt` | PBT trio matrix | `pbt.{postfix}.mt` |
| `--phase-multi-families` | exploded PBT MT | `pbt_multi_families.{postfix}.mt` |
| `--derive-trio-vps` | PBT MT (probands), variant filter HT | `trio_variant_pairs.{postfix}.ht` |
| `--call-trio-chet` | PBT MT (probands), trio VPs | `trio_phase_counts.{postfix}.ht` (per-pair `n_same_hap`/`n_chet`) |

Comparison side (gnomAD = `release_only`, **all PBT trio members removed**):

| flag | input | output |
|---|---|---|
| `--gnomad-counts-no-pbt` | gnomAD VDS, trio VPs | `gnomad_no_pbt.genotype_counts.{postfix}.ht` |
| `--phase-gnomad-counts` | gnomAD-no-PBT counts | `gnomad_no_pbt.phased.{postfix}.ht` |
| `--export-comparison` | trio truth + gnomAD phase | `trio_comparison.{postfix}.ht` + `.tsv` |

(All trio-phasing outputs are additionally suffixed with `.{trio_set}` — see `--trio-set` below — e.g. `pbt_trio_matrix.{postfix}.pedigree.mt`. The shared variant-filter HT input keeps the plain `{postfix}`.)

Notes:
- The v4 VDS is read split, so no `split_multi_hts` pass is needed (unlike v2). Uses Hail built-ins `phase_trio_matrix_by_transmission` + `explode_trio_matrix` directly. This matches `gnomad_qc.v4` `generate_variant_qc_annotations.run_generate_trio_stats` (densify → adj → `trio_matrix(complete_trios=True)`); that path drops multiallelics whereas we split them.
- **`--trio-set {pedigree, trios}`** (default `pedigree`) selects the pedigree resource: `pedigree()` = all trios, multiple offspring per family (n=15061, the v2 default `fam_path` and what the existing trio stats were run on) vs `trios()` = one random trio per family (n=12731, the v2 `true_trios=True` analog). The value is folded into every trio output path so the two sets never overwrite each other. There is **no separate parent-pair dedup** — family deduplication is expressed solely by choosing `trios`; probands are simply the `s == source_trio.id` columns of the chosen set.
- The finalized pedigree is **exomes-derived** in gnomad_qc v4.
- **No persisted dense trio MT / PBT trio matrix exists in v4 to reuse** — `identify_trios.run_mendel_errors` and `generate_variant_qc_annotations.run_generate_trio_stats` both densify the trio MT transiently and only write aggregate HTs (`ped_mendel_errors`, `trio_stats.ht`). `--create-pbt-trio-matrix` is the first persisted one.
- PBT-sample exclusion from the gnomAD counts is the **direct-filter** approach (build the dense MT `release_only` then `hl.vds.filter_samples(..., keep=False)` on the PBT member set). If that proves too expensive, the alternative is to count PBT∩release and subtract from the production counts — but element-wise subtraction is only exact for the 8 non-AABB cells (AABB needs per-variant `n_callable` bookkeeping; see the AABB gotcha below).
- Per-pop stratification is **not** implemented yet (full-cohort only); the schema/flow leaves room to add it.
- `--phase-multi-families` builds its consensus over `PBT_GT`; the v2 original collected the unphased `GT`, which made its phased-vote sort a no-op.
- `hl.vds.filter_samples` must **not** be called with `remove_dead_alleles=True` on a split VDS — split data has no `LA` field, so it raises `AttributeError: ... 'LA'`.

`v4/resources.py` is the single source of truth for paths and constants:
- `DEFAULT_TMP_DIR = gs://gnomad-tmp-30day` (scratch)
- `VARIANT_COOCCURRENCE_ROOT = gs://gnomad/v4.1/variant_cooccurrence` (production)
- `_get_resource_path` builds `{output_dir}/{data_type}.{resource_name}.{output_postfix}.{ext}`. The `{output_postfix}` argument is **prepended verbatim** — if you pass `--output-postfix sgca_v9_oe_extra_6_12`, you get e.g. `exomes.variant_pairs.sgca_v9_oe_extra_6_12.ht`. The script does **not** auto-prepend a gene name; the gene goes in the postfix.
- `TEST_INTERVALS` maps gene symbol → interval string (used by `--gene <SYMBOL>` to scope all HTs/VDSes).

## Code conventions

Match the existing style in each file rather than imposing new ones.

- Utils functions are pure transformations: HTs in → HTs out. All file I/O lives in CLI scripts.
- Utils functions taking a single Hail Table name the parameter `ht`; use descriptive names only when taking multiple HTs.
- Top-level imports only. Don't add lazy imports unless needed to break a circular import. Right now there are none — `v4/in_trans_oe.py` no longer imports from `v4/create_vp_list.py`; the OE variant-filter-HT builders live in `create_vp_list.py`.
- Don't add docstrings, type annotations, or comments to code you didn't touch. Comments should explain *why*, not what.
- Don't add backwards-compatibility shims. If something is removed, delete it cleanly.

## Hail / pipeline gotchas

These have all bitten the v4 pipeline at some point.

- **`gt_counts_adj` is `array<int64>`**, but `hl.experimental.haplotype_freq_em` requires `array<int32>`. Cast explicitly: `gt_counts.map(hl.int32)`.
- **`p_chet` can be NaN** when EM denominators collapse (e.g., zero double-het pairs). Guard with `hl.if_else(hl.is_nan(p_chet), hl.missing(...), p_chet)`, or equivalently set `o_pair = 0` when `double_carriers == 0`.
- **Pair list is keyed `(locus1, alleles1, locus2, alleles2)` with `v1 ≤ v2`**. The candidate may appear as either `v1` or `v2`. Aggregations need to handle both orientations (`view_a ∪ view_b`).
- **`filter_pair_ht_to_in_trans_oe_pairs`** drops pairs where **both** sides are OE-candidate-only (i.e., source ⊆ `{in_trans_oe_candidate, in_trans_oe_intronic_padding}`). The intent is to avoid candidate × candidate explosion while preserving (OE-candidate × P/LP-partner) pairs.
- **Cache before group-by shuffles**. `pairs = pairs.cache()` before `pairs.group_by(...)` prevents intermittent shuffle errors on the Spark backend.
- **Hail expression source-tracking after `.filter()`**: when you filter a Table and then refer to its fields, bind the filtered Table to a new variable — `ht_sub = ht.filter(...); ht_sub.select(x=ht_sub.foo)`. Re-using `ht.foo` after a filter triggers `Cannot combine expressions from different source objects`.
- **Table indexing requires expressions, not Python literals**. `ht[locus, alleles]` works only if `alleles` is an `hl.literal([...])` or a column expression — passing a plain Python list raises `Cannot index with a scalar expression`.
- **The count steps read the encoded intermediates + size-info from disk; they don't re-encode.** `--encode-genotypes` overwrites `var_idx.ht` + `encoded_gt_sets_by_var_idx.ht` under `genotype_count_intermediates.{postfix}/` on every run (`overwrite=True`), so there's no stale-cache trap within a step. But if you change the pair list or `--min-an-pct`, you must re-run `--encode-genotypes` (and `--build-variant-size-info`) before re-running counts — otherwise `--compute-counts-{light,heavy}` count against the previous run's variant set. Symptom: high-AF candidates with hundreds of pair-list entries but zero gt_counts entries.
- **AABB cell uses `min(n_callable_v1, n_callable_v2)`, not the global cohort size.** The encoder (`_encode_genotype_sets_by_var_idx`) tracks `n_callable` per variant — count of samples with a GT entry at that variant. The AABB formula in `_count_from_sets` (the shared light/heavy/per-pop count kernel) is `min(n_callable_v1, n_callable_v2) - v1_n - v2_n + overlap`. Using the global cohort here would silently bucket missing-GT samples into AABB, inflating the cell at low-coverage positions (3' UTR etc.) and inflating downstream `sum(gt_counts)` (used as `n_pair`) by 3-4× for those variants. If you change the encoder or `_count_from_sets`, keep them in sync.
- **AN_pct floor (`--min-an-pct`) is set at the encode step, applied at consumption, and may only be raised downstream.** The pair list stays the complete raw artifact (no floor baked in). `--encode-genotypes` filters pairs by `--min-an-pct` before densifying and stamps the floor directly onto the encoded HT's globals; each count step re-applies the filter (`filter_pairs_by_an_pct`) and asserts its floor is ≥ the encoded HT's stamped floor (`_assert_min_an_pct_not_lowered`). A lower downstream floor would count pairs whose variants were never densified — the guard turns that into a `ValueError`. Pass `--min-an-pct` *consistently* across the encode/count invocations (like `--output-postfix`). Default `0` drops only uncallable (`an_pct==0`) endpoints; pre-guard artifacts carry no global and are treated as floor `-1` (permissive).
- **Pairs whose variants are missing from the encoded MT are not real comparators** across counting methods. The pair list (`vp_ht`) is built upstream of the dense MT; some pair endpoints can be absent from the MT (different multi-allelic split, filtered out by `--min-an-pct`, etc.) so the `var_idx_ht[locus, alleles].var_idx` lookup returns NULL for those pairs. `compute_counts_light` / `compute_counts_heavy` now drop those pairs upfront and log the count via `_drop_pairs_missing_v_idx` — but the other benchmark methods (`approach_b`, `per_sample`, `v4_original` in `analysis/benchmarks/benchmark_single_gene.py`) each invent a *different* default cell when their inner lookup is missing: approach_b's set algebra propagates NULL into some cells but evaluates others to 0; per_sample collapses to `[NULL, 0, 0, …]`; v4_original treats every sample as filtered-out and dumps them all into AABB (≈ `n_samples` even though no GT data exists). Symptom seen on ANO5: `approach_b vs v4_original` reports `223 raw mism / 223 adj mism` on the full 22,155 pairs, but **0 mismatches** on the 21,932 pairs where both v_idx are defined. **When validating a new counting method, only diff cells on pairs where both variants are in `var_idx_ht`.** light_heavy's dropped-row count is the cleanest signal — if it doesn't match `vp_ht.count()`, the difference is the absent-from-MT pair count for that gene.
- **Complement-form set storage must store the PROPER complement, not just the dominant category.** The decoder identity `|pos ∩ A_pos| = |pos| − |pos ∩ stored|` (used in `_encode_genotype_sets_by_var_idx` when `use_complement=True` for `all_samples`) is only valid when `stored = N \ A_pos` exactly — i.e. `cat_2 ∪ cat_7`, not just `cat_7` (adj-PASS-hom-ref majority). A previous bug stored only `cat_7`, which silently over-counted `|pos ∩ A_pos|` by `|pos ∩ cat_2|` and drove downstream gt_counts cells negative wherever cat_2 is non-trivial (low-AN regions). Fixed in the encoder (`_encode_genotype_sets_by_var_idx`) in `v4/compute_vp_counts.py`; regression tests live in `tests/v4/test_create_vp_matrix_unit.py::TestCountFromSetsRealData` (real chr19 fixtures). Future edits to the encoder **must** preserve the proper-complement invariant; the cheap CI-style sanity check is `sum(gt_counts_adj) ≤ n_samples` on any output HT.
- **`hl.experimental.haplotype_freq_em` will infinite-loop on negative gt_counts.** EM iterations are uncapped, and negative input cells make the likelihood oscillate rather than converge — the symptom is `phase_gnomad.py` stalling for hours on a single partition with no error. The primary defense is the encoder complement invariant above (correct counts converge in ~2 minutes on ~1k partitions); the secondary defense is asserting `sum(gt_counts_adj) ≤ n_samples` before feeding EM; `phase_gnomad.py --em-partitions N` bounds worst-case per-partition wall time if a stall slips through.

## Dataproc submission

Generic invocation pattern (cluster name and SDK / conda paths are in `CLAUDE.local.md`):

```
hailctl dataproc submit <CLUSTER> \
  --pyfiles=/absolute/path/to/gnomad_chets \
  /absolute/path/to/gnomad_chets/v4/<script>.py -- \
  --gene SGCA --output-postfix sgca_v9_oe_extra_6_12 --overwrite
```

`--pyfiles` zips the repo into the executor classpath so `import gnomad_chets.v4.*` works. **Use absolute script paths** — `hailctl` zips up the repo and changes cwd, so relative paths break.

## Testing

The pytest suite lives in a top-level `tests/` module mirroring the package (like `gnomad_methods/tests/`): `tests/conftest.py` holds a session-scoped Hail fixture and tests are under `tests/v4/`. Run with `pytest tests/` (needs the `gnomad_chets` conda env and `gnomad_chets` importable — i.e. `PYTHONPATH=/path/to/PycharmProjects`, since tests do `from gnomad_chets.v4... import ...`).

- `tests/v4/test_utils.py`, `tests/v4/test_create_vp_list.py` — pure-transform unit tests on tiny `hl.Table.parallelize` inputs; run in seconds against local Spark. `test_create_vp_list.py` covers `assemble_sites_ht` (the `--preprocess-sites-ht` step), incl. the drop-on-InbreedingCoeff-agreement gate.
- `tests/v4/test_create_vp_matrix_unit.py` — pytest units but loads real chr19 fixtures at import (needs GCS/Hail).
- `tests/v4/test_create_vp_matrix.py`, `tests/v4/test_in_trans_oe.py` — `__main__` integration / runnable-script harnesses (not pure unit tests).

## Browser integration

The in-trans depletion panel lives in the gnomAD browser repo. The TypeScript types in `browser/src/VariantPage/VariantInTransDepletion.tsx` mirror the JSON shape produced by `run_in_trans_oe.py --output-json`; help text in `browser/help/topics/in-trans-depletion.md` mirrors the formula in `v4/in_trans_oe.py:NULL_MODEL_POISSON_HWE`. Keep them in sync.
