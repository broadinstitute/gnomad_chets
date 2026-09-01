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
| compute_vp_counts | `--merge-subset-counts` | per-subset counts HTs | `variant_pairs.genotype_counts.{postfix}.ht` |

`--encode-genotypes` densifies only the pair-list variants out of the VDS in-memory (transient scratch checkpoint, no persisted dense MT) before encoding — there is no separate dense-MT step. It also keeps read-backed phase (`PGT` + `PID`) and stores a per-variant `phased_het` sidecar (`sample_idx → {pid, gt0}`); the full-cohort count steps then emit `n_phased_cis` / `n_phased_trans` (default on via `--emit-phase-counts`), the physical-phase split of the double-het `AaBb` cell — the same signal the gnomAD MNV pipeline derives. Per-pop counts (`--stratify-by-pop`) don't carry phase yet, and phase is not yet consumed by `phase_gnomad.py` / `in_trans_oe.py`. `--encode-genotypes` also applies the **v4 high-AB het → hom-alt correction** (see gotcha below), so it now joins the release freq HT (`get_freq`, per-variant `af`) and the meta HT (per-sample `fixed_homalt_model`).

The full-cohort count steps also accept **`--emit-no-pbt-counts`** (with `--trio-set {pedigree,trios}`, matching `trio_phasing.py`): they additionally emit `gt_counts_raw_no_pbt` / `gt_counts_adj_no_pbt` — the 9-cell counts with the PBT (trio) members subtracted out (`release \ PBT`). Mechanism: the count step **restricts the same encode to `PBT∩cohort` once** (`restrict_encoded_to_samples` — re-indexes the sets to a dense PBT space, checkpointed) and counts it to `counts_pbt.ht`; `--combine-counts` then subtracts element-wise (`_subtract_pbt_counts`). This is **exact for all 9 cells, AABB included** (disjoint-cohort additivity — see the PBT-exclusion gotcha); both arrays come from the same encode so they share the identical adj + high-AB classification. **Why restrict-once, not per-pair**: the earlier count-time version embedded an `hl.literal(PBT∩release index set)` in the per-pair count expression, which Hail replicated into every light/heavy task and blew past `spark.rpc.message.maxSize` (835 MB/task). The keep-set literal now lives only in the one checkpointed re-index transform (`_restrict_encoded_to_indices`, shared with `restrict_encoded_to_pops`). `trio_phasing.py --gnomad-counts-no-pbt --gnomad-counts-path <that HT>` then reuses these columns directly — no separate gnomAD-no-PBT densify/count pass. Full-cohort only (incompatible with `--stratify-by-pop`).

**By-sample subsetting (`--vds-subset` / `--merge-subset-counts`)** splits the expensive densify+encode+count across a disjoint by-sample partition of the cohort and sums the results back. Julia pre-split the raw v4.0 exomes VDS (`gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds` — what `get_gnomad_v4_vds` reads) into 10 column subsets: `non_ukb` + `ukb.<group>` for the 9 v4 gen-anc groups (`resources.get_count_subsets` / `get_count_subset_vds_path`). `compute_vp_counts --vds-subset <name> --encode-genotypes …` reads that subset VDS by path (`_read_count_subset_vds`: same release/HQ meta semi-join + interval + chr19-multiallelic-drop as the full loader, handed to `densify_encode_input_mt` via its new `vds=` param), and writes **subset-qualified** count-group outputs (`genotype_count_intermediates.{postfix}.{subset}/`, `variant_pairs.genotype_counts.{postfix}.{subset}.ht`) via `get_variant_pair_resources(subset=…)` — the pair-list/filter **inputs stay the full-dataset artifacts** (subsetting is counts-only; the pair list is always built on everything). `--merge-subset-counts [--merge-subsets a,b,…|all]` then element-wise-sums the per-subset count HTs (`merge_subset_counts`) into the plain full-cohort path — **exact for all 9 cells incl. AABB** by disjoint-cohort additivity, GIVEN two things: (1) every subset carries every pair-list variant as a row (the v4.0 splits preserve all variant rows — see `analysis/ukb_vds_split_runs.md`; the `ukb.vds` intermediate has since been deleted), and (2) **the release/high-quality restriction is applied AFTER densify, not on the sparse VDS.** That second point is the subtle part: `hl.vds.filter_samples` on the sparse VDS PRUNES any variant row with zero entries among the kept samples, so restricting the subset to release *before* densify drops the row (and hence the hom-ref baseline) of every pair-list variant monomorphic within that subset — which for a small subset is most of them. `hl.vds.to_dense_mt` instead PRESERVES those rows and fills them hom-ref from the reference blocks, so `_read_count_subset_vds` returns ALL cohort samples and `densify_encode_input_mt(restrict_samples_ht=…)` restricts to release as a dense-MT COLUMN filter (rows kept). Diagnosed on CAPN3 (2026-08-31): the pre-fix version undercounted AABB (merged ≈ ¼ of full on rare pairs, `gt_counts_adj` matched only 10,060/70,492) because `_read_count_subset_vds` did `filter_samples(release)` up front — staged row counts showed raw `ukb.mid`=11,699 rows over CAPN3, `to_dense_mt`=11,699 (preserved), but `filter_samples(release)`=290. **Validated after the fix (2026-09-01):** all 10 subsets → merge → diff vs a full-cohort CAPN3 run is **bit-identical on every one of the 70,510 pairs the full run emits** (raw + adj + phase, 0 mismatches), no subset dropping any pair. The merge additionally keeps 744 pairs the **full run itself drops** — `get_gnomad_v4_vds`'s pre-densify `filter_samples(release)` prunes pair-list variants that have zero release carriers in the VDS (freq-based pair list vs VDS-release cohort mismatch), which the densify-then-restrict subset path retains as inert all-hom-ref pairs; so the merge is exact on the full run's universe and a strict superset of it. Exomes only; incompatible with `--stratify-by-pop` / `--pops` / `--emit-no-pbt-counts` (run those on the merged full-cohort counts).

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
| `--gnomad-counts-no-pbt` | production counts HT w/ `_no_pbt` cols (`--gnomad-counts-path`) **or** gnomAD VDS + trio VPs (from scratch) | `gnomad_no_pbt.genotype_counts.{postfix}.ht` |
| `--phase-gnomad-counts` | gnomAD-no-PBT counts | `gnomad_no_pbt.phased.{postfix}.ht` |
| `--export-comparison` | trio truth + gnomAD phase | `trio_comparison.{postfix}.ht` + `.tsv` |

(All trio-phasing outputs are additionally suffixed with `.{trio_set}` — see `--trio-set` below — e.g. `pbt_trio_matrix.{postfix}.pedigree.mt`. The shared variant-filter HT input keeps the plain `{postfix}`.)

Notes:
- The v4 VDS is read split, so no `split_multi_hts` pass is needed (unlike v2). Uses Hail built-ins `phase_trio_matrix_by_transmission` + `explode_trio_matrix` directly. This matches `gnomad_qc.v4` `generate_variant_qc_annotations.run_generate_trio_stats` (densify → adj → `trio_matrix(complete_trios=True)`); that path drops multiallelics whereas we split them.
- **`--trio-set {pedigree, trios}`** (default `pedigree`) selects the pedigree resource: `pedigree()` = all trios, multiple offspring per family (n=15061, the v2 default `fam_path` and what the existing trio stats were run on) vs `trios()` = one random trio per family (n=12731, the v2 `true_trios=True` analog). The value is folded into every trio output path so the two sets never overwrite each other. There is **no separate parent-pair dedup** — family deduplication is expressed solely by choosing `trios`; probands are simply the `s == source_trio.id` columns of the chosen set.
- The finalized pedigree is **exomes-derived** in gnomad_qc v4.
- **No persisted dense trio MT / PBT trio matrix exists in v4 to reuse** — `identify_trios.run_mendel_errors` and `generate_variant_qc_annotations.run_generate_trio_stats` both densify the trio MT transiently and only write aggregate HTs (`ped_mendel_errors`, `trio_stats.ht`). `--create-pbt-trio-matrix` is the first persisted one.
- PBT-sample exclusion from the gnomAD counts uses **subtraction from the same release encode**: `release \ PBT = release − (PBT∩release)`, element-wise on the 9-cell arrays. This is **exact for all 9 cells, AABB included** — every cell is an additive per-sample count over the disjoint cohorts (`release = (release\PBT) ⊎ (PBT∩release)`), and the AABB `n_samples` baseline (`AABB = n − |D'_v1 ∪ D'_v2|`) also splits additively across them. The one requirement is that both arrays use the **same genotype classification** (same adj + high-AB correction), so a sample lands in the same cell in both — guaranteed by computing both from the *same* release encode. `--emit-no-pbt-counts` therefore restricts that encode to the (small) `PBT∩release` cohort **once** (`restrict_encoded_to_samples`, re-indexed + checkpointed), counts it to `counts_pbt.ht`, and subtracts at `--combine-counts` (`_subtract_pbt_counts`); it must NOT be sourced from a separately-built MT with a different adj/correction. **Do not** put the PBT-index-set as an `hl.literal` in the per-pair count expression — Hail replicates it into every count task and it blows past `spark.rpc.message.maxSize` (that footgun is why the restriction is a one-shot re-index, not a per-pair intersect). The production run (`compute_vp_counts --emit-no-pbt-counts`) does this; `trio_phasing.py --gnomad-counts-no-pbt --gnomad-counts-path <that HT>` then reuses the emitted `_no_pbt` columns directly (restrict to trio VPs, promote to `gt_counts_{raw,adj}`) — no densify. Its from-scratch fallback (omit `--gnomad-counts-path`) instead densifies release-minus-PBT via `hl.vds.filter_samples(keep=False)`.
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
- **No-entry (uncallable) samples must stay out of the AABB (hom-ref/hom-ref) cell.** `_count_from_sets` (the shared light/heavy/per-pop kernel) derives AABB from the global `n_samples` as `n_samples − |D'_v1 ∪ D'_v2|` (plus F terms for the raw cell), where `D'_v = all_samples ∪ raw_hr_adj_missing = cats 1-6`. Correctness hinges on the encoder putting **cat 1 (no-entry) samples in `all_samples`**: because they're in `D'_v` they fall outside `H_v = N \ D'_v` (the adj-PASS-0/0 hom-ref set), so a sample with no GT at either variant is never counted as hom-ref. If `_encode_genotype_sets_by_var_idx` ever stopped including no-entry samples in `all_samples`, they'd leak into AABB and inflate it (and `sum(gt_counts)`, used as `n_pair`) at low-coverage positions (3' UTR etc.). Cheap CI check: `sum(gt_counts_adj) ≤ n_samples`.
- **AN_pct floor (`--min-an-pct`) is set at the encode step, applied at consumption, and may only be raised downstream.** The pair list stays the complete raw artifact (no floor baked in). `--encode-genotypes` filters pairs by `--min-an-pct` before densifying and stamps the floor directly onto the encoded HT's globals; each count step re-applies the filter (`filter_pairs_by_an_pct`) and asserts its floor is ≥ the encoded HT's stamped floor (`_assert_min_an_pct_not_lowered`). A lower downstream floor would count pairs whose variants were never densified — the guard turns that into a `ValueError`. Pass `--min-an-pct` *consistently* across the encode/count invocations (like `--output-postfix`). Default `0` drops only uncallable (`an_pct==0`) endpoints; pre-guard artifacts carry no global and are treated as floor `-1` (permissive).
- **Pairs whose variants are missing from the encoded MT are not real comparators** across counting methods. The pair list (`vp_ht`) is built upstream of the dense MT; some pair endpoints can be absent from the MT (different multi-allelic split, filtered out by `--min-an-pct`, etc.) so the `var_idx_ht[locus, alleles].var_idx` lookup returns NULL for those pairs. `compute_counts_light` / `compute_counts_heavy` now drop those pairs upfront and log the count via `_drop_pairs_missing_v_idx` — but the other benchmark methods (`approach_b`, `per_sample`, `v4_original` in `analysis/benchmarks/benchmark_single_gene.py`) each invent a *different* default cell when their inner lookup is missing: approach_b's set algebra propagates NULL into some cells but evaluates others to 0; per_sample collapses to `[NULL, 0, 0, …]`; v4_original treats every sample as filtered-out and dumps them all into AABB (≈ `n_samples` even though no GT data exists). Symptom seen on ANO5: `approach_b vs v4_original` reports `223 raw mism / 223 adj mism` on the full 22,155 pairs, but **0 mismatches** on the 21,932 pairs where both v_idx are defined. **When validating a new counting method, only diff cells on pairs where both variants are in `var_idx_ht`.** light_heavy's dropped-row count is the cleanest signal — if it doesn't match `vp_ht.count()`, the difference is the absent-from-MT pair count for that gene.
- **Encoded sets are stored positive form (no complement encoding).** `_encode_genotype_sets_by_var_idx` stores `all_samples` (cats 1,3-6) and `raw_hr_adj_missing` (cat 2) directly; cat 7 (adj-PASS-0/0 majority) is never stored — `_count_from_sets` reconstructs it as `N − (cats 1-6)` from `n_samples` via plain set-intersection algebra. The old "store the smaller of set/complement + `*_is_complement` flag" optimisation was **removed** (commit on `jg/v4-pipeline`): it halved worst-case row size but its proper-complement invariant (`stored = cat_2 ∪ cat_7`, not just `cat_7`) was a subtle, once-buggy footgun. Tradeoff: at very-low-AN variants the positive `all_samples` (large cat-1 no-entry) can approach `n_samples`, so run the count steps with a **non-zero `--min-an-pct`** to bound row/shuffle size. Invariant to keep: `sum(gt_counts_adj) ≤ n_samples` on any output HT (regression tests: `tests/v4/test_create_vp_matrix_unit.py::TestCountFromSetsRealData`, real chr19 fixtures normalised to positive form in the harness).
- **High-AB het → hom-alt correction (v4 consistency).** GATK <4.1.4.1 mis-called some true hom-alts as high-allele-balance hets; gnomAD v4 release freq corrects these ([`gnomad_qc` `generate_freq.py`](https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v4/annotations/generate_freq.py)). `_encode_genotype_sets_by_var_idx` mirrors it: a call is reclassified het→hom-var (`1→2`) when `GT.is_het_ref() & adj & (AD[1]/DP > HIGH_AB_CUTOFF=0.9) & ~_het_non_ref & ~fixed_homalt_model & (af > HIGH_AB_AF_THRESHOLD=0.01)`. **Applied to BOTH raw and adj**, matching the release: the finalized `freq` is `ab_adjusted_freq`, whose `correct_for_high_ab_hets` adds the (adj-determined) high-AB hom-alt count to *every* stratum including `freq[1]` = raw (the raw-group aggregate of the adj-gated `high_ab_het`). So a corrected call sits in `raw_hv` AND `adj_hv`. The raw reclassification is **gated on adj-pass** (`high_ab_homalt & adj_pass_expr`) — a non-adj high-AB het is not in gnomAD's correction set, so it stays het in raw. (`gnomad_qc`'s docstring says "raw is not adjusted" but the code + released data adjust it; verify against released `freq[1].AC` if in doubt.) Requires three inputs wired onto the dense MT by the encode step: per-variant `af` (release `get_freq().freq[0].AF`, same source as `create_vp_list`), per-sample `fixed_homalt_model` (`meta().project_meta`), and per-entry `_het_non_ref` (`LGT.is_het_non_ref()` captured **pre-split** in `_split_variant_data_keeping_phase`). The correction auto-skips when any of those fields is absent (e.g. the trio PBT MT), so `use_precomputed_adj` callers are unaffected. The `fix_freq_an.py` v4.0→v4.1 AN fix is **not** replicated — it was specific to per-subset (UKB/non-UKB) freq computation; we densify the full release cohort in one pass and use the corrected `all_sites_an` for `an_pct`.
- **Sex-ploidy adjustment (v4 consistency, sex chromosomes only).** `densify_encode_input_mt` applies `adjusted_sex_ploidy_expr(locus, GT, sex_karyotype)` right after densify — matching `gnomad_qc` `generate_freq` `densify_and_prep_vds_for_freq`, so adj + the high-AB correction are computed on the sex-adjusted GT (order: `_het_non_ref` pre-split → densify → sex ploidy → adj → high-AB). It's a **strict no-op on autosomes** (its first branch is `in_autosome`), so per-chromosome autosomal runs (the default plan) are unaffected; only chrX/chrY jobs change. On non-PAR X/Y: an XY call becomes haploid (`Call(gt[0])`) and an XY het is dropped to missing; an XX call on Y becomes missing. The existing classifier handles all three with no change: a **hemizygous alt (haploid `Call(1)`) → `is_hom_var` → the `aa`/`bb` hom cell** (the chosen co-occurrence modeling; note this makes its per-variant AC 2 not 1, so revisit if AC-matching on chrX is ever needed), a haploid ref → hom-ref, and a dropped het (missing GT, defined entry) → `raw_gt=0` → in `all_samples`/`n_with_data` but no genotype cell, so it's uncallable and excluded from AABB. Needs the per-sample `sex_karyotype` col (`meta().sex_imputation`). Sex chromosomes may still warrant bespoke handling (gnomAD's MNV pipeline punts on hemizygous modeling) — decide per chrX/chrY run.
- **`hl.experimental.haplotype_freq_em` will infinite-loop on negative gt_counts.** EM iterations are uncapped, and negative input cells make the likelihood oscillate rather than converge — the symptom is `phase_gnomad.py` stalling for hours on a single partition with no error. Positive-form set algebra can't produce the old complement-bug negatives, but assert `sum(gt_counts_adj) ≤ n_samples` before feeding EM as a cheap guard; `phase_gnomad.py --em-partitions N` bounds worst-case per-partition wall time if a stall slips through.
- **Keeping read-backed phase (`PGT`/`PID`) requires our own split — `gnomad_qc`'s loader drops it.** `_split_and_filter_variant_data_for_loading` only remaps `{GT, AD, PL}` → local names, so `LPGT` is neither selected pre-split nor produced post-split (and `get_gnomad_v4_vds` also rejects `filter_variant_ht` on unsplit reads). `hl.experimental.sparse_split_multi` itself *does* downcode `LPGT`→`PGT` and pass `PID` through — so `--encode-genotypes` reads `split=False` (samples + intervals only) and calls `_split_variant_data_keeping_phase`, which adds `PGT` to the remap set (matching gnomad_mnv), applies the variant filter itself (locus pre-filter + post-split semi-join), and reconstructs the VDS before `to_dense_mt`. If you switch the encode densify back to `get_gnomad_v4_vds(split=True, ...)`, phase is silently lost and `phased_het` goes empty.

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
