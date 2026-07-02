# Variant Co-occurrence Pipeline: v2 to v4 Code Changes

This document summarizes the changes in the variant co-occurrence (chet) pipeline from gnomAD v2 to v4, how they improve efficiency, and how they accommodate the much larger variant and variant-pair scales in v4.

---

## 1. Pipeline Structure

### v2 Pipeline (3 main steps)

| Step | Flag | Output | Notes |
|------|------|--------|------|
| 1 | `--create_vp_list` | Variant pair list Table | Optional `--vp_list_by_chrom` to compute by chromosome and union. |
| 2 | `--create_full_vp` | **Full variant-pair MatrixTable** | Rows = variant pairs, cols = samples; **caused memory challenges** (see code TODO). |
| 3 | `--create_vp_summary` | Summary Table with genotype counts | Reads full VP MT, aggregates by population. |

Additional steps: `--create_vp_ann` (annotations), PBT-specific steps (`--create_pbt_summary`, `--create_pbt_trio_ht`).

### v4 Pipeline (6 discrete steps)

| Step | Flag | Output | Notes |
|------|------|--------|------|
| 1 | `--create-variant-filter-ht` | Variant filter Table | QC, consequence, AF; no genotype data. |
| 2 | `--filter-vmt` | Filtered VDS/MatrixTable | Full sample set, filtered variants only. |
| 3 | `--create-variant-pair-list-ht` | Variant pair list Table | Unique ordered pairs per gene/sample. |
| 4 | `--create-dense-filtered-mt` | Dense filtered MatrixTable | **Variants only** (rows = variants in any pair), not variant pairs. |
| 5 | `--create-variant-pair-genotype-ht` | Variant pair genotype Table | Per-pair genotype info in Table form. |
| 6 | `--create-variant-pair-genotype-counts-ht` | Genotype counts Table | Raw and adj count arrays per pair. |

v4 **never** builds a variant-pair × sample MatrixTable. Genotype information is carried in Tables keyed by variant pair, with localized encodings and aggregations.

---

## 2. Key Code and Architectural Changes

### 2.1 Variant pair list creation

**v2** (`create_variant_pair_ht`):

- Uses **entries table**: `mt.select_cols().select_rows(*row_groups).entries()`.
- Groups by `*row_groups, *mt.col_key` (e.g. gene + sample), collects variant structs per cell, then generates pairs with `flatmap`/`map` over indices and explodes.
- Single aggregation over the full MT entries; pair explosion happens in one pass.

**v4** (`create_variant_pair_ht`):

- **Gene–sample grouping first**: converts to entries, explodes on `gene_id`, then `group_by("gene_id", "s").aggregate(variants=collect_as_set(...))`.
- **Checkpoint** after grouping to avoid huge lineage and re-computation.
- Filters to `variants.length() >= 2`, then generates ordered pairs via `flatmap`/`map`, explodes, keys by pair, `distinct()`.
- Adds **unique index** (`vp_ht_idx`) and keys by it for downstream joins.
- Uses `hl._set_flags(use_new_shuffle="1")` for shuffle stability at scale.

**Efficiency:** Checkpointing after group-by and indexing by pair + `vp_ht_idx` keeps the pair list manageable and makes later steps (dense MT, genotype annotation) partition-friendly and join-efficient.

### 2.2 From “full VP MatrixTable” to “dense variant MT + Tables”

**v2** (`create_full_vp`):

- Builds a **variant-pair × sample MatrixTable**:
  - Indexes variant list by `locus2, alleles2` with `all_matches=True`, filters rows with at least one match, explodes so each row is one (v1, v2) pair.
  - Joins with full MT to get GT/adj for both v1 and v2 per sample.
  - Multiple checkpoints and a large repartition (e.g. 10000) to cope with size.
- **Problem:** The resulting MT has one row per variant pair and one column per sample. Size is **O(variant_pairs × samples)**. With millions of pairs and hundreds of thousands of samples, this does not scale and led to the documented “memory challenges.”

**v4** (no full VP MT):

- **Dense MT** (`create_dense_filtered_mt`): rows = **variants** (union of v1 and v2 in the pair list), columns = samples. Size is **O(variants × samples)**, not O(pairs × samples).
- **Variant pair genotype Table** (`create_variant_pair_genotype_ht`): each row is a variant pair; genotype info is stored in a **localized** form (per-sample arrays of encoded GTs) in a Table, not as a giant matrix.
- **Variant pair genotype counts Table** (`create_variant_pair_genotype_counts_ht`): aggregates those per-sample arrays into raw/adj count arrays per pair; no VP MT ever materialized.

**Efficiency:** Avoiding the full VP MatrixTable removes the main memory and shuffle bottleneck. Data flow is: variant-level MT → encoded Table → pair-level Table with localized genotypes → aggregated counts.

### 2.3 Genotype encoding and localization

**v2:**

- Keeps full `GT` (and adj) in the VP MT entries; `create_vp_summary` uses `get_counts_agg_expr(mt)` over entries to compute 9-bin counts per population.

**v4** (`_encode_and_localize_genotypes`):

- **Encodes** genotypes as small integers: missing = hom_ref (stored as missing for space), 0 = missing data, 1 = het, 2 = hom_var (and similarly for adj).
- **Localizes** the dense MT to a Table: one row per variant, with an array of (sample_index, raw_gt, adj_gt).
- **Filters** to “called” variants only (drops hom_ref/missing for that variant), so each row’s array is much shorter than full sample count.

**Efficiency:** Smaller types and dropping hom_ref from the per-variant arrays greatly reduce storage and transfer. Counts are then computed by aggregation over these arrays in the pair Table, not over a full VP matrix.

### 2.4 Annotating pairs with genotype info (v4-only)

**v4** uses a **split-and-union** strategy so it never holds both variants’ full genotype vectors in one huge structure:

- `_prepare_variant_pair_index`: from the pair Table (keyed by `vp_ht_idx`), build two Tables keyed by (locus, alleles) for v1 and v2, each with `vp_ht_idx` and a “which variant” tag; **union** and `collect_by_key` so each (locus, alleles) maps to a collection of (vp_ht_idx, vp). Repartition and checkpoint.
- `_annotate_variant_pairs_with_genotypes`: join this unioned Table with the **localized genotype Table** (keyed by locus, alleles), then group by `vp_ht_idx` and collect (vp, gt_info) per pair.

**Efficiency:** Joins are by variant key; each variant’s genotype vector is read once and attached to all pairs that use it, then grouped by pair. This avoids the v2 pattern of exploding to a full VP MT and keeps memory and shuffle size bounded by variant count and pair count separately, not by their product with sample count.

### 2.5 Count aggregation

**v2** (`create_vp_summary`, `get_counts_agg_expr`):

- Operates on the full VP MT entries: for each (pair, sample) cell, a 9-element count vector is computed from GT1/GT2/missing/adj; then `hl.agg.group_by(pop, hl.agg.array_agg(sum, ...))` over entries.

**v4** (`_calculate_genotype_counts`, `_convert_gt_info_to_counts`):

- Input is the **Table** of variant pairs with per-sample genotype arrays (already encoded).
- For each pair, per-sample (v1_gt, v2_gt) is extracted, filtered/remapped (hom_ref/missing handling), then `hl.agg.counter` over the list of [v1_gt, v2_gt].
- Counter dict is turned into the same 9-bin count array via `_convert_gt_info_to_counts` (including `n_samples_filtered_out` for hom_ref/hom_ref). Done for both raw and adj.

**Efficiency:** No VP MatrixTable scan. All work is per-row aggregation on the pair Table, which is partition-friendly and scales with number of pairs and length of per-pair arrays (bounded by “called” samples per pair).

---

## 3. Data and Scale

### 3.1 v2 (gnomAD v2 exomes, reference)

Numbers below are from Guo et al., *Nat Genet.* 2024 (inferring compound heterozygosity from gnomAD exomes) and the gnomAD browser variant co-occurrence release (August 2021).

- **Samples:** 125,748 exomes after quality control.
- **Genes:** 19,877 genes (paper); 19,685 genes in the browser release.
- **Unique variant pairs (total):** 5,320,037,963 pairs meeting the criteria (both variants global AF ≤5%, coding/flanking intronic/UTR, same gene).
- **Variant pairs carried by same individual (≥1 sample):** 11,786,014 pairs.
- **Variant pairs in browser release:** 20,921,100 pairs across 19,685 genes.
- **Singleton pairs (no phase prediction):** 105,322 pairs (both variants singleton, same individual).

The v2 pipeline was limited by the need to materialize and shuffle a variant-pair × sample MatrixTable; at 11.8M+ pairs × 125K samples the full VP MT was already challenging.

### 3.2 Why v2 design does not scale to v4 sizes

- **Full VP MatrixTable:**  
  - v2: 1.18×10⁷ pairs (carried by same individual) × 1.26×10⁵ samples ≈ 1.5×10¹² cells.  
  - v4 (projected): co-occurring pair count will grow with more variants, genes, and samples—not as a simple multiple of v2. A full VP matrix would be **~10⁸–10⁹ pairs** × **730,947 samples** ≈ **10¹³–10¹⁴ cells**.  
  Materializing and shuffling this is not feasible.

- **v4 design:**  
  - Dense MT: variants × samples (**variants** × **730,947**), not pairs × samples.  
  - Pair Table: **~10⁸–10⁹ rows** with compact per-pair genotype arrays and count arrays; no (pairs × samples) matrix.

**References for v2 numbers:** Guo et al., Inferring compound heterozygosity from large-scale exome sequencing data, *Nat Genet.* 56:152–161 (2024); gnomAD browser, “Variant Co-Occurrence (Phasing) Information in gnomAD” (August 2021).

**References for v4 numbers:** gnomAD v4.0 (November 2023) and v4.1: **730,947** individuals with exome sequencing (release samples); total 807,162 individuals (exome + genome). Variant and pair counts are projected ranges from scaling with sample size and variant discovery, not yet from a full v4 exome run.

---

## 4. Size and Runtime Projections (v2 → v4)

### 4.1 Approximate scale factors (v2 → v4 projected)

These are rough projections, not exact multiples; actual v4 numbers depend on variant/sample growth and filters.

- **Samples:** v2 125,748 → v4 **730,947** exomes (**~5.8×**).
- **Variants (filtered):** v2 full exome ~hundreds of thousands across 19,877 genes → v4 full exome **~1–2 million (projected)**.
- **Variant pairs:** v2 11.8M (carried by same individual; 5.32B total possible) → v4 full exome co-occurring pairs **~100M–500M+ (projected)**; scaling is not a simple multiple of v2 because pair count depends on variant count, sample count, and co-occurrence.

### 4.2 What grows in v4 (projected)

| Quantity | v2 | v4 (projected) |
|----------|----|-----------------|
| Samples | 125,748 | **730,947** |
| Genes | 19,877 | full exome |
| Filtered variants | (hundreds of thousands) | **~1–2 million** |
| Variant pairs (co-occur in ≥1 sample) | 11,786,014 | **~100M–500M+** |
| Variant pairs (total possible, same criteria) | 5,320,037,963 | — |
| VP matrix (if built) | pairs × 125,748 | N/A (not built) |
| Dense variant MT | variants × 125,748 | variants × **730,947** |
| Pair Table + counts | summary only | **~100M–500M+ rows** |

### 4.3 Runtime (projected)

- Step 2 (variant pair list) and step 3 (dense MT) will dominate runtime and scale with variant count and pair count (and number of partitions).
- Step 4 (genotype + counts) will scale with pair count and per-pair sample list length; design is intended to keep this as Table aggregation without ever building the full VP matrix.

---

## 5. Summary Table

| Aspect | v2 | v4 |
|--------|----|----|
| **Core bottleneck** | Full variant-pair × sample MatrixTable | Avoided; no VP MT |
| **Genotype storage** | Full GT in VP MT entries | Encoded, localized, filtered to “called” only |
| **Pair genotype annotation** | Index + explode in MT | Split by variant, join by (locus, alleles), then group by pair index |
| **Counts** | Entry aggregation on VP MT | Table aggregation on pair Table with small arrays |
| **Scalability** | Limited by O(pairs × samples) | O(variants × samples) for dense MT; O(pairs) for pair Tables |
| **Pipeline steps** | 3 main (list → full VP → summary) | 6 (filter → VDS → list → dense MT → genotype HT → counts HT) |
| **Checkpointing** | In create_full_vp and summary | After group-by in pair list; after encode/localize; after prepare_variant_pair_index |

These changes allow the v4 pipeline to handle the larger variant sets and variant pair sets (and sample sizes) of gnomAD v4 while staying within feasible memory and shuffle sizes by never materializing the variant-pair × sample matrix and by encoding and localizing genotype data in Tables.

---

## 6. Summary of Required Updates

*(This section aligns with grant-language summaries of the required updates.)*

We had on the order of **~5.2 billion** potential variant pairs in the gnomAD v2 dataset and anticipate scale to increase by **~6×** when using the gnomAD v4 dataset (e.g., ~6× more samples and many more variant pairs). The codebase, written in Hail, needed to be substantially updated to accommodate this increase.

The v2 design used a **variant-pair × sample matrix**: one row per variant pair and one column per sample. At v4 scale (~6× as many samples and variant pairs), that matrix would be impractically large to build and shuffle. We replaced the need for this large matrix with a **multi-stage pipeline** that: (1) first filters variants; (2) creates smaller **variant-level** genotype representations—e.g., compressed arrays of only the samples that carry the alternative allele—stored in Tables; (3) annotates these onto a table of unique variant pairs via indexed joins and group-by aggregations. This allows us to replace the global aggregation over the full matrix with **many, smaller aggregations** over per-pair data, keeping memory sizes reasonable. Importantly, this approach also **scales with the data**, making it better able to handle larger future datasets.

Implementing this redesign demanded substantial developer effort. It required designing and implementing six discrete pipeline steps with clear interfaces and checkpointing between stages; refactoring genotype handling (encoding hom_ref as missing, filtering to called variants only) to reduce storage and I/O; implementing a split-and-union strategy for annotating pairs with genotype data that avoids exploding to a full matrix; and adapting to Hail and Spark behavior at scale (e.g., shuffle flags, partition counts, and join patterns). The work also included writing tests and documentation, tuning for different Hail versions and cluster types (e.g., standard vs high-memory), and validating correctness of the new count logic against the previous pipeline’s semantics. Delivering a production-ready v4 pipeline that can handle **~100M+ variant pairs** and **730K samples** therefore required dedicated software engineering and computational genomics expertise, not just re-running existing code on larger hardware.
