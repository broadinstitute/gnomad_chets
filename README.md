# gnomad_chets
Phase of compound heterozygotes in gnomAD

This repository serves as a home for the pipeline used to infer the phase of rare variants in the gnomAD v2 exomes, as reported in our corresponding manuscript (see https://www.biorxiv.org/content/10.1101/2023.03.19.533370v2), and is coded in Hail 0.2.

The main components of the pipeline can be found in “phasing.py”. Briefly, to infer variant phase, we generate haplotype frequency estimates from genotype counts by applying the expectation-maximization (EM) algorithm (see “get_em_expressions” function) and calculate the probability of two variants being in trans (compound heterozygous, “p_chet”), using the haplotype frequency estimates in a simple equation (“get_em_expr” function).

The remaining scripts in the repository serve to compute the phase of rare variant pairs specifically in the gnomAD and Center for Mendelian Genetics rare disease datasets, and to generate the gnomAD variant co-occurrence look-up tool (see https://gnomad.broadinstitute.org/variant-cooccurrence) and variant co-occurrence counts by gene resource (see https://gnomad.broadinstitute.org/news/2023-03-variant-co-occurrence-counts-by-gene-in-gnomad/). These scripts cannot be run outside of the gnomAD team, as they require access to the individual level data, and are provided for reference only.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10034663.svg)](https://doi.org/10.5281/zenodo.10034663)
