# Maize nitrogen RNA-seq catalogue — README

Built 16 Oct 2026. Purpose: list published RNA-seq experiments in maize (Zea) that contrast nitrogen supply, with the setup and data accessions needed to download and re-quantify the reads.

## Files
- `maize_N_rnaseq_studies.tsv` — 52 core studies (one row per paper). Main input for the download step.
- `maize_N_related.tsv` — 36 related items: N microarray and small-RNA studies, reprocessed/meta-analysis papers (pointers to datasets), and RNA-seq of N physiology without an N-supply contrast.
- `maize_N_to_check.tsv` — 46 items not yet included: papers that need full-text reading and repository IDs that came up in searches without a visible title.

## Inclusion
New RNA-seq (bulk, lncRNA/total, 3' tag/QuantSeq, LCM, single-cell, Iso-Seq + RNA-seq) of Zea tissue under two or more experimentally controlled N supply conditions: N levels, starvation/resupply time courses, N forms, split-root/local nitrate, field N rates, or N crossed with genotype, mutant or another stress.

Dropped from the main list (16 Oct 2026, not core): Griffiths2021 and Sugimura2022, Fudjoe2025, Li2026b, Lopinski2025 (N confounded with other nutrients or inputs, not set as a treatment, or host genes not analysed); Zhang2024, McLoughlin2018, Xie2025 (RNA-seq or N design not confirmed from the abstract).

## Key columns
- `design_category`: LEVEL, GRADIENT, TIMECOURSE, FORM, LOCAL, FIELD, NxGENO, NxSTRESS, MUTANT (joined with +).
- `accession_status`: confirmed (number in the paper), repository match (repository title matches the paper), lead (ID from searches, unconfirmed; see `accession_leads`), to resolve (data statement not readable here), on request, not deposited.
- `maizegdb`: yes = already reprocessed on B73 v5 in MaizeGDB/qTeller (Hayford et al. 2024).
- `details_level`: full text vs abstract only (6 studies meet the criteria from the abstract but lack full setup details).

## Counts
20 studies with an accession (17 confirmed, 3 repository matches); 3 leads; 27 to resolve; 1 on request; 1 not deposited.

## Known gaps
- The PubMed Central full-text service drops hyperlinked accession numbers and some data-availability sections, so many deposited datasets show "to resolve".
- NCBI E-utilities and ENA/Europe PMC APIs were not reachable from the cloud workspace; web search budget ran out before every lead was checked.

## Next steps (planned)
1. Resolve the 27 "to resolve" accessions (NCBI PubMed-to-BioProject/SRA/GEO links and title searches, run where NCBI is reachable, e.g. the HPC).
2. Expand each BioProject/GSE into run tables (SRR lists, sample-to-condition mapping).
3. Download and quantify uniformly against one reference (B73 v5) to build a queryable expression matrix.
