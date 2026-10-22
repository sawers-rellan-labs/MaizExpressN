# Handoff: MaizExpressN

**For:** Claude Code
**From:** Rubén Rellán-Álvarez (GEMMA Lab, NC State), prepared with Claude (Cowork) on 16 Oct 2026
**Goal:** Build an R/Shiny dashboard where lab members look up candidate maize genes and see how their expression changes between nitrogen (N) levels in published RNA-seq studies, with results pooled across studies by tissue, contrast type, growth system and other factors.

Working rules and conventions are in `CLAUDE.md`. This file is the spec.

---

## 1. What is already in this folder

| Path | What it is |
|---|---|
| `catalogue/maize_N_rnaseq_studies.tsv` | **Primary input.** 52 curated studies, one row per paper: design, genotypes, tissue, stage, growth system, N treatments, timing, replication, platform, accessions and accession status. |
| `catalogue/catalogue.json` | The same 52 studies as JSON (used by the catalogue page). |
| `catalogue/README_maize_N_catalogue.md` | How the catalogue was built, column meanings, counts, studies dropped. |
| `catalogue/maize_N_related.tsv` | 36 related items: N microarray and small-RNA studies, meta-analyses that point to datasets, N RNA-seq without an N contrast. Out of scope for v1. |
| `catalogue/maize_N_to_check.tsv` | 46 unverified candidates (papers to read in full, repository IDs without visible titles). Out of scope for v1. |
| `catalogue/build/` | Provenance: search log, the 276 PubMed candidates, raw extraction tables, extraction rules and the merge script. |
| `reference/maize_n_catalogue.html` | Standalone catalogue browser (open in a browser). Use it as the model for the app's Studies page. |

### Catalogue columns that matter most

- `study_id`: unique key (FirstAuthorYear, with a/b where needed).
- `design_category`: `+`-joined tags from LEVEL, TIMECOURSE, FORM, LOCAL, FIELD, NxGENO, NxSTRESS, MUTANT (GRADIENT is in the vocabulary but no core study uses it).
- `tissue_group` (Root, Leaf, Shoot, Ear/kernel, Whole seedling) plus free-text `tissue` and `stage_age`.
- `system_group` (Hydroponic, Solid substrate, Soil/pots, Field, NR) plus free-text `growth_system`.
- `N_treatments`, `other_factors`, `timepoints_duration`, `bio_reps`, `n_libraries`, `platform_layout`: free text taken from the papers. Use them to curate sample sheets by hand; don't parse them blindly.
- `accession`, `accession_status`, `accession_leads`, `accession_note`, `repository`.
- `maizegdb = yes`: MaizeGDB already reprocessed this dataset on B73 v5 ([Hayford et al. 2024](https://doi.org/10.1186/s12864-024-10443-7)); useful for validating our numbers.
- `details_level = abstract only`: 6 studies whose setup details are incomplete.

### Accession status today

| Status | Studies | Meaning |
|---|---|---|
| confirmed | 17 | Number stated in the paper |
| repository match | 3 | Repository record whose title matches the paper (Ravazzolo2020 GSE141860, Cao2024 GSE227472, He2016 PRJNA304223) |
| lead | 3 | ID from searches, unconfirmed (see `accession_leads`) |
| to resolve | 27 | Data statement couldn't be read; number still needs a lookup |
| on request | 1 | Li2023b |
| not deposited | 1 | Fang2024 (data only in the supplement) |

The cloud session that built the catalogue couldn't reach NCBI E-utilities or the ENA API, which is why 27 are unresolved. Hazel or Rubén's Mac should reach them.

## 2. Scope

Two parts:

- **A. Processing pipeline (on Hazel).** Resolve accessions → run tables → curated sample sheets → download → QC → salmon quantification against B73 v5 → gene-level counts and TPM → per-study contrasts with DESeq2 → harmonized Parquet tables.
- **B. Shiny app (on Rubén's Mac first).** Reads the Parquet tables through DuckDB: single-gene lookup, gene-list heatmaps, cross-study pooling and genome-wide ranking.

Out of scope for v1: the related and to-check lists; single-cell data (Li2022) and 21-nt DGE tag data (Pan2015), which need separate handling (keep them as stubs); re-mapping to non-B73 references.

## 3. Environment and constraints

- **HPC:** NC State Hazel, LSF scheduler ([LSF docs](https://hpc.ncsu.edu/Documents/LSF.php), [quick start](https://hpc.ncsu.edu/QuickStart/QuickStart.php)). Group `maize`, Unity ID `rrellan`.
- **Software location:** `/usr/local/usrapps/maize` (quota 100 GB / 250K files). Conda environments use up the file quota quickly, so prefer one pinned environment or Apptainer containers.
- **Raw data location:** FASTQ and salmon outputs must not go in the software directory. Ask Rubén for the group scratch or research-storage path ([HPC data storage](https://hpc.ncsu.edu/Documents/Storage.php)).
- **Volume:** expect roughly 1,000–1,500 runs, several TB if all FASTQ were kept at once. Process study by study, delete FASTQ after quantification, and keep run lists and checksums so anything can be re-downloaded.
- **Network:** check which Hazel nodes reach NCBI and ENA. Download there (login or transfer nodes) and compute on batch nodes.
- **Local:** macOS on Apple silicon, R (Positron is installed), DuckDB.
- **Ask first** before large downloads, deleting anything you didn't create, or changing which studies are in scope.

## 4. Milestones

Do these in order. ⏸ marks a checkpoint where you stop and show Rubén the output.

**Build a vertical slice early.** Run five accession-ready studies end to end (M2–M7) before scaling up. They span root and leaf, hydroponics and field, and level, form and time-course designs:

- Mu2017: B73 leaf, 4 mM vs 40 µM nitrate; GSE107562; in MaizeGDB.
- Ravazzolo2020: B73 root apex, −N vs nitrate vs ammonium; GSE141860; in MaizeGDB.
- Mager2018: B73 roots, −N (and −P); SRP135301.
- Wang2022: B73 primary root, low vs normal nitrate over growth stages; PRJNA858579.
- Chen2015: field ear leaf, two inbreds, low vs normal N; SRP052559.

### M0: Repo scaffold

- Create the layout in §8, `renv` for R, one environment file or container set for CLI tools (sra-tools, salmon, fastp, FastQC, MultiQC, pigz, EDirect).
- Choose the workflow engine. Prefer [nf-core/fetchngs](https://nf-co.re/fetchngs) (downloads SRA, ENA and GEO runs with metadata and writes samplesheets) and [nf-core/rnaseq](https://nf-co.re/rnaseq) in salmon pseudo-alignment mode, if Nextflow and Apptainer run on Hazel under LSF. Otherwise write a small Snakemake pipeline with an LSF profile.
- Add notes to `CLAUDE.md` as decisions are made.

### M1: Resolve accessions ⏸

Output: `metadata/accessions_resolved.tsv` with study_id, accession, accession_type, evidence (link or record title), status (resolved / unavailable / ambiguous), notes.

- For each study, follow NCBI links from the PMID to BioProject, SRA and GEO ([EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/) `elink`), then search BioProject, GEO and SRA by title words and authors.
- Map GEO series to SRA runs with [pysradb](https://github.com/saketkc/pysradb) or [ffq](https://github.com/pachterlab/ffq).
- Verify the leads: Pan2015 (GSE59613, PRJNA255755), Yu2015 (SRX1020537, PRJNA283053), Ravazzolo2021b (GSE173102).
- Check the repository matches: He2016 PRJNA304223 (confirm the rtcs mutant runs are there, not only B73 roots), Cao2024 GSE227472, Ravazzolo2020 GSE141860. Also confirm Ravazzolo2021a's GSE162411 is public; the accepted manuscript cited a reviewer token.
- Jia2023 is in CNCB GSA (CRA012157). Download it from [GSA](https://ngdc.cncb.ac.cn/gsa/), not SRA.
- Huang2025 lists its accessions in supplementary tables; Cheng2021 says GEO. If a supplement is behind a paywall, ask Rubén for the file.

### M2: Run tables and curated sample sheets ⏸

Outputs: `metadata/runs.tsv`, `metadata/samples.tsv`, `metadata/needs_review.tsv` (schemas in §5).

- Pull run metadata (SRA RunInfo plus BioSample attributes; ENA filereport for FASTQ URLs and MD5s).
- Assign every run to the paper's design using the catalogue's free-text fields. Put anything ambiguous in `needs_review.tsv` rather than guessing.
- Drop runs that aren't maize RNA-seq: sorghum runs in Du2020 and Braynen2025; Arabidopsis and tomato in Li2026a; ChIP-seq in Ning2023 (PRJNA837716); small-RNA-seq in Ma2021a (one of PRJNA661965 / PRJNA662035); PacBio Iso-Seq in Wang2020a (keep the Illumina libraries).

### M3: Download, QC and quantification

- FASTQ from ENA where possible (usually faster), otherwise `prefetch` + `fasterq-dump`. Check MD5s.
- QC with fastp/FastQC and MultiQC per study; record read counts and mapping rates.
- Reference: B73 RefGen_v5 (Zm-B73-REFERENCE-NAM-5.0, gene models Zm00001eb.1) from [MaizeGDB downloads](https://download.maizegdb.org/Genomes/NAM_Founders/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.README) (also on [Ensembl Plants](https://plants.ensembl.org/Zea_mays)). Build a decoy-aware [salmon](https://salmon.readthedocs.io/) index (transcripts plus genome as decoys) and summarize to genes with [tximport](https://bioconductor.org/packages/tximport/).
- Salmon with `--libType A` and `--validateMappings` (the default in recent versions).
- 3′ tag-seq (Huang2025, Lexogen QuantSeq with UMIs): add `--noLengthCorrection` and import without length offsets ([QuantSeq with salmon](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6906367/)). Trim the UMI and spacer from read starts, or align with STAR and deduplicate with UMI-tools if exact counts matter.
- Deferred: Pan2015 (21-nt NlaIII tags need a short-k index or an aligner) and Li2022 (single cell). Yu2016 (laser-capture pericycle cells) is standard bulk RNA-seq.
- Non-B73 genotypes are quantified against B73. Record mapping rate by genotype and keep a genotype field so users can restrict to B73.
- Per-study outputs: salmon `quant.sf` per sample, then gene counts and TPM.

### M4: Per-study contrasts and differential expression ⏸

Write `metadata/contrasts.tsv` first (schema in §5) and get it approved.

**Contrast classes**

| Class | Meaning |
|---|---|
| `deficiency` | Steady-state low vs sufficient N |
| `starvation` | −N vs +N after short-term removal |
| `resupply` | N added after starvation vs time-matched −N control |
| `excess` | High vs sufficient N |
| `field_rate` | Low vs high N fertilizer rate in the field |
| `local` | Split-root: local high-nitrate patch vs homogeneous low |
| `form` | Alternative N form vs nitrate (ammonium, urea, mixed) |

**Sign convention:** log2FC = log2(lower N / higher N), so positive means higher under N limitation. This also applies to resupply and local contrasts (lower-N side in the numerator). For `form`, log2FC = log2(alternative form / nitrate). Store a `sign_note` per contrast. The app offers a toggle to flip the sign so users can read results as a response to N addition.

**Differential expression**

- Run [DESeq2](https://bioconductor.org/packages/DESeq2/) per study, with a design matching what the study varies (often a single group factor plus genotype or batch). Use Wald tests per contrast and BH-adjusted p within each contrast.
- Keep unshrunk log2FC and lfcSE for pooling; store [apeglm](https://bioconductor.org/packages/apeglm/)-shrunk log2FC for display.
- Unreplicated studies (Bi2014, Du2020, Singh2022; check Lv2016): store log2FC of normalized counts with SE missing and `replicated = FALSE`. Show them in plots but leave them out of inverse-variance pooling.
- Interaction designs (N × genotype, mutant, P or drought): emit the N effect within each level as its own contrast; store interaction terms in a separate table.
- Time courses: contrast each time point against the time-matched control; if there isn't one, against time 0, and say so in `notes`.

**Validation:** compare DEG counts with numbers reported in the papers. Matching the order of magnitude is enough, since references and pipelines differ. For example, Ravazzolo2020 reported 1,349 nitrate- and 1,264 ammonium-responsive genes vs −N; Cao2024 reported 891 nitrate-induced and 840 nitrate-repressed genes in wild-type roots; Li2026a reported 661 DEGs between split nitrate and control nitrate roots. Also compare with MaizeGDB's processed values for the four datasets flagged `maizegdb = yes`.

### M5: Gene annotation and ID mapping

Output: `annotation/genes.parquet` (schema in §5).

- Source the files from MaizeGDB: gene model correspondence between v3, v4 and v5, gene symbols and descriptions, and GO terms. Start from the [MaizeGDB gene center](https://maizegdb.org/gene_center) and [download.maizegdb.org](https://download.maizegdb.org/). Check the actual files there rather than assuming names.
- Add the Arabidopsis best hit so users can search by Arabidopsis gene names (NRT1.1, NLP7 and so on).
- Search must accept any identifier or symbol, case-insensitive, and show every match when one ID maps to several genes.

### M6: Pooling across studies

Outputs: `results/de.parquet` plus precomputed `results/meta_*.parquet` (schemas in §5).

- **Within a study:** combine contrasts that share tissue group and class (inverse-variance) so each study contributes one estimate per tissue group. This stops long time courses or many-genotype studies from dominating. For the single-gene view, a multilevel model (contrasts nested in studies, [metafor](https://www.metafor-project.org/) `rma.mv`) is a good alternative.
- **Across studies:** random-effects meta-analysis. Use vectorized DerSimonian–Laird for genome-wide tables and REML in metafor for the single-gene view. Report pooled log2FC, 95% CI, p, τ², I², and the number of studies and contrasts.
- **Vote counts:** number of contrasts and studies up, down and not significant (default padj < 0.05 and |log2FC| ≥ 1, both adjustable).
- **Precompute** tables for tissue group × contrast class and tissue group × all level classes, so ranking is fast. Recompute on the fly when the user changes filters.

### M7: Shiny app ⏸ (show the vertical slice first)

**Pages**

1. **Gene lookup.** Search box with autocomplete on v5, v4 and v3 IDs, symbols and Arabidopsis names. A gene card with IDs, description and a link to the MaizeGDB gene page. Then:
   - Forest plot: every contrast's log2FC ± 95% CI, faceted by tissue group and coloured by contrast class, with hollow markers for unreplicated studies. Tooltips show study, genotype, tissue, stage, N levels, time point, padj and the DOI.
   - Pooled panel: one estimate per tissue group (and per contrast class) with CI, I² and the number of studies.
   - Expression panel: TPM by condition for a chosen study, as points plus means, with a log-scale toggle.
   - Time-course panel: TPM or log2FC against time for time-course studies.
   - Results table with CSV download.
2. **Gene list.** Paste or upload IDs or symbols. Show a heatmap of genes × contrasts or genes × pooled estimates per tissue, with optional clustering, plus a summary table and downloads.
3. **Rank genes.** Genome-wide ranking under the current filters (pooled log2FC, p, consistency) and a meta-volcano plot. Clicking a gene opens Gene lookup.
4. **Studies.** The catalogue (filterable, with DOI and accession links), per-study QC (libraries, mapping rate) and the study's contrasts.
5. **Methods.** Pipeline and reference versions, conventions, caveats and citations.

**Global filters** (sidebar, apply to every page): tissue group and detail, contrast class, growth system, stage class (seedling, vegetative, reproductive), genotype (with a B73-only toggle), design tags, time since the N change, replicated studies only, include or exclude individual studies, and the sign toggle.

**Implementation**

- Query Parquet through DuckDB (via [arrow](https://arrow.apache.org/docs/r/) or the [duckdb](https://duckdb.org/docs/current/clients/r) R package) so a gene lookup only touches that gene's rows. Target under 1 s per gene.
- Pick one interactive plotting approach (plotly or ggiraph) and use it throughout. Use a colour-blind-safe palette.
- Run locally first (`shiny::runApp("app")`). Keep the data path configurable so the app can later go to Posit Connect, shinyapps.io or an NC State server.

### M8: Extensibility and tests

- Adding a study, public or the lab's own: new rows in the catalogue, `samples.tsv` and `contrasts.tsv`, then the same pipeline. Mark unpublished lab data `private = TRUE` and keep it out of any public deployment.
- Tests with testthat: sign convention, contrast construction for each design type, ID search, and pooling maths checked against metafor on a few genes.
- Positive controls the app should reproduce: nitrate-assimilation genes (nitrate reductase, nitrite reductase, the NRT2 high-affinity nitrate transporters and their NAR2 partners) should come out strongly nitrate-induced in roots in resupply studies, which means negative log2(lower / higher N).

## 5. Data model

All tables use `study_id` from the catalogue and v5 gene IDs.

**`metadata/runs.tsv`** (one row per sequencing run)
run_accession, study_id, bioproject, biosample, experiment, sample_title, organism, library_strategy, library_source, library_layout, read_length, instrument, n_spots, n_bases, fastq_urls, fastq_md5, biosample_attributes (JSON), include, exclude_reason

**`metadata/samples.tsv`** (one row per biological sample; several runs can belong to one sample)
sample_id, study_id, run_accessions (`;`-joined), genotype, genotype_class (B73 / other inbred / hybrid / mutant / transgenic), background, tissue_group, tissue_detail, po_term (optional Plant Ontology ID), stage, stage_class, age_days, growth_system, N_label_paper (as written in the paper), N_form, N_amount, N_unit (mM or kg N/ha …), N_rank_in_study (1 = lowest), phase (steady / starvation / resupply), time_h (time since the N change; empty for steady state), other_factor, other_factor_level, replicate, pooled_plants, include, curation_note

**`metadata/contrasts.tsv`** (one row per comparison)
contrast_id, study_id, contrast_class, numerator_group, denominator_group, held_constant (genotype, tissue, time, other factor), N_numerator, N_denominator, time_h, replicated, n_numerator, n_denominator, sign_note, notes

**`metadata/needs_review.tsv`**
study_id, item, question, proposed_answer, evidence

**`results/de.parquet`** (one row per gene × contrast)
gene_id, contrast_id, study_id, base_mean, log2fc, lfc_se, log2fc_shrunk, stat, pvalue, padj, single_library, tpm_numerator_mean, tpm_denominator_mean, replicated

**`results/tpm.parquet`** (long format, partitioned by study)
gene_id, sample_id, study_id, tpm, est_count

**`annotation/genes.parquet`**
gene_id, v4_ids, v3_ids, symbol, synonyms, description, chrom, start, end, strand, go_terms, ath_best_hit, ath_symbol

**`results/meta_*.parquet`** (precomputed pooling)
gene_id, grouping keys (tissue_group, contrast_class, filter set), k_studies, k_contrasts, pooled_log2fc, ci_low, ci_high, pvalue, tau2, i2, n_up, n_down, n_ns

## 6. Caveats to keep visible in the app

- Studies differ in lab, genotype, tissue, stage and N units. Never pool counts across studies, only effect sizes.
- N levels are relative within a study (for example 0.04 vs 4 mM nitrate, or 0 vs 240 kg N/ha). Store the actual amounts and units and show them.
- Reads from non-B73 genotypes are mapped to B73, which biases expression for divergent genes.
- Some designs are unreplicated or pooled; 3′ tag and laser-capture libraries behave differently from standard mRNA-seq.
- Some datasets are unavailable (on request, not deposited, or still unresolved). The Studies page should say which.

## 7. Study notes worth knowing up front

- AlonsoNieves2022 and TorresRodriguez2021 share one BioProject (PRJNA669055 / SRP287300). Wild-type libraries are shared; the mutant libraries belong to AlonsoNieves2022. Avoid counting the shared libraries twice when pooling.
- Ravazzolo2021b probably reuses Ravazzolo2020's −N and +nitrate libraries (GSE141860) and adds new +nitrate + cPTIO libraries. Same rule.
- Ma2021b (circRNA) may reuse Ma2021a's libraries (same group, same P178 HN vs LN design).
- Braynen2025: raw reads in PRJNA883413 (which may also hold the sorghum reads); the GEO series hold processed TPM only.
- Cheng2021: field trials in Urbana, IL (2014–2016), many inbreds and hybrids, leaf at R1, 0 vs 200 kg N/ha. Huang2025 (same lab) is a fine B73 time course (0–120 min after N resupply, shoot and root) using QuantSeq.
- He2016 compares the rtcs mutant with wild type; Ge2020 (zmnlp5), Meng2020, Cao2024 (zmnrt1.1b) and Wang2020b (zmga3ox) are also mutant × N designs.
- Yu2016 used laser-capture pericycle cells; Yu2015 used root stele tissue.
- Jia2023 is a deprivation then resupply series (8 time points × leaf and root, inbred Xu178) in CNCB GSA.
- Ning2023: B73 endosperm in the field at 0 vs 240 kg N/ha (RNA-seq in PRJNA602742).

## 8. Proposed repo layout

```
MaizExpressN/
  HANDOFF.md            this spec
  CLAUDE.md             working notes and rules
  catalogue/            curated study list (inputs; don't edit by hand without asking)
  reference/            catalogue browser page
  metadata/             accessions_resolved, runs, samples, contrasts, needs_review
  pipeline/             Nextflow or Snakemake configs and LSF profile for Hazel
  analysis/             R: tximport, DESeq2 contrasts, pooling, QC reports
  annotation/           gene IDs, symbols, GO, Arabidopsis hits
  results/              Parquet outputs synced from Hazel (keep raw data on Hazel)
  app/                  Shiny app
  tests/                testthat tests
```

## 9. Kickoff prompt

> Read HANDOFF.md and CLAUDE.md. Start with M0 and M1: propose the repo scaffold and how you'll resolve accessions, then run the accession resolution and stop with `metadata/accessions_resolved.tsv` for my review. Don't download sequence data yet.

## 10. Key links

- MaizeGDB: [gene center](https://maizegdb.org/gene_center), [downloads](https://download.maizegdb.org/), [B73 v5 README](https://download.maizegdb.org/Genomes/NAM_Founders/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.README); [Ensembl Plants maize](https://plants.ensembl.org/Zea_mays)
- Data access: [EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/), [pysradb](https://github.com/saketkc/pysradb), [ffq](https://github.com/pachterlab/ffq), [nf-core/fetchngs](https://nf-co.re/fetchngs), [CNCB GSA](https://ngdc.cncb.ac.cn/gsa/)
- Processing: [nf-core/rnaseq](https://nf-co.re/rnaseq), [salmon](https://salmon.readthedocs.io/), [tximport](https://bioconductor.org/packages/tximport/), [DESeq2](https://bioconductor.org/packages/DESeq2/), [apeglm](https://bioconductor.org/packages/apeglm/), [metafor](https://www.metafor-project.org/)
- App: [Shiny](https://shiny.posit.co/), [DuckDB R API](https://duckdb.org/docs/current/clients/r), [Apache Arrow for R](https://arrow.apache.org/docs/r/)
- Hazel: [LSF](https://hpc.ncsu.edu/Documents/LSF.php), [quick start](https://hpc.ncsu.edu/QuickStart/QuickStart.php), [data storage](https://hpc.ncsu.edu/Documents/Storage.php)
- Reprocessed maize stress datasets: [Hayford et al. 2024](https://doi.org/10.1186/s12864-024-10443-7)
