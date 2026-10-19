CONTEXT
You are helping a maize geneticist (NC State, GEMMA lab) build a catalogue of published RNA-seq experiments in maize grown under different nitrogen (N) conditions. The catalogue will drive bulk download of raw reads next, so DATA ACCESSIONS and an accurate experimental design matter most. Accuracy over speed: never guess.

STEP 1 - TRIAGE
Get metadata/abstracts with mcp__PubMed__get_article_metadata (about 10-12 PMIDs per call).
INCLUDE a paper only if BOTH:
 (a) it generated NEW transcriptome sequencing data: bulk mRNA-seq, total/strand-specific RNA-seq (incl. lncRNA/circRNA/alternative-splicing studies), 3' tag-seq/QuantSeq, cell-type/LCM RNA-seq, single-cell/single-nucleus or spatial RNA-seq, or Iso-Seq/full-length combined with expression; AND
 (b) maize (Zea mays; teosinte/other Zea is OK but flag it) was sampled under >=2 N supply conditions: different N levels (deficiency/low vs sufficient/high; N gradient; N starvation), N starvation -> resupply or nitrate-induction time course, N forms (nitrate vs ammonium vs urea/mixed), local/heterogeneous N (split-root, nitrate patch), or N crossed with other factors (genotypes, hybrids, mutants/transgenics, drought, P, density, microbes) PROVIDED N supply itself was varied.
Put in RELATED (not the main table): microarray-only; small RNA/miRNA/degradome-only; ribosome-profiling-only; papers that only REANALYSE public RNA-seq data (record which datasets they used, if stated).
EXCLUDE: qPCR-only, proteomics/metabolomics-only, reviews, no N treatment (e.g., biostimulant/biochar/bacterial inoculation without N levels), non-Zea species, microbe/soil-only transcriptomes.
If an abstract is ambiguous (e.g., "transcriptome" without saying RNA-seq vs microarray), check the full text before deciding.

STEP 2 - EXTRACT (each included study)
Get the PMCID (mcp__PubMed__convert_article_ids) and read the full text with mcp__PubMed__get_full_text_article, ONE article per call: Methods, figure legends, Data Availability / accession numbers. If there is no PMC full text, use the abstract plus WebSearch (e.g., the paper title + "BioProject" / "PRJNA" / "GSE" / "SRA" / "CRA") to find the accession.
WebFetch is rate-limited in this session: use it at most ~3 times in total, only for a publisher page when essential. If any fetch is refused or rate-limited, do not retry it and do not try to get that content any other way (no curl, no Python, no mirrors or caches): mark the field 'not found' and move on. Never use Bash for web access.
NEVER invent or infer accessions, numbers, units or genotypes. Use 'NR' when a detail is not reported and 'not found' for an accession you could not locate. Distinguish 'not stated' (no data statement) from 'on request'.

STEP 3 - OUTPUT
Write a TSV (tab-separated, UTF-8, header row, no tabs or newlines inside fields; use '; ' inside fields) with the Write tool to the path you were given. Do not modify any other files. Columns, in this order:
study_id	pmid	pmcid	doi_url	first_author	year	journal	title	design_category	genotypes	n_genotypes	tissue	stage_age	growth_system	N_treatments	other_factors	timepoints_duration	bio_reps	n_libraries	platform_layout	accession	repository	other_omics	key_finding	source_of_details	notes
Field guidance:
- study_id: FirstAuthorSurnameYear (e.g., Bi2014); add a/b if needed.
- doi_url: https://doi.org/<doi>
- design_category: one or more of LEVEL (two N levels), GRADIENT (>=3 levels), TIMECOURSE (starvation/resupply or induction kinetics), FORM (N forms), LOCAL (split-root/patch), FIELD (field-grown with N rates), NxGENO (contrasting genotypes/hybrids/panels), NxSTRESS (with drought/P/other), MUTANT (mutant/transgenic x N); join with '+'.
- genotypes: names (e.g., 'B73; Mo17; B73xMo17'); note contrasting NUE / low-N tolerance if claimed; flag teosinte.
- tissue: e.g., 'seedling root; shoot', 'ear leaf', 'developing kernels 15 DAP', 'lateral-root pericycle (LCM)'.
- stage_age: e.g., 'V3', '14 d after sowing', 'silking (R1)', '7 d after treatment'.
- growth_system: hydroponic / sand / vermiculite / soil pots / field (+ location, years); growth chamber vs greenhouse.
- N_treatments: concentrations or rates with units and N form, e.g., 'LN 0.04 mM vs HN 4 mM NO3-' or '0 vs 240 kg N/ha (urea)'.
- other_factors: e.g., 'drought', 'P levels', 'mutant zmnlp5', 'planting density', or 'none'.
- timepoints_duration: treatment duration and sampling times.
- bio_reps: biological replicates per condition (as stated). n_libraries: total RNA-seq libraries if stated or computable.
- platform_layout: e.g., 'Illumina HiSeq 2500, PE150', 'BGISEQ-500 SE50', 'Iso-Seq + Illumina'.
- accession: ALL accessions (BioProject PRJNA/PRJEB/PRJDB/PRJCA, SRA SRP/SRR, GEO GSE, CNCB GSA CRA, ArrayExpress E-MTAB), '; '-separated; or 'on request' / 'not stated' / 'not found'.
- repository: NCBI SRA / NCBI GEO / CNCB GSA / ENA / DDBJ / Figshare / other.
- other_omics: metabolomics, proteomics, small RNA, ionomics, physiology... or 'none'.
- key_finding: <=20 words.
- source_of_details: 'PMC full text', 'abstract + web search', 'publisher page', etc.
- notes: caveats (e.g., 'pooled samples, no replicates', 'same dataset as <study_id>').

FINAL REPLY (compact - the TSV file is the main deliverable):
 1. Path of the TSV and number of included studies.
 2. RELATED: one line each - 'PMID | FirstAuthorYear | data type (microarray/sRNA/degradome/reanalysis) | accession(s) if stated | <=10-word description'.
 3. EXCLUDED: comma-separated 'PMID (<=4-word reason)'.
 4. Problems (tools that failed, papers you could not read).
If writing the file fails, paste the TSV in the reply instead.
