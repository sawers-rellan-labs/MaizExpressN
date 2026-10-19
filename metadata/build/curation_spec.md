# Sample-sheet curation spec (M2)

Goal: for each study, map every run in `metadata/runs.tsv` (rows with that `study_id`) to a biological sample with its experimental design. One study's output goes in `metadata/build/samples/<study_id>.tsv`; open questions go in `metadata/build/samples/<study_id>.review.tsv`. Write only those two files for your studies.

Accuracy over speed. Never guess a design detail: if the run metadata and the paper don't settle it, use `NR` (not reported) or flag it in the review file. Don't download sequence data.

## Sources, in order of trust
1. Run and sample metadata in `metadata/runs.tsv`: `sample_title`, and in `biosample_attributes` (JSON) the submitter `alias`, `title` and attributes (treatment, tissue, cultivar, time, replicate...). These say which run is which sample.
2. The paper. For PMC papers, the web page has full text including back matter: `curl -sL -A "Mozilla/5.0" https://pmc.ncbi.nlm.nih.gov/articles/<PMCID>/` (strip tags). Use Methods, figure legends and tables. Supplements are linked from that page. GEO sample pages (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM...&targ=self&form=text`) often add characteristics.
3. The catalogue row (`catalogue/maize_N_rnaseq_studies.tsv`): a summary written by an earlier pass. It's helpful but may be wrong; the paper wins.
4. `metadata/accessions_resolved.tsv` notes for that study (duplicates, shared runs, exclusions).

Be polite to servers: at most one request per second, and don't hammer one site.

## Output: `<study_id>.tsv` (tab-separated, UTF-8, header row, no tabs or newlines in fields)

One row per biological sample (library). Several runs of the same library (lanes, re-sequencing) belong to one sample. Include every run of the study that has `include = TRUE` in runs.tsv, exactly once. Runs with `include = FALSE` are left out.

| Column | Content |
|---|---|
| sample_id | `<study_id>_<label>`, e.g. `Mu2017_LN_R1`. Unique, no spaces. |
| study_id | as in the catalogue |
| run_accessions | `;`-joined run accessions for this sample |
| genotype | as named in the paper, e.g. `B73`, `Mo17`, `zmnrt1.1b-1`, `B73xMo17` |
| genotype_class | one of `B73`, `other inbred`, `hybrid`, `mutant`, `transgenic` (mutant = induced or natural mutant line incl. CRISPR knockouts and insertion lines; transgenic = overexpression / RNAi lines) |
| background | genetic background for mutants/transgenics (e.g. `W22`, `B73`, `ND101`); the line itself for inbreds/hybrids |
| tissue_group | one of `Root`, `Leaf`, `Shoot`, `Ear/kernel`, `Whole seedling` (Shoot = whole shoot or stem without leaf-only sampling; husk leaf = Leaf; endosperm, cob, florets, ears = Ear/kernel) |
| tissue_detail | as specific as the paper allows, e.g. `primary root apex (0-1 cm)`, `ear leaf`, `leaf elongation zone`, `phloem-pole pericycle (LCM), crown root` |
| po_term | Plant Ontology ID if obvious (e.g. `PO:0009005` root, `PO:0025034` leaf, `PO:0009025` vascular leaf, `PO:0020127` primary root, `PO:0025142` leaf lamina, `PO:0009089` endosperm); else empty |
| stage | as written, e.g. `V3`, `14 d after sowing`, `R1`, `15 DAP` |
| stage_class | `seedling` (up to ~V4 or <= 21 d old, hydroponic seedlings), `vegetative` (later vegetative), `reproductive` (VT/R stages, ears, kernels), or `NR` if the stage is unknown |
| age_days | plant age at sampling in days if stated or computable, else empty |
| growth_system | free text, e.g. `hydroponic, growth chamber`, `field, Urbana IL 2015`, `soil pots, greenhouse` |
| N_label_paper | label used in the paper/metadata, e.g. `LN`, `HN`, `-N`, `+NO3`, `N0`, `N240` |
| N_form | `nitrate`, `ammonium`, `ammonium nitrate`, `urea`, `mixed` (state the mix in curation_note), `none` (zero-N control), or `NR` |
| N_amount | number only, the amount of N itself (e.g. `0.04`, `4`, `240`); `0` for no N |
| N_unit | `mM N` for solutions: salt molarity x N atoms (KNO3 x1, Ca(NO3)2 x2, NH4NO3 x2, (NH4)2SO4 x2, urea x2), so 2.5 mM Ca(NO3)2 = 5 mM N; `kg N/ha` for field rates (convert product rates: 450 kg urea/ha = 207 kg N/ha); other units only if neither applies. Put the salt and the original figure in curation_note. |
| N_rank_in_study | integer, 1 = lowest N supply within the study (ties share a rank: e.g. nitrate vs ammonium at the same amount) |
| phase | what happened to this sample's own N supply at time 0: `steady` (unchanged), `starvation` (N removed or reduced, incl. a KCl mock after N-containing pre-growth or -N continued after starvation), `resupply` (N added or increased) |
| time_h | hours since the experiment's N change, for every sample of a time-resolved design (controls included, so they can be time-matched; 0 for time-zero samples); empty for steady-state designs |
| other_factor | non-N factor varied in the design, e.g. `P supply`, `water`, `inoculation`, `inhibitor`, `hormone treatment`, `leaf zone`; empty if none |
| other_factor_level | e.g. `low P`, `water stress`, `H. seropedicae ZAE94`, `cPTIO`, `control` |
| replicate | biological replicate number within its group (1, 2, 3...) |
| pooled_plants | plants pooled per sample if stated (a single number; ranges go in curation_note), else empty |
| include | `TRUE`, or `FALSE` with the reason in curation_note (e.g. a sample outside any N contrast) |
| curation_note | short provenance: how the run was mapped (e.g. `alias LN-R1`), plus anything unusual. Mention batch differences (different instrument/lane/submission between groups). |

Notes on the design fields:
- Every sample must be placeable in an N contrast later. Time-course controls (e.g. -N at the same time point, or time 0) are samples like any other.
- For field N rates use the fertilizer rate (kg N/ha) and the form (`urea` etc.) if stated.
- For mutant x N designs, both genotypes get their own samples; the N levels stay in the N fields.
- Shared runs: if a run's `also_used_by` lists another study, curate it under its owning `study_id` only.
- If a study has samples that are not part of any N contrast (e.g. a hormone-only treatment at one N level), still list them with `include = TRUE` and explain in curation_note; contrasts will decide later. Use `include = FALSE` only for samples that clearly can't be used (failed library, wrong organism or tissue mislabelled with no fix).

## Output: `<study_id>.review.tsv`

Columns: `study_id`, `item`, `question`, `proposed_answer`, `evidence`. One row per open issue: an ambiguous run-to-sample mapping, a mismatch between paper and repository (e.g. replicate counts, genotype labels), missing N amounts, unexpected extra runs, batch confounding, anything a human must decide. Write the header even if there are no issues.

## Finish
Check your file: every included run of the study appears exactly once; columns are exactly as listed, in that order; controlled vocabularies are respected. Then reply with a 3-5 line summary per study (samples, design, open issues).
