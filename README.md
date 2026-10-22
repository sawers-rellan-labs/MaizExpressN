# MaizExpressN
[DOI Zenodo](https://zenodo.org)(https://doi.org/10.5281/zenodo.23577994) 
Look up a maize gene and see how its expression responds to nitrogen across published RNA-seq
studies, study by study and pooled across studies by tissue, contrast class, growth system,
stage and genotype. GEMMA Lab, NC State (Rubén Rellán-Álvarez).

The full specification is `HANDOFF.md`; working notes and every decision taken are in `CLAUDE.md`.

Large files (Parquet tables, per-study gene tables, MultiQC reports, annotation sources) are stored with
[Git LFS](https://git-lfs.com). Install it before cloning (`brew install git-lfs`, then `git lfs install`), or clone
with `GIT_LFS_SKIP_SMUDGE=1` to get only the code and pointer files.

## Run the app

```bash
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

Run from the project root (the renv library activates there). The app reads `app/data/`, built by
`app/prepare_data.R` and `app/prepare_partners.R`; set `MAIZEXPRESSN_DATA` to use a data folder elsewhere.

## Deploy to Posit Connect Cloud

Live app: <https://01a155ac-51db-275e-103a-4c139985dd6d.share.connect.posit.cloud/> (public link; managed at
<https://connect.posit.cloud/gemmalab/content/01a155ac-51db-275e-103a-4c139985dd6d>, account `gemmalab`).

From an R session started in the project root (so the renv library is active):

```r
rsconnect::connectCloudUser()                               # once: authorise this computer in the browser
rsconnect::deployApp(appDir = "app", appName = "maizexpressn")  # upload app/ (about 430 MB) and publish
```

Re-run `deployApp()` after `Rscript app/prepare_data.R` and `Rscript app/prepare_partners.R` to update the live app. For publishing from
GitHub instead, refresh `app/manifest.json` with `rsconnect::writeManifest(appDir = "app")` and commit it
with the `app/` folder. Do a public build (`MAIZEXPRESSN_PUBLIC=1`) if any study is marked private.

## Layout

| Folder | What is in it |
|---|---|
| `catalogue/` | The study catalogue (`maize_N_rnaseq_studies.tsv`), the source of truth for which studies are in scope |
| `metadata/` | Curated tables: `accessions_resolved.tsv`, `runs.tsv`, `samples.tsv`, `contrasts.tsv`, `needs_review.tsv`, `library_qc.tsv`, `library_exclusions.tsv`; builders and per-study sample sheets in `metadata/build/` (see its README) |
| `pipeline/` | Download and quantification on NC State's Hazel cluster (Slurm, nf-core/rnaseq); see `pipeline/README.md` |
| `results/` | Per-study gene tables copied from Hazel (`results/<study>/`), `de.parquet` (M4), `meta_*.parquet` (M6), validation notes |
| `annotation/` | Gene annotation (M5): sources in `annotation/source/`, `genes.parquet`, `gene_aliases.parquet`, `go_terms.parquet` |
| `analysis/` | `de_per_study.R` (DESeq2 per study) and `pool_meta.R` (pooling across studies) |
| `app/` | The Shiny app (`app.R`, `R/`, `www/`, `methods.md`) `prepare_data.R` and `prepare_partners.R`. Pages: Gene lookup, Gene list, Co-expression, Rank genes, Studies, Methods |
| `tests/` | testthat suite: `Rscript tests/testthat.R` |

## Rebuild everything after a change

From the project root, in this order (each step reads the outputs of the one before):

```bash
python3 metadata/build/build_runs.py
python3 metadata/build/merge_samples.py --write
python3 metadata/build/build_contrasts.py
python3 metadata/build/validate_contrasts.py
python3 metadata/build/library_exclusions.py
Rscript analysis/de_per_study.R
Rscript analysis/pool_meta.R
python3 annotation/build_genes.py && Rscript annotation/build_parquet.R
Rscript app/prepare_data.R
Rscript app/prepare_partners.R   # co-expression partner lists, ~15 min
Rscript tests/testthat.R
```

## Add a study (public or the lab's own)

1. **Catalogue.** Add a row to `catalogue/maize_N_rnaseq_studies.tsv` with a new `study_id`
   (FirstAuthorYear). For unpublished lab data add a column `private` and set it to `TRUE`;
   `MAIZEXPRESSN_PUBLIC=1 Rscript analysis/pool_meta.R` and `MAIZEXPRESSN_PUBLIC=1 Rscript app/prepare_data.R`
   (then `Rscript app/prepare_partners.R`, which reads the app data)
   then leave it out of a public deployment.
2. **Accessions and runs.** Record the accession in `metadata/accessions_resolved.tsv`, then run
   `build_runs.py`. Lab data that is not in a public repository needs its runs added by hand
   (run ID, file paths, MD5s) in the same columns; any run-level exclusion goes in a `RULES`
   entry in `build_runs.py` with a reason.
3. **Samples.** Write `metadata/build/samples/<study_id>.tsv` (and a `.review.tsv`) following
   `metadata/build/curation_spec.md`: genotype, tissue, N amount in `mM N` or `kg N/ha`, `N_rank_in_study`
   (1 = lowest), phase and time. `merge_samples.py --write` validates it.
4. **Contrasts.** `build_contrasts.py` writes the study's contrasts from the sample groups using the
   sign convention below; check them with `validate_contrasts.py`.
5. **Quantify on Hazel.** `python3 pipeline/make_manifests.py <study_id>`, sync `pipeline/` to Hazel, then
   `sbatch --job-name=dl_<study_id> download_study.sbatch <study_id>` and
   `sbatch --job-name=<study_id> run_study.sh <study_id>` (details in `pipeline/README.md`). Copy
   `salmon.merged.gene_counts.tsv`, `gene_tpm.tsv`, `gene_lengths.tsv` and the MultiQC report to `results/<study_id>/`.
6. **QC.** Run `pipeline/library_qc.py` on Hazel into `metadata/library_qc.tsv`, then `library_exclusions.py`.
   Drop the study if it fails the reliability rules in `CLAUDE.md`.
7. **Analyse and publish.** `de_per_study.R`, `pool_meta.R`, `app/prepare_data.R`, `app/prepare_partners.R`, then the tests.

## Conventions

- Gene key: B73 RefGen_v5 gene IDs (`Zm00001eb…`); v4, v3, symbols and Arabidopsis names are search aliases.
- Effect sign: log2FC = log2(lower N / higher N), so positive means higher under N limitation; N-form
  contrasts are log2(alternative form / nitrate). The app's sign toggle flips it.
- Counts are never pooled across studies; only log2FC with its SE.
- Library QC is strict (`metadata/build/library_exclusions.py`); unreplicated contrasts are shown but never pooled.

## License

The code is released under the [MIT License](LICENSE) (© 2026 Rubén Rellán-Álvarez and the GEMMA Lab).
The expression results are derived from published studies whose raw reads are in public repositories;
please cite the original studies (listed with their DOIs in `catalogue/maize_N_rnaseq_studies.tsv`) when using them.
