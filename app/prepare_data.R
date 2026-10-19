#!/usr/bin/env Rscript
# M7: build the app's data folder from the project outputs.
#   Rscript app/prepare_data.R            (from the project root; writes app/data/)
#   MAIZEXPRESSN_PUBLIC=1 Rscript app/prepare_data.R   public build without private (unpublished) studies
#
# Small tables are written whole and loaded at start-up. Per-gene tables (DE results and TPM)
# are split into ~440 files by gene-ID bucket (the 3 digits after "Zm00001eb"), so a gene lookup
# reads one small file.
#   app/data/genes.parquet, gene_aliases.parquet, go_terms.parquet      (annotation, M5)
#   app/data/contrasts.parquet   one row per contrast in the analysis (tested or display only), with study,
#                                class, tissue, system, stage, genotype class and samples used
#   app/data/samples.parquet     samples in the gene tables, with condition labels and QC
#   app/data/studies.parquet     the studies in the analysis (at least one tested or display-only contrast)
# Only studies used in the analysis reach the app (Rubén, 19 Oct 2026): studies that were dropped, deferred,
# have no data or no contrasts are not listed, and no reasons for dropping a study are shipped. Library QC
# exclusions are not reported per study or per contrast either: library counts and mapping rates describe
# the libraries in use, and Methods gives only the QC rules.
#   app/data/meta_tissue_class.parquet, meta_tissue_level.parquet       (pooling, M6)
#   app/data/de/<bucket>.parquet   results/de.parquet split by gene bucket
#   app/data/tpm/<bucket>.parquet  gene x sample TPM and estimated counts (included, non-excluded samples)

suppressPackageStartupMessages({ library(data.table); library(arrow) })
out <- "app/data"
dir.create(file.path(out, "de"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tpm"), recursive = TRUE, showWarnings = FALSE)
bucket <- function(g) substr(g, 10, 12)
unlink(Sys.glob(file.path(out, c("de", "tpm"), "*.parquet")))   # no stale buckets from an earlier build

# annotation
for (f in c("genes", "gene_aliases", "go_terms"))
  file.copy(file.path("annotation", paste0(f, ".parquet")), file.path(out, paste0(f, ".parquet")), overwrite = TRUE)

# catalogue, samples, QC
cat_ <- fread("catalogue/maize_N_rnaseq_studies.tsv", colClasses = "character", na.strings = NULL)
# public build (MAIZEXPRESSN_PUBLIC=1): leave out studies marked private = TRUE (unpublished lab data);
# run analysis/pool_meta.R with the same setting first so the pooled tables match
PRIVATE <- if ("private" %in% names(cat_)) cat_[toupper(private) == "TRUE", study_id] else character()
PUBLIC <- nzchar(Sys.getenv("MAIZEXPRESSN_PUBLIC"))
drop <- if (PUBLIC) PRIVATE else character()
if (length(drop)) message("public build: leaving out private studies ", paste(drop, collapse = ", "))
smp <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
qc <- fread("metadata/library_qc.tsv", colClasses = "character")
excl <- fread("metadata/library_exclusions.tsv", colClasses = "character")
cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)[!(study_id %in% drop)]
summ <- fread("results/de_summary.tsv", colClasses = "character", na.strings = NULL)

# contrasts
first <- function(x) sub(";.*$", "", x)
cm <- copy(cons)
cm[, first_num := first(numerator_samples)]
cm <- smp[, .(first_num = sample_id, stage_class, genotype_class)][cm, on = "first_num"]
cm <- cat_[, .(study_id, system_group, doi_url, year, first_author, title)][cm, on = "study_id"]
cm <- summ[, .(contrast_id, n_num_used = n_numerator, n_den_used = n_denominator, tested_replicated = replicated,
               n_up, n_down, note)][cm, on = "contrast_id"]
cm[, status := fifelse(!is.na(note) & grepl("no samples", note), "skipped",
                fifelse(is.na(tested_replicated), "not run",
                fifelse(tested_replicated == "TRUE", "tested", "display only")))]
cm <- cm[status %in% c("tested", "display only")]   # skipped / not run contrasts contribute nothing
used <- sort(unique(cm$study_id))
cm[, `:=`(time_h = suppressWarnings(as.numeric(time_h)), pool = pool == "TRUE",
          replicated = tested_replicated == "TRUE",
          system_group = fifelse(system_group %in% c("", "NR"), "Not reported", system_group),
          stage_class = fifelse(stage_class %in% c("", "NR"), "Not reported", stage_class))]
cm <- cm[, .(contrast_id, study_id, contrast_class, tissue_group, tissue_detail, genotype, genotype_class,
             numerator_group, denominator_group, N_numerator, N_denominator, time_h, other_factor_level,
             tags, pool, replicated, status, n_num_used, n_den_used, n_up, n_down,
             system_group, stage_class, doi_url, first_author, year, title, sign_note, notes,
             numerator_samples, denominator_samples)]
write_parquet(cm, file.path(out, "contrasts.parquet"))

# studies in the analysis
st <- copy(cat_)[study_id %in% used]
in_use <- smp[include == "TRUE" & !(sample_id %in% excl$sample_id), sample_id]   # libraries in use only
q <- qc[sample_id %in% in_use, .(libraries = .N, median_mapped_pct = median(as.numeric(percent_mapped))), by = study_id]
st <- q[st, on = "study_id"]
st[, n_contrasts := cm[, .N, by = study_id][st, on = "study_id", x.N]]
write_parquet(st[, .(study_id, first_author, year, title, journal, doi_url, design_category, genotypes, tissue,
                     tissue_group, growth_system, system_group, N_treatments, accession,
                     libraries, median_mapped_pct, n_contrasts)], file.path(out, "studies.parquet"))

# samples and TPM (studies with gene tables; included and not excluded)
studies_q <- intersect(basename(dirname(Sys.glob("results/*/salmon.merged.gene_tpm.tsv"))), used)
keep_s <- smp[include == "TRUE" & !(sample_id %in% excl$sample_id) & study_id %in% studies_q]
keep_s <- qc[, .(sample_id, reads_mapped = as.numeric(reads_mapped), percent_mapped = as.numeric(percent_mapped))][keep_s, on = "sample_id"]
keep_s[, group := sub("_R[0-9]+[a-z]?$", "", sub(paste0("^", study_id, "_"), "", sample_id)), by = study_id]
write_parquet(keep_s[, .(sample_id, study_id, group, genotype, genotype_class, tissue_group, tissue_detail,
                         N_label_paper, N_form, N_amount, N_unit, N_rank_in_study, phase, time_h,
                         other_factor_level, replicate, reads_mapped, percent_mapped)],
              file.path(out, "samples.parquet"))

tpm_long <- rbindlist(lapply(studies_q, function(s) {
  t <- fread(file.path("results", s, "salmon.merged.gene_tpm.tsv"))
  k <- fread(file.path("results", s, "salmon.merged.gene_counts.tsv"))
  cols <- intersect(setdiff(names(t), c("gene_id", "gene_name")), keep_s$sample_id)
  if (!length(cols)) return(NULL)
  lt <- melt(t[, c("gene_id", cols), with = FALSE], id.vars = "gene_id", variable.name = "sample_id", value.name = "tpm")
  lk <- melt(k[, c("gene_id", cols), with = FALSE], id.vars = "gene_id", variable.name = "sample_id", value.name = "est_count")
  d <- lt[lk, on = c("gene_id", "sample_id")]
  d[, `:=`(sample_id = as.character(sample_id), study_id = s)]
  d[tpm > 0 | est_count > 0]
}))
tpm_long[, b := bucket(gene_id)]
for (bb in unique(tpm_long$b))
  write_parquet(tpm_long[b == bb, .(gene_id, sample_id, study_id, tpm = round(tpm, 3), est_count = round(est_count, 1))],
                file.path(out, "tpm", paste0(bb, ".parquet")))
message(sprintf("tpm: %d rows in %d buckets", nrow(tpm_long), uniqueN(tpm_long$b)))
rm(tpm_long); invisible(gc())

# DE results by bucket
de <- as.data.table(read_parquet("results/de.parquet"))[contrast_id %in% cm$contrast_id]
de[, b := bucket(gene_id)]
for (bb in unique(de$b))
  write_parquet(de[b == bb][, b := NULL], file.path(out, "de", paste0(bb, ".parquet")))
message(sprintf("de: %d rows in %d buckets", nrow(de), uniqueN(de$b)))

for (f in c("meta_tissue_class", "meta_tissue_level"))
  file.copy(file.path("results", paste0(f, ".parquet")), file.path(out, paste0(f, ".parquet")), overwrite = TRUE)
writeLines(format(Sys.time(), "%d %b %Y %H:%M"), file.path(out, "BUILT"))
message("app data written to ", out)
