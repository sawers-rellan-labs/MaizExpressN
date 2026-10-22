#!/usr/bin/env Rscript
# M6: pool per-contrast effects across studies.
#   Rscript analysis/pool_meta.R
#
# Input: results/de.parquet (analysis/de_per_study.R) and metadata/contrasts.tsv.
# Only replicated contrasts with pool = TRUE and a finite lfcSE enter pooling (unreplicated and
# display-only contrasts are shown in the app but never pooled; counts are never pooled). Within them,
# genes without a p-value are left out: fold changes carried by one library (single_library) and
# DESeq2 Cook's outliers.
#  1. Within a study, contrasts sharing the grouping (tissue group x class) are combined by
#     inverse variance, so a long time course or many genotypes count as one study estimate.
#  2. Across studies, DerSimonian-Laird random effects (vectorised over genes).
# Groupings, each for filter_set "all" and "B73" (contrasts whose samples are all genotype class B73):
#   results/meta_tissue_class.parquet   tissue_group x contrast_class
#   results/meta_tissue_level.parquet   tissue_group x all N-level classes together (contrast_class =
#                                       "all_levels": deficiency, starvation, resupply, excess,
#                                       field_rate, local; form contrasts are not N levels)
# Private lab data: a catalogue column private = TRUE marks unpublished studies; with
# MAIZEXPRESSN_PUBLIC=1 they are left out of the pooled tables (use the same setting for
# app/prepare_data.R).
# Columns: gene_id, tissue_group, contrast_class, filter_set, k_studies, k_contrasts, pooled_log2fc,
#   se, ci_low, ci_high, pvalue, padj (BH within each grouping and filter set), tau2, i2,
#   n_up, n_down, n_ns (contrasts; padj < 0.05 and |log2FC| >= 1), k_studies_up, k_studies_down.
# tau2 and i2 are NA when only one study contributes.

suppressPackageStartupMessages({ library(data.table); library(arrow) })
PADJ <- 0.05; LFC <- 1
LEVEL_CLASSES <- c("deficiency", "starvation", "resupply", "excess", "field_rate", "local")

cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)
# public build (MAIZEXPRESSN_PUBLIC=1): leave out studies marked private = TRUE in the catalogue
cat_ <- fread("catalogue/maize_N_rnaseq_studies.tsv", colClasses = "character", na.strings = NULL)
PRIVATE <- if ("private" %in% names(cat_)) cat_[toupper(private) == "TRUE", study_id] else character()
if (nzchar(Sys.getenv("MAIZEXPRESSN_PUBLIC")) && length(PRIVATE)) {
  cons <- cons[!(study_id %in% PRIVATE)]
  message("public build: leaving out private studies ", paste(PRIVATE, collapse = ", "))
}
smp <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
b73 <- smp[genotype_class == "B73", sample_id]
cons[, all_b73 := vapply(paste(numerator_samples, denominator_samples, sep = ";"),
                         function(x) all(setdiff(strsplit(x, ";", fixed = TRUE)[[1]], "") %in% b73), logical(1))]

de <- as.data.table(read_parquet("results/de.parquet",
                                 col_select = c("gene_id", "contrast_id", "study_id", "log2fc", "lfc_se", "pvalue", "padj",
                                                "single_library")))
de <- de[cons[, .(contrast_id, contrast_class, tissue_group, pool, all_b73)], on = "contrast_id", nomatch = 0]
de <- de[pool == "TRUE" & is.finite(lfc_se) & lfc_se > 0 & is.finite(log2fc) & !is.na(pvalue) & !(single_library %in% TRUE)]
de[, dir := fifelse(!is.na(padj) & padj < PADJ & log2fc >= LFC, 1L,
               fifelse(!is.na(padj) & padj < PADJ & log2fc <= -LFC, -1L, 0L))]
message(sprintf("%d gene x contrast rows from %d contrasts in %d studies enter pooling",
                nrow(de), uniqueN(de$contrast_id), uniqueN(de$study_id)))

pool_groups <- function(d, keys) {
  by_s <- c(keys, "study_id", "gene_id"); by_g <- c(keys, "gene_id")
  d[, w := 1 / lfc_se^2]
  s <- d[, .(sw = sum(w), swy = sum(w * log2fc), nc = .N, up = sum(dir == 1L), down = sum(dir == -1L),
             ns = sum(dir == 0L)), by = by_s]
  s[, `:=`(y = swy / sw, v = 1 / sw)]
  s[, w := 1 / v]
  a <- s[, .(k_studies = .N, W = sum(w), Sw2 = sum(w^2), Swy = sum(w * y), Swy2 = sum(w * y^2),
             k_contrasts = sum(nc), n_up = sum(up), n_down = sum(down), n_ns = sum(ns),
             k_studies_up = sum(up > 0), k_studies_down = sum(down > 0)), by = by_g]
  a[, Q := pmax(0, Swy2 - Swy^2 / W)]
  a[, tau2 := fifelse(k_studies > 1, pmax(0, (Q - (k_studies - 1)) / (W - Sw2 / W)), 0)]
  s <- s[a[, c(by_g, "tau2"), with = FALSE], on = by_g]
  s[, ws := 1 / (v + tau2)]
  r <- s[, .(sws = sum(ws), swsy = sum(ws * y)), by = by_g]
  a <- a[r, on = by_g]
  a[, `:=`(pooled_log2fc = swsy / sws, se = sqrt(1 / sws))]
  a[, `:=`(ci_low = pooled_log2fc - 1.96 * se, ci_high = pooled_log2fc + 1.96 * se,
           pvalue = 2 * pnorm(-abs(pooled_log2fc / se)),
           i2 = fifelse(k_studies > 1, fifelse(Q > 0, pmax(0, (Q - (k_studies - 1)) / Q), 0), NA_real_),
           tau2 = fifelse(k_studies > 1, tau2, NA_real_))]
  a[, padj := p.adjust(pvalue, "BH"), by = keys]
  a[, c(by_g, "k_studies", "k_contrasts", "pooled_log2fc", "se", "ci_low", "ci_high", "pvalue", "padj",
        "tau2", "i2", "n_up", "n_down", "n_ns", "k_studies_up", "k_studies_down"), with = FALSE]
}

run <- function(d, level) {
  d <- copy(d)
  if (level) d <- d[contrast_class %in% LEVEL_CLASSES][, contrast_class := "all_levels"]
  out <- rbind(pool_groups(copy(d), c("tissue_group", "contrast_class"))[, filter_set := "all"],
               pool_groups(copy(d[all_b73 == TRUE]), c("tissue_group", "contrast_class"))[, filter_set := "B73"])
  setcolorder(out, c("gene_id", "tissue_group", "contrast_class", "filter_set"))
  out[order(tissue_group, contrast_class, filter_set, pvalue)]
}
tc <- run(de, FALSE); write_parquet(tc, "results/meta_tissue_class.parquet")
tl <- run(de, TRUE);  write_parquet(tl, "results/meta_tissue_level.parquet")

summ <- function(x) x[, .(genes = .N, max_k_studies = max(k_studies), genes_k2plus = sum(k_studies >= 2),
                          sig_k2plus = sum(k_studies >= 2 & padj < PADJ & abs(pooled_log2fc) >= LFC)),
                      by = .(filter_set, tissue_group, contrast_class)][order(filter_set, tissue_group, contrast_class)]
fwrite(rbind(summ(tc), summ(tl)), "results/meta_summary.tsv", sep = "\t")
message(sprintf("meta_tissue_class: %d rows; meta_tissue_level: %d rows", nrow(tc), nrow(tl)))
print(summ(tc)[filter_set == "all"])
print(summ(tl)[filter_set == "all"])
