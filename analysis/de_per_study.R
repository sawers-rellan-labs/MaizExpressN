#!/usr/bin/env Rscript
# M4: per-study differential expression for every contrast in metadata/contrasts.tsv.
#
#   Rscript analysis/de_per_study.R              # all studies with gene tables in results/
#   Rscript analysis/de_per_study.R Mu2017 ...   # selected studies (writes results/de_<study>.parquet only)
#
# One DESeq2 model per study: design ~ 0 + grp, where grp is the sample group (the same
# group key build_contrasts.py uses) and the model holds every included sample of the study
# plus any samples its contrasts borrow from another study. Libraries in
# metadata/library_exclusions.tsv are left out. Each contrast is a Wald test of the mean of
# its numerator groups against the mean of its denominator groups, so
# log2FC = log2(lower N / higher N) (form contrasts: alternative form / nitrate), as in
# contrasts.tsv. Counts come from salmon via tximport with average-transcript-length offsets;
# 3' tag studies (Huang2025) use counts without length offsets.
# Unshrunk log2FC and lfcSE are kept for pooling; ashr-shrunk log2FC is for display (apeglm
# cannot shrink arbitrary contrasts). Unreplicated contrasts get the log2 ratio of mean
# normalised counts (+0.5) and no SE, p or padj.
#
# Outputs: results/de.parquet (gene x contrast; schema in HANDOFF section 5) and
# results/de_summary.tsv (one row per contrast: samples used, genes tested, DEG counts).

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(arrow)
})
stopifnot(requireNamespace("ashr", quietly = TRUE))  # used by lfcShrink(type = "ashr")
options(warn = 1)

NO_LENGTH <- c("Huang2025")       # 3' tag-seq: no length offsets
PADJ_DEG <- 0.05; LFC_DEG <- 1      # DEG call for the summary table only
MIN_COUNT <- 10                    # a gene is tested if >= MIN_COUNT reads in >= MIN_SAMPLES samples
GROUP_KEY <- c("study_id", "genotype", "tissue_group", "tissue_detail", "N_label_paper",
               "N_form", "N_amount", "N_unit", "N_rank_in_study", "phase", "time_h",
               "other_factor", "other_factor_level")

samples <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)
excl <- fread("metadata/library_exclusions.tsv", colClasses = "character")
samples[, group := do.call(paste, c(.SD, sep = "|")), .SDcols = GROUP_KEY]
use <- samples[include == "TRUE" & !(sample_id %in% excl$sample_id)]

split_ids <- function(x) setdiff(strsplit(x, ";", fixed = TRUE)[[1]], "")

# gene x sample matrices from each owning study's nf-core tables, cached
FILES <- c(counts = "salmon.merged.gene_counts.tsv", lengths = "salmon.merged.gene_lengths.tsv",
           tpm = "salmon.merged.gene_tpm.tsv")
cache <- new.env()
read_mat <- function(owner, kind) {
  key <- paste(owner, kind)
  if (is.null(cache[[key]])) {
    d <- fread(file.path("results", owner, FILES[[kind]]))
    m <- as.matrix(d[, -(1:2)]); rownames(m) <- d$gene_id
    cache[[key]] <- m
  }
  cache[[key]]
}
get_mat <- function(sids, kind) {
  owners <- samples[match(sids, sample_id), study_id]
  parts <- lapply(unique(owners), function(o) read_mat(o, kind)[, sids[owners == o], drop = FALSE])
  genes <- rownames(parts[[1]])
  m <- do.call(cbind, lapply(parts, function(p) p[genes, , drop = FALSE]))
  m[, sids, drop = FALSE]
}

run_study <- function(st) {
  cs <- cons[study_id == st]
  sids <- unique(c(unlist(lapply(c(cs$numerator_samples, cs$denominator_samples), split_ids)),
                   use[study_id == st, sample_id]))
  sids <- sids[sids %in% use$sample_id]
  owners <- unique(samples[match(sids, sample_id), study_id])
  missing <- owners[!file.exists(file.path("results", owners, FILES[["counts"]]))]
  if (length(missing)) return(list(note = paste("no gene table for", paste(missing, collapse = ", "))))

  meta <- as.data.frame(use[match(sids, sample_id)])
  lev <- unique(meta$group)
  meta$grp <- factor(paste0("L", match(meta$group, lev)), levels = paste0("L", seq_along(lev)))
  rownames(meta) <- sids
  cnt <- get_mat(sids, "counts"); tpm <- get_mat(sids, "tpm")
  if (st %in% NO_LENGTH) {
    dds <- DESeqDataSetFromMatrix(round(cnt), meta, ~ 0 + grp)
  } else {
    len <- get_mat(sids, "lengths")
    len[!is.finite(len) | len <= 0] <- NA
    rm <- rowMeans(len, na.rm = TRUE); rm[!is.finite(rm)] <- 1
    idx <- which(is.na(len), arr.ind = TRUE); len[idx] <- rm[idx[, 1]]
    txi <- list(counts = cnt, abundance = tpm, length = len, countsFromAbundance = "no")
    dds <- suppressMessages(DESeqDataSetFromTximport(txi, meta, ~ 0 + grp))
  }
  replicated_fit <- ncol(dds) > nlevels(meta$grp)
  min_samples <- if (replicated_fit) 2 else 1
  dds <- dds[rowSums(counts(dds) >= MIN_COUNT) >= min_samples, ]
  if (replicated_fit) dds <- DESeq(dds, quiet = TRUE) else dds <- estimateSizeFactors(dds)
  ncounts <- counts(dds, normalized = TRUE)
  gmap <- setNames(as.character(meta$grp), sids)

  out <- list(); summ <- list()
  for (i in seq_len(nrow(cs))) {
    c1 <- cs[i]
    num <- intersect(split_ids(c1$numerator_samples), sids)
    den <- intersect(split_ids(c1$denominator_samples), sids)
    dropped <- setdiff(c(split_ids(c1$numerator_samples), split_ids(c1$denominator_samples)), sids)
    if (!length(num) || !length(den)) {
      summ[[i]] <- data.table(contrast_id = c1$contrast_id, study_id = st, note = "no samples left after exclusions")
      next
    }
    rep_ok <- c1$replicated == "TRUE" && replicated_fit && length(num) >= 2 && length(den) >= 2
    if (rep_ok) {
      ng <- unique(gmap[num]); dg <- unique(gmap[den])
      vec <- setNames(numeric(length(resultsNames(dds))), resultsNames(dds))
      vec[paste0("grp", ng)] <- 1 / length(ng); vec[paste0("grp", dg)] <- -1 / length(dg)
      res <- results(dds, contrast = vec, alpha = PADJ_DEG)
      sh <- suppressMessages(lfcShrink(dds, res = res, type = "ashr", quiet = TRUE))
      d <- data.table(gene_id = rownames(res), base_mean = res$baseMean, log2fc = res$log2FoldChange,
                      lfc_se = res$lfcSE, log2fc_shrunk = sh$log2FoldChange, stat = res$stat,
                      pvalue = res$pvalue, padj = res$padj)
    } else {
      a <- rowMeans(ncounts[, num, drop = FALSE]); b <- rowMeans(ncounts[, den, drop = FALSE])
      d <- data.table(gene_id = rownames(ncounts), base_mean = rowMeans(ncounts[, c(num, den), drop = FALSE]),
                      log2fc = log2((a + 0.5) / (b + 0.5)), lfc_se = NA_real_, log2fc_shrunk = NA_real_,
                      stat = NA_real_, pvalue = NA_real_, padj = NA_real_)
    }
    tp <- tpm[d$gene_id, , drop = FALSE]
    d[, `:=`(contrast_id = c1$contrast_id, study_id = st,
             tpm_numerator_mean = rowMeans(tp[, num, drop = FALSE]),
             tpm_denominator_mean = rowMeans(tp[, den, drop = FALSE]),
             replicated = rep_ok)]
    out[[i]] <- d
    summ[[i]] <- data.table(
      contrast_id = c1$contrast_id, study_id = st, contrast_class = c1$contrast_class,
      replicated = rep_ok, n_numerator = length(num), n_denominator = length(den),
      excluded_libraries = paste(dropped, collapse = ";"), genes_tested = nrow(d),
      n_up = if (rep_ok) d[padj < PADJ_DEG & log2fc >= LFC_DEG, .N] else NA_integer_,
      n_down = if (rep_ok) d[padj < PADJ_DEG & log2fc <= -LFC_DEG, .N] else NA_integer_,
      note = if (replicated_fit) "" else "study has no replicates: log2 ratio of normalised counts")
  }
  list(de = rbindlist(out), summary = rbindlist(summ, fill = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
studies <- if (length(args)) args else sort(unique(cons$study_id))
res <- list(); summ <- list()
for (st in studies) {
  t0 <- Sys.time()
  r <- tryCatch(run_study(st), error = function(e) list(note = paste("error:", conditionMessage(e))))
  if (!is.null(r$note)) {
    message(sprintf("%-20s skipped: %s", st, r$note))
    summ[[st]] <- data.table(study_id = st, note = r$note)
    next
  }
  res[[st]] <- r$de; summ[[st]] <- r$summary
  message(sprintf("%-20s %3d contrasts  %6.1f s", st, nrow(r$summary), as.numeric(Sys.time() - t0, units = "secs")))
}
de <- rbindlist(res, use.names = TRUE)
summary <- rbindlist(summ, fill = TRUE)
setcolorder(de, c("gene_id", "contrast_id", "study_id", "base_mean", "log2fc", "lfc_se", "log2fc_shrunk",
                  "stat", "pvalue", "padj", "tpm_numerator_mean", "tpm_denominator_mean", "replicated"))
if (length(args)) {
  for (st in names(res)) write_parquet(res[[st]], file.path("results", paste0("de_", st, ".parquet")))
} else {
  write_parquet(de, "results/de.parquet")
  fwrite(summary, "results/de_summary.tsv", sep = "\t")
}
message(sprintf("%d gene x contrast rows from %d studies", nrow(de), length(res)))
