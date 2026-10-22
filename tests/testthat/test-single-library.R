# Fold changes carried by one library (analysis/de_per_study.R, single_library; Rubén, 22 Oct 2026):
# flagged, no p-value, log2FC kept for display, never pooled.

# the rule itself, taken from analysis/de_per_study.R without running the script
rule <- local({
  env <- new.env()
  for (e in parse(file.path(ROOT, "analysis", "de_per_study.R")))
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("SPIKE_SHARE", "SPIKE_PADJ", "is_single_library")) eval(e, env)
  env$is_single_library
})

test_that("the single-library rule flags spikes and leaves consistent responses alone", {
  m <- function(...) matrix(c(...), nrow = 1)
  expect_true(rule(m(0, 0, 96), m(0, 0, 0), 22, 1e-10))          # all reads in one replicate
  expect_true(rule(m(0, 0, 96), m(0, 0, 0), 22, NA))             # ... whatever the p-value
  expect_true(rule(m(1, 2, 300), m(3, 2, 1), 5, 0.2))            # one library > 90% of the counts
  expect_false(rule(m(30, 32, 35), m(0, 0, 0), 20, 1e-8))        # on/off, but in every replicate
  expect_true(rule(m(100, 5), m(6, 6), 3, 0.01))                 # significant, gone without the top library
  expect_false(rule(m(100, 5), m(6, 6), 3, 0.2))                 # ... the leave-one-out part needs padj < 0.05
  expect_true(rule(m(6, 6), m(100, 5), -3, 0.01))                # same, higher in the denominator
  expect_false(rule(m(40, 50, 60), m(10, 12, 9), 2.2, 1e-5))     # ordinary response
  expect_false(rule(m(0, 0), m(0, 0), NA, NA))                   # nothing to judge
})

test_that("single-library genes in results/de.parquet have no p-value and keep their log2FC", {
  de <- pq("results", "de.parquet")
  # Wang2019 shoot, ammonium vs nitrate: each gene's reads sit in one of three replicates
  x <- de[contrast_id == "Wang2019_form_Shoot_ammonium" &
            gene_id %in% c("Zm00001eb378180", "Zm00001eb439650", "Zm00001eb436810")]
  expect_equal(nrow(x), 3)
  expect_true(all(x$single_library))
  expect_true(all(is.na(x$pvalue) & is.na(x$padj)))
  expect_true(all(abs(x$log2fc) > 20))
  expect_false(anyNA(de[replicated == TRUE, single_library]))
  expect_true(all(is.na(de[replicated == FALSE, single_library])))
  expect_true(all(is.na(de[single_library == TRUE, pvalue])))
  share <- de[replicated == TRUE & is.finite(log2fc), mean(single_library)]
  expect_gt(share, 0); expect_lt(share, 0.03)                    # a small minority of tested genes
})

test_that("the counts behind a flagged example sit in one library", {
  cons <- tsv("metadata", "contrasts.tsv")[contrast_id == "Wang2019_form_Shoot_ammonium"]
  libs <- c(split_ids(cons$numerator_samples), split_ids(cons$denominator_samples))
  k <- fread(file.path(ROOT, "results", "Wang2019", "salmon.merged.gene_counts.tsv"))[gene_id == "Zm00001eb378180"]
  cnt <- unlist(k[, libs, with = FALSE])
  expect_equal(sum(cnt > 0), 1)
})

test_that("no estimate without a p-value enters the pooled tables", {
  cons <- tsv("metadata", "contrasts.tsv")
  de <- pq("results", "de.parquet")[, .(gene_id, contrast_id, log2fc, lfc_se, pvalue, single_library)]
  de <- de[cons[pool == "TRUE", .(contrast_id, contrast_class, tissue_group)], on = "contrast_id", nomatch = 0]
  de <- de[is.finite(lfc_se) & lfc_se > 0 & is.finite(log2fc)]
  expect_gt(de[single_library == TRUE, .N], 0)
  want <- de[!is.na(pvalue), .(k_contrasts = .N), by = .(gene_id, tissue_group, contrast_class)]
  meta <- pq("results", "meta_tissue_class.parquet")[filter_set == "all", .(gene_id, tissue_group, contrast_class, k_contrasts)]
  x <- merge(want, meta, by = c("gene_id", "tissue_group", "contrast_class"), all = TRUE)
  expect_equal(x$k_contrasts.x, x$k_contrasts.y)
})
