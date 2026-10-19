# Pooling maths (HANDOFF M6/M8): analysis/pool_meta.R must match metafor::rma(method = "DL").

test_that("pooled estimates, SE, tau2 and I2 match metafor on sample genes", {
  skip_if_not_installed("metafor")
  cons <- tsv("metadata", "contrasts.tsv")
  de <- pq("results", "de.parquet")[, .(gene_id, contrast_id, study_id, log2fc, lfc_se)]
  de <- de[cons[pool == "TRUE", .(contrast_id, contrast_class, tissue_group)], on = "contrast_id", nomatch = 0]
  de <- de[is.finite(lfc_se) & lfc_se > 0 & is.finite(log2fc)]
  meta <- pq("results", "meta_tissue_class.parquet")[filter_set == "all"]
  set.seed(1)
  genes <- c(unname(GENE), meta[k_studies >= 3, .N, by = gene_id][sample(.N, 6), gene_id])
  n <- 0
  for (g in genes) for (grp in unique(meta[gene_id == g & k_studies >= 2, paste(tissue_group, contrast_class, sep = "|")])) {
    tg <- strsplit(grp, "|", fixed = TRUE)[[1]]
    s <- de[gene_id == g & tissue_group == tg[1] & contrast_class == tg[2],
            .(y = sum(log2fc / lfc_se^2) / sum(1 / lfc_se^2), v = 1 / sum(1 / lfc_se^2)), by = study_id]
    fit <- metafor::rma(yi = s$y, vi = s$v, method = "DL")
    m <- meta[gene_id == g & tissue_group == tg[1] & contrast_class == tg[2]]
    expect_equal(m$pooled_log2fc, fit$b[1], tolerance = 1e-6, ignore_attr = TRUE, label = paste(g, grp))
    expect_equal(m$se, fit$se, tolerance = 1e-6, label = paste(g, grp))
    expect_equal(m$tau2, fit$tau2, tolerance = 1e-6, label = paste(g, grp))
    expect_equal(m$i2, fit$I2 / 100, tolerance = 1e-6, label = paste(g, grp))
    expect_equal(m$k_studies, fit$k)
    n <- n + 1
  }
  expect_gt(n, 30)
})

test_that("the app's on-the-fly pooling equals the precomputed tables", {
  source(file.path(ROOT, "app", "R", "pool.R"), local = TRUE)
  cons <- tsv("metadata", "contrasts.tsv")
  de <- pq("results", "de.parquet")[gene_id == GENE[["nii1"]]]
  d <- de[cons[pool == "TRUE", .(contrast_id, contrast_class, tissue_group)], on = "contrast_id", nomatch = 0]
  P <- pool_dl(d[is.finite(lfc_se)], c("tissue_group", "contrast_class"))
  meta <- pq("results", "meta_tissue_class.parquet")[filter_set == "all" & gene_id == GENE[["nii1"]]]
  x <- P[meta, on = c("tissue_group", "contrast_class"), nomatch = 0]
  expect_equal(nrow(x), nrow(meta))
  expect_equal(x$est, x$pooled_log2fc, tolerance = 1e-9)
  expect_equal(x$se, x$i.se, tolerance = 1e-9)
})
