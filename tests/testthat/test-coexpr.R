# Co-expression page (app/R/coexpr.R): correlations, within-condition partial correlations,
# pooling and co-response, checked against base R, lm() and metafor.

app_env <- function() {
  skip_if_not(dir.exists(file.path(ROOT, "app", "data", "tpm")), "app data not built")
  e <- new.env()
  old <- Sys.getenv("MAIZEXPRESSN_DATA"); Sys.setenv(MAIZEXPRESSN_DATA = file.path(ROOT, "app", "data"))
  on.exit(Sys.setenv(MAIZEXPRESSN_DATA = old))
  for (f in c("data.R", "plots.R", "pool.R", "coexpr.R")) sys.source(file.path(ROOT, "app", "R", f), envir = e)
  e
}

test_that("per-study r is the correlation of log2(TPM + 1), zeros included", {
  e <- app_env()
  ta <- e$get_tpm(GENE[["nii1"]]); tb <- e$get_tpm(GENE[["nnr1"]])
  u <- e$coexpr_units(ta, tb, e$SAMPLES, within = FALSE, min_tpm = 1)
  for (i in which(is.finite(u$r))[1:5]) {
    s <- e$SAMPLES[study_id == u$study_id[i] & tissue_group == u$tissue_group[i], .(sample_id, group)]
    s[ta, a := i.tpm, on = "sample_id"]; s[tb, b := i.tpm, on = "sample_id"]
    s[is.na(a), a := 0]; s[is.na(b), b := 0]
    expect_equal(u$n[i], nrow(s))
    expect_equal(u$r[i], cor(log2(s$a + 1), log2(s$b + 1)))
  }
})

test_that("within-condition r equals the partial correlation given the condition group", {
  e <- app_env()
  ta <- e$get_tpm(GENE[["nii1"]]); tb <- e$get_tpm(GENE[["nnr1"]])
  u <- e$coexpr_units(ta, tb, e$SAMPLES, within = TRUE, min_tpm = 1)
  for (i in which(is.finite(u$r))[1:4]) {
    s <- e$SAMPLES[study_id == u$study_id[i] & tissue_group == u$tissue_group[i], .(sample_id, group)]
    s[ta, a := i.tpm, on = "sample_id"]; s[tb, b := i.tpm, on = "sample_id"]
    s[is.na(a), a := 0]; s[is.na(b), b := 0]
    pr <- cor(resid(lm(log2(a + 1) ~ group, s)), resid(lm(log2(b + 1) ~ group, s)))
    expect_equal(u$r[i], pr)
    expect_equal(u$se[i], 1 / sqrt(nrow(s) - uniqueN(s$group) - 2))   # partial-correlation Fisher z SE
  }
})

test_that("pooled correlation matches metafor DerSimonian-Laird on Fisher z", {
  skip_if_not_installed("metafor")
  e <- app_env()
  u <- e$coexpr_units(e$get_tpm(GENE[["nii1"]]), e$get_tpm(GENE[["nnr1"]]), e$SAMPLES, FALSE, 1)
  P <- e$pool_r(u); ok <- u[is.finite(r)]
  m <- metafor::rma(yi = z, sei = se, data = ok, method = "DL")
  expect_equal(P[tissue_group == "All", r], tanh(as.numeric(m$b)))
  expect_equal(P[tissue_group == "All", lo], tanh(m$ci.lb), tolerance = 1e-4)   # app uses 1.96, metafor qnorm(0.975)
  expect_equal(P[tissue_group == "All", k], nrow(ok))
  mr <- metafor::rma(yi = z, sei = se, data = ok[tissue_group == "Root"], method = "DL")
  expect_equal(P[tissue_group == "Root", r], tanh(as.numeric(mr$b)))
})

test_that("a gene with itself correlates perfectly and nii1 co-expresses with nnr1 in roots", {
  e <- app_env()
  t1 <- e$get_tpm(GENE[["nii1"]])
  u <- e$coexpr_units(t1, t1, e$SAMPLES, FALSE, 1)
  expect_true(all(abs(u[is.finite(r), r] - 1) < 1e-9))
  P <- e$pool_r(e$coexpr_units(t1, e$get_tpm(GENE[["nnr1"]]), e$SAMPLES, FALSE, 1))
  expect_gt(P[tissue_group == "Root", lo], 0.5)   # nitrate assimilation genes move together
})

test_that("co-response pairs the two genes contrast by contrast", {
  e <- app_env()
  da <- e$get_de(GENE[["nii1"]]); db <- e$get_de(GENE[["nnr1"]])
  x <- e$coresponse(da, db, e$CONTRASTS[status %in% c("tested", "display only")], replicated_only = TRUE)
  expect_true(all(x$replicated))
  i <- x$contrast_id[1]
  expect_equal(x[contrast_id == i, a], da[contrast_id == i, log2fc])
  expect_equal(x[contrast_id == i, b], db[contrast_id == i, log2fc])
  s <- e$coresponse_summary(x)
  expect_equal(s[tissue_group == "All", pearson], cor(x$a, x$b))
})

test_that("precomputed partner lists match the two-gene view and are well formed", {
  e <- app_env()
  skip_if_not(dir.exists(file.path(ROOT, "app", "data", "partners")), "partner lists not built")
  p <- e$get_partners(GENE[["nii1"]])
  expect_setequal(unique(p$scope), c("All", "Root", "Leaf"))
  expect_false(any(p$partner_id == GENE[["nii1"]]))
  for (s in split(p, by = c("scope", "view", "direction"))) {
    expect_equal(s$rank[order(s$rank)], seq_len(nrow(s)))
    s <- s[order(rank)]
    expect_true(all(diff(if (s$direction[1] == "+") -s$r else s$r) >= -1e-9))   # ranked by pooled r
    expect_lte(nrow(s), 25)
  }
  # the top entries agree with coexpr_units + pool_r on the same libraries
  chk <- p[rank <= 2]
  for (i in seq_len(nrow(chk))) {
    smp <- if (chk$scope[i] == "All") e$SAMPLES else e$SAMPLES[tissue_group == chk$scope[i]]
    u <- e$coexpr_units(e$get_tpm(GENE[["nii1"]]), e$get_tpm(chk$partner_id[i]), smp, chk$view[i] == "within", 1)
    P <- e$pool_r(u)[tissue_group == "All"]
    expect_equal(P$k, chk$units[i])
    expect_equal(P$r, chk$r[i], tolerance = 1e-3)
    expect_equal(P$lo, chk$lo[i], tolerance = 1e-3)
  }
  # nitrate assimilation partners: nii1 with nii2 (paralog) among its strongest root partners across libraries
  expect_true("Zm00001eb193660" %in% p[scope == "Root" & view == "across" & direction == "+", partner_id])
})
