# Pooling on the fly for the current filters: within a study, contrasts sharing the grouping
# are combined by inverse variance; across studies, DerSimonian-Laird random effects.
# Same maths as analysis/pool_meta.R (checked against metafor in tests/test_pooling.R).

pool_dl <- function(d, keys) {
  d <- d[is.finite(lfc_se) & lfc_se > 0 & is.finite(log2fc)]
  if (!nrow(d)) return(data.table())
  by_s <- c(keys, "study_id"); d <- copy(d); d[, w := 1 / lfc_se^2]
  s <- d[, .(y = sum(w * log2fc) / sum(w), v = 1 / sum(w), nc = .N), by = by_s]
  s[, w := 1 / v]
  a <- s[, .(k = .N, W = sum(w), Sw2 = sum(w^2), Swy = sum(w * y), Swy2 = sum(w * y^2), nc = sum(nc)), by = keys]
  a[, Q := pmax(0, Swy2 - Swy^2 / W)]
  a[, tau2 := fifelse(k > 1, pmax(0, (Q - (k - 1)) / (W - Sw2 / W)), 0)]
  s <- s[a[, c(keys, "tau2"), with = FALSE], on = keys]
  s[, ws := 1 / (v + tau2)]
  r <- s[, .(est = sum(ws * y) / sum(ws), se = sqrt(1 / sum(ws))), by = keys]
  a <- a[r, on = keys]
  a[, `:=`(p = 2 * pnorm(-abs(est / se)),
           i2 = fifelse(k > 1, fifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q), 0), NA_real_),
           tau2 = fifelse(k > 1, tau2, NA_real_))]
  a[, c(keys, "k", "nc", "est", "se", "p", "tau2", "i2"), with = FALSE]
}
