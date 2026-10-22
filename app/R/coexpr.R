# Co-expression of two genes (Co-expression page).
#
# 1. Across samples, study by study: Pearson correlation of log2(TPM + 1) between the two genes
#    over the libraries of each study x tissue group (the unit), optionally after removing each
#    condition group's mean ("within conditions": the partial correlation given the condition,
#    so a shared response to the N treatment does not count). Fisher z, SE 1/sqrt(df - 1) with
#    df = n - 2 (across all libraries) or n - g - 1 (within, g condition groups; one df per group
#    mean). Units where either gene is barely expressed or too few degrees of freedom are left out.
#    Pooling: DerSimonian-Laird random effects on Fisher z, by tissue group and overall.
# 2. Across contrasts: the two genes' log2FC in every contrast where both were tested (each point is
#    a contrast), with Pearson and Spearman correlations. Contrasts within a study share libraries,
#    so these points are not independent: the correlation is descriptive.

TISSUE_COL <- c(Root = "#7A5C3A", Leaf = "#009E73", Shoot = "#0072B2", `Ear/kernel` = "#E69F00",
                `Whole seedling` = "#CC79A7", All = "#3C4650")
MIN_DF <- 3   # at least 3 residual degrees of freedom (n >= 5 across all libraries)

coexpr_units <- function(ta, tb, samples, within = FALSE, min_tpm = 1) {
  if (!nrow(ta) || !nrow(tb)) return(data.table())
  # every library of the studies where either gene appears; the TPM tables leave out zeros,
  # so a library with no row for a gene has TPM 0
  x <- samples[study_id %in% union(ta$study_id, tb$study_id), .(sample_id, study_id, tissue_group, group)]
  x[ta, a := i.tpm, on = "sample_id"]; x[tb, b := i.tpm, on = "sample_id"]
  x[is.na(a), a := 0]; x[is.na(b), b := 0]
  x[, `:=`(la = log2(a + 1), lb = log2(b + 1))]
  u <- x[, {
    g <- uniqueN(group)
    ra <- if (within) la - ave(la, group) else la - mean(la)
    rb <- if (within) lb - ave(lb, group) else lb - mean(lb)
    df <- if (within) .N - g - 1 else .N - 2
    expressed <- mean(a) >= min_tpm && mean(b) >= min_tpm
    ok <- expressed && df >= MIN_DF && sd(ra) > 0 && sd(rb) > 0
    r <- if (ok) sum(ra * rb) / sqrt(sum(ra^2) * sum(rb^2)) else NA_real_
    rho <- if (ok && !within) suppressWarnings(cor(la, lb, method = "spearman")) else NA_real_
    list(n = .N, groups = g, df = df, mean_tpm_a = mean(a), mean_tpm_b = mean(b), r = r, rho = rho,
         reason = if (ok) "" else if (!expressed) sprintf("mean TPM < %g", min_tpm) else if (df < MIN_DF) "too few libraries" else "no variation")
  }, by = .(study_id, tissue_group)]
  u[, z := atanh(pmin(pmax(r, -0.9999), 0.9999))]
  u[, se := NA_real_]; u[!is.na(r), se := 1 / sqrt(df - 1)]
  u[, `:=`(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se), p = 2 * pnorm(-abs(z / se)))]
  u[order(match(tissue_group, TISSUES), study_id)]
}

# DerSimonian-Laird on Fisher z (same estimator as pool_dl), by tissue group and overall
pool_r <- function(u) {
  u <- u[is.finite(z) & is.finite(se)]
  if (!nrow(u)) return(data.table())
  one <- function(d) {
    w <- 1 / d$se^2; k <- nrow(d)
    zf <- sum(w * d$z) / sum(w); Q <- sum(w * (d$z - zf)^2)
    tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0
    ws <- 1 / (d$se^2 + tau2); zr <- sum(ws * d$z) / sum(ws); s <- sqrt(1 / sum(ws))
    list(k = k, n = sum(d$n), r = tanh(zr), lo = tanh(zr - 1.96 * s), hi = tanh(zr + 1.96 * s),
         p = 2 * pnorm(-abs(zr / s)), tau2 = if (k > 1) tau2 else NA_real_,
         i2 = if (k > 1) (if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0) else NA_real_)
  }
  rbind(u[, one(.SD), by = tissue_group],
        cbind(tissue_group = "All", as.data.table(one(u))))
}

coresponse <- function(da, db, allowed, replicated_only = TRUE) {
  if (!nrow(da) || !nrow(db)) return(data.table())
  x <- merge(da[, .(contrast_id, a = log2fc, a_se = lfc_se, a_padj = padj)],
             db[, .(contrast_id, b = log2fc, b_se = lfc_se, b_padj = padj)], by = "contrast_id")
  x <- allowed[x, on = "contrast_id", nomatch = 0][is.finite(a) & is.finite(b)]
  if (replicated_only) x <- x[replicated == TRUE]
  x
}
coresponse_summary <- function(x) {
  if (nrow(x) < 3) return(data.table())
  s <- function(d) {
    if (nrow(d) < 3 || sd(d$a) == 0 || sd(d$b) == 0)
      return(list(contrasts = nrow(d), studies = uniqueN(d$study_id), pearson = NA_real_, spearman = NA_real_, same_sign = NA_real_))
    list(contrasts = nrow(d), studies = uniqueN(d$study_id), pearson = cor(d$a, d$b),
         spearman = suppressWarnings(cor(d$a, d$b, method = "spearman")),
         same_sign = mean(sign(d$a) == sign(d$b)))
  }
  rbind(cbind(tissue_group = "All", as.data.table(s(x))), x[, s(.SD), by = tissue_group])
}

# ---- plots

coexpr_forest <- function(u, P, la, lb) {
  if (!nrow(u)) return(empty_plot("Neither gene has expression in the filtered studies"))
  d <- copy(u[is.finite(r)])
  if (!nrow(d)) return(empty_plot("No study has both genes expressed with enough libraries"))
  d[, tis := factor(tissue_group, levels = TISSUES)]
  setorder(d, tis, study_id)
  P <- copy(P); P[, tis := factor(tissue_group, levels = c(TISSUES, "All"))]; setorder(P, tis)
  rows <- rbind(d[, .(kind = "study", tissue_group, lab = paste0(study_id, " · ", tissue_group), r, lo, hi,
                      hover = paste0("<b>", study_id, "</b> · ", tissue_group, "<br>r = ", sprintf("%.2f", r),
                                     " [", sprintf("%.2f", lo), ", ", sprintf("%.2f", hi), "]<br>", n, " libraries, ", groups,
                                     " condition groups<br>mean TPM ", signif(mean_tpm_a, 3), " (", la, ") / ", signif(mean_tpm_b, 3), " (", lb, ")"))],
                P[, .(kind = "pooled", tissue_group, lab = paste0("Pooled · ", tissue_group, " (k=", k, ")"), r, lo, hi,
                      hover = paste0("<b>Pooled · ", tissue_group, "</b><br>r = ", sprintf("%.2f", r), " [", sprintf("%.2f", lo), ", ",
                                     sprintf("%.2f", hi), "]<br>", k, ifelse(k == 1, " study", " studies"), ", ", n, " libraries",
                                     ifelse(is.na(i2), "", paste0(" · I² = ", round(100 * i2), "%")), "<br>p = ", fmtp(p)))])
  rows[, row := .I]
  p <- plot_ly(height = max(240, 20 * nrow(rows) + 70))
  for (tg in unique(rows$tissue_group)) for (kd in c("study", "pooled")) {
    q <- rows[tissue_group == tg & kind == kd]
    if (!nrow(q)) next
    col <- TISSUE_COL[[tg]]
    p <- add_trace(p, data = q, x = ~r, y = ~row, type = "scatter", mode = "markers", showlegend = FALSE,
                   error_x = list(type = "data", symmetric = FALSE, array = q$hi - q$r, arrayminus = q$r - q$lo,
                                  color = col, thickness = if (kd == "pooled") 2.2 else 1.4, width = 3),
                   marker = list(color = col, size = if (kd == "pooled") 13 else 8, symbol = if (kd == "pooled") "diamond" else "circle",
                                 line = list(color = "white", width = 1)),
                   text = ~hover, hoverinfo = "text")
  }
  sep <- rows[kind == "pooled", min(row)] - 0.5
  base_layout(p,
    xaxis = list(title = "correlation of log2(TPM + 1)", range = c(-1.05, 1.05), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
    yaxis = list(title = "", tickvals = rows$row, ticktext = rows$lab, autorange = "reversed", gridcolor = "rgba(0,0,0,0)"),
    shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = sep, y1 = sep, line = list(color = "#C9D1D8", width = 1))),
    margin = list(l = 10, r = 10, t = 10, b = 40))
}

coexpr_scatter <- function(x, la, lb, within) {
  if (!nrow(x)) return(empty_plot("Pick a study"))
  x <- copy(x)
  ranks <- sort(unique(suppressWarnings(as.integer(x$N_rank_in_study))))
  pal <- if (length(ranks) > 1) setNames(grDevices::colorRampPalette(c("#B2182B", "#BFA6A0", "#2166AC"))(length(ranks)), ranks) else c(`1` = "#3C4650")
  x[, col := pal[as.character(N_rank_in_study)]]; x[is.na(col), col := "#8A949C"]
  if (within) x[, `:=`(px = la_ - ave(la_, group), py = lb_ - ave(lb_, group))] else x[, `:=`(px = la_, py = lb_)]
  sym <- c(Root = "circle", Leaf = "diamond", Shoot = "square", `Ear/kernel` = "triangle-up", `Whole seedling` = "cross")
  plot_ly(x, x = ~px, y = ~py, type = "scatter", mode = "markers", height = 380,
          symbol = ~tissue_group, symbols = unname(sym[sort(unique(x$tissue_group))]),
          marker = list(color = ~col, size = 10, line = list(color = "white", width = 1)),
          hovertext = ~paste0(sample_id, "<br>", group, " · ", N_label_paper, "<br>", la, " TPM ", signif(a, 3),
                              " · ", lb, " TPM ", signif(b, 3)), hoverinfo = "text") |>
    base_layout(xaxis = list(title = paste0(la, if (within) " (log2 TPM+1, centred within condition)" else " log2(TPM + 1)"), gridcolor = GRID, zeroline = within),
                yaxis = list(title = paste0(lb, if (within) " (centred)" else " log2(TPM + 1)"), gridcolor = GRID, zeroline = within),
                legend = list(orientation = "h", y = -0.2), margin = list(l = 55, r = 10, t = 10, b = 50))
}

coresponse_plot <- function(x, la, lb, sign) {
  if (nrow(x) < 1) return(empty_plot("No contrasts where both genes were tested"))
  x <- copy(x); x[, `:=`(xa = sign * a, xb = sign * b)]
  lim <- min(10, max(1.5, ceiling(max(abs(c(x$xa, x$xb)), na.rm = TRUE) + 0.3)))
  p <- plot_ly(height = 420)
  for (cl in intersect(CLASSES, unique(x$contrast_class))) {
    d <- x[contrast_class == cl]
    p <- add_trace(p, data = d, x = ~xa, y = ~xb, type = "scatter", mode = "markers", name = CLASS_LABEL[[cl]],
                   marker = list(color = CLASS_COL[[cl]], size = 8, opacity = 0.85, line = list(color = "white", width = 0.8)),
                   hovertext = ~paste0("<b>", study_id, "</b> · ", CLASS_LABEL[contrast_class], " · ", tissue_group,
                                       "<br>", numerator_group, " vs ", denominator_group, ifelse(is.na(time_h), "", paste0(" · ", hrs(time_h))),
                                       "<br>", la, " ", fmt(xa), " (padj ", fmtp(a_padj), ")<br>", lb, " ", fmt(xb), " (padj ", fmtp(b_padj), ")"),
                   hoverinfo = "text")
  }
  base_layout(p,
    xaxis = list(title = paste(la, sign_lab(sign)), range = c(-lim, lim), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
    yaxis = list(title = paste(lb, sign_lab(sign)), range = c(-lim, lim), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID, scaleanchor = "x"),
    shapes = list(list(type = "line", x0 = -lim, x1 = lim, y0 = -lim, y1 = lim, line = list(color = "#C9D1D8", dash = "dot", width = 1))),
    legend = list(orientation = "h", y = -0.18), margin = list(l = 55, r = 10, t = 10, b = 60))
}
