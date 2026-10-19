# Plotly figures for the MaizExpressN app.

INK <- "#22303A"; MUTED <- "#6B7884"; GRID <- "#E6E9EC"
FONT <- list(family = "IBM Plex Sans, system-ui, sans-serif", size = 12, color = INK)
sign_lab <- function(sign) if (sign > 0) "log2(low N / high N)" else "log2(high N / low N)"
base_layout <- function(p, ...) {
  layout(p, font = FONT, paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
         hoverlabel = list(bgcolor = "white", bordercolor = "#C9D1D8", font = list(color = INK, size = 12)),
         ...) |> config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))
}
hrs <- function(t) ifelse(is.na(t), "", ifelse(t < 1, paste0(round(t * 60), " min"), paste0(signif(t, 3), " h")))
fmt <- function(x, d = 2) ifelse(is.na(x), "–", formatC(x, format = "f", digits = d, flag = "+"))
fmtp <- function(p) ifelse(is.na(p), "–", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 4)))
empty_plot <- function(msg) {
  plot_ly(type = "scatter", mode = "markers") |>
    base_layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                annotations = list(list(text = msg, showarrow = FALSE, font = list(color = MUTED, size = 13))))
}

short_label <- function(e) {
  parts <- e$study_id
  geno <- ifelse(e$genotype_class != "B73" & nzchar(e$genotype) & !grepl(" / ", e$genotype),
                 substr(sub("\\s*\\(.*\\)", "", e$genotype), 1, 16), "")
  ofl <- ifelse(nzchar(e$other_factor_level) & !grepl(" / ", e$other_factor_level) & e$contrast_class != "local",
                substr(sub("\\s*\\(.*\\)", "", e$other_factor_level), 1, 16), "")
  tim <- hrs(e$time_h)
  lab <- paste(parts, geno, ofl, tim, sep = " · ")
  gsub("( · )+", " · ", sub("( · )+$", "", lab))
}

forest_plot <- function(e, sign) {
  if (!nrow(e)) return(empty_plot("No contrasts under the current filters"))
  e <- copy(e)
  e[, `:=`(x = sign * log2fc, ci = fifelse(is.finite(lfc_se), 1.96 * lfc_se, NA_real_),
           cls = factor(contrast_class, levels = CLASSES), tis = factor(tissue_group, levels = TISSUES))]
  setorder(e, tis, cls, study_id, time_h)
  e[, row := .I]
  e[, lab := make.unique(short_label(.SD), sep = " #")]
  e[, hover := paste0("<b>", study_id, "</b> · ", CLASS_LABEL[contrast_class], ifelse(replicated, "", " · unreplicated"),
                      "<br>", tissue_group, " · ", tissue_detail, "<br>", genotype,
                      "<br>", numerator_group, " (", N_numerator, ")<br>vs ", denominator_group, " (", N_denominator, ")",
                      ifelse(is.na(time_h), "", paste0("<br>time ", hrs(time_h))),
                      "<br>", sign_lab(sign), " = ", fmt(x),
                      ifelse(is.na(ci), " (no SE)", paste0(" [", fmt(x - ci), ", ", fmt(x + ci), "]")),
                      "<br>padj = ", fmtp(padj), " · n = ", n_num_used, " vs ", n_den_used)]
  lim <- min(10, max(2, ceiling(max(abs(e$x) + ifelse(is.na(e$ci), 0, e$ci), na.rm = TRUE) + 0.3)))
  p <- plot_ly(height = max(220, 18 * nrow(e) + 60))
  for (cl in intersect(CLASSES, unique(e$contrast_class))) for (rep in c(TRUE, FALSE)) {
    d <- e[contrast_class == cl & replicated == rep]
    if (!nrow(d)) next
    p <- add_trace(p, data = d, x = ~x, y = ~row, type = "scatter", mode = "markers", name = CLASS_LABEL[cl],
                   showlegend = FALSE,
                   error_x = list(type = "data", array = d$ci, color = CLASS_COL[[cl]], thickness = 1.6, width = 3),
                   marker = list(color = if (rep) CLASS_COL[[cl]] else "white", size = 9, symbol = "circle",
                                 line = list(color = CLASS_COL[[cl]], width = 1.6)),
                   text = ~hover, hoverinfo = "text")
  }
  tis_breaks <- e[, .(start = min(row), end = max(row)), by = tis]
  shapes <- lapply(tis_breaks$start[-1] - 0.5, function(y) list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = y, y1 = y,
                                                                line = list(color = "#C9D1D8", width = 1, dash = "dot")))
  shapes <- c(shapes, list(list(type = "line", x0 = 0, x1 = 0, yref = "paper", y0 = 0, y1 = 1, line = list(color = MUTED, width = 1))))
  ann <- lapply(seq_len(nrow(tis_breaks)), function(i) list(x = 1, xref = "paper", xanchor = "right", y = tis_breaks$start[i] - 0.4,
                 yanchor = "bottom", text = paste0("<b>", toupper(tis_breaks$tis[i]), "</b>"), showarrow = FALSE,
                 font = list(size = 11, color = MUTED)))
  base_layout(p,
    xaxis = list(title = sign_lab(sign), range = c(-lim, lim), zeroline = FALSE, gridcolor = GRID),
    yaxis = list(title = "", tickvals = e$row, ticktext = e$lab, autorange = "reversed", gridcolor = "rgba(0,0,0,0)",
                 tickfont = list(size = 11)),
    shapes = shapes, annotations = ann, margin = list(l = 10, r = 10, t = 10, b = 40))
}

pooled_plot <- function(P, sign) {
  if (!nrow(P)) return(empty_plot("Nothing to pool: needs replicated, poolable contrasts"))
  P <- copy(P)
  P[, `:=`(x = sign * est, ci = 1.96 * se, tis = factor(tissue_group, levels = TISSUES),
           cls = factor(contrast_class, levels = c(CLASSES, "all_levels")))]
  setorder(P, tis, cls); P[, row := .I]
  P[, lab := paste0(tissue_group, " · ", CLASS_LABEL[contrast_class])]
  P[, hover := paste0("<b>", lab, "</b><br>", sign_lab(sign), " = ", fmt(x), " [", fmt(x - ci), ", ", fmt(x + ci), "]",
                      "<br>", k, ifelse(k == 1, " study", " studies"), ", ", nc, " contrasts",
                      ifelse(is.na(i2), "", paste0(" · I² = ", round(100 * i2), "% · τ² = ", round(tau2, 2))),
                      "<br>p = ", fmtp(p))]
  lim <- min(8, max(1.5, ceiling(max(abs(P$x) + P$ci) + 0.2)))
  p <- plot_ly(height = max(220, 28 * nrow(P) + 80))
  for (cl in unique(as.character(P$cls))) {
    d <- P[contrast_class == cl]
    p <- add_trace(p, data = d, x = ~x, y = ~row, type = "scatter", mode = "markers", name = CLASS_LABEL[cl],
                   error_x = list(type = "data", array = d$ci, color = CLASS_COL[[cl]], thickness = 2, width = 4),
                   marker = list(symbol = "diamond", size = ifelse(d$k >= 2, 13, 9), color = CLASS_COL[[cl]],
                                 opacity = ifelse(d$k >= 2, 1, 0.45), line = list(color = "white", width = 1.5)),
                   text = ~hover, hoverinfo = "text", showlegend = FALSE)
  }
  base_layout(p,
    xaxis = list(title = paste("pooled", sign_lab(sign)), range = c(-lim, lim), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
    yaxis = list(title = "", tickvals = P$row, ticktext = P$lab, autorange = "reversed", gridcolor = "rgba(0,0,0,0)"),
    annotations = lapply(seq_len(nrow(P)), function(i) list(x = 1, xref = "paper", xanchor = "left", y = P$row[i], text = paste0("k=", P$k[i]),
                                                          showarrow = FALSE, font = list(size = 11, color = MUTED))),
    margin = list(l = 10, r = 50, t = 10, b = 40))
}

expr_plot <- function(x, logy) {
  if (!nrow(x)) return(empty_plot("No expression for this gene in this study"))
  x <- copy(x)
  ranks <- sort(unique(suppressWarnings(as.integer(x$N_rank_in_study))))
  pal <- if (length(ranks) > 1) setNames(grDevices::colorRampPalette(c("#B2182B", "#BFA6A0", "#2166AC"))(length(ranks)), ranks) else c(`1` = "#3C4650")
  x[, col := pal[as.character(N_rank_in_study)]]; x[is.na(col), col := "#8A949C"]
  x[, y := if (logy) pmax(tpm, 0.01) else tpm]
  m <- x[, .(mean = mean(tpm), col = col[1]), by = cond]
  m[, ym := if (logy) pmax(mean, 0.01) else mean]
  lv <- unique(x$cond)
  p <- plot_ly(height = 360) |>
    add_bars(data = m, x = ~factor(cond, levels = lv), y = ~ym, marker = list(color = ~col, opacity = 0.18, line = list(color = ~col, width = 1.4)),
             hovertext = ~paste0(cond, "<br>mean TPM ", signif(mean, 3)), hoverinfo = "text", showlegend = FALSE) |>
    add_markers(data = x, x = ~factor(cond, levels = lv), y = ~y, marker = list(color = ~col, size = 9, line = list(color = "white", width = 1)),
                hovertext = ~paste0(sample_id, "<br>TPM ", signif(tpm, 3), "<br>", N_label_paper), hoverinfo = "text", showlegend = FALSE)
  base_layout(p, bargap = 0.35,
    xaxis = list(title = "", tickangle = -30, tickfont = list(size = 11)),
    yaxis = list(title = "TPM", type = if (logy) "log" else "linear", gridcolor = GRID, dtick = if (logy) 1 else NULL),
    margin = list(l = 50, r = 10, t = 10, b = 110))
}

tc_plot <- function(e, sign) {
  e <- e[!is.na(time_h)]
  if (!nrow(e)) return(empty_plot("No time-resolved contrasts for this gene in this study"))
  e <- copy(e); e[, `:=`(x = sign * log2fc, ci = fifelse(is.finite(lfc_se), 1.96 * lfc_se, 0))]
  e[, line := paste(tissue_group, sub("\\s*\\(.*\\)", "", genotype), CLASS_LABEL[contrast_class], sep = " · ")]
  pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  p <- plot_ly(height = 340)
  for (i in seq_along(unique(e$line))) {
    l <- unique(e$line)[i]; d <- e[line == l][order(time_h)]
    p <- add_trace(p, data = d, x = ~time_h, y = ~x, type = "scatter", mode = "lines+markers", name = l,
                   line = list(color = pal[(i - 1) %% 6 + 1], width = 2),
                   marker = list(color = pal[(i - 1) %% 6 + 1], size = 8, symbol = ifelse(d$replicated, "circle", "circle-open")),
                   error_y = list(type = "data", array = d$ci, color = pal[(i - 1) %% 6 + 1], thickness = 1, width = 3),
                   hovertext = ~paste0(line, "<br>", hrs(time_h), " · ", fmt(x), "<br>padj ", fmtp(padj)), hoverinfo = "text")
  }
  tv <- c(5/60, 15/60, 0.5, 1, 2, 4, 8, 24, 48, 96, 168, 336, 720, 1440)
  tv <- tv[tv >= min(e$time_h) / 1.5 & tv <= max(e$time_h) * 1.5]
  base_layout(p, xaxis = list(title = "time since the N change (log scale)", type = "log", gridcolor = GRID,
                              tickvals = tv, ticktext = hrs(tv)),
              yaxis = list(title = sign_lab(sign), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
              legend = list(orientation = "h", y = -0.25), margin = list(l = 55, r = 10, t = 10, b = 60))
}

heat_plot <- function(Z, hov, sign) {
  if (!nrow(Z)) return(empty_plot("No matched genes"))
  z <- sign * as.matrix(Z[, -1, with = FALSE])
  plot_ly(x = colnames(Z)[-1], y = Z$gene, z = z, type = "heatmap", zmin = -3, zmax = 3, zmid = 0,
          colorscale = list(c(0, "#2166AC"), c(0.5, "#F2F2F0"), c(1, "#B2182B")), xgap = 2, ygap = 2,
          text = hov, hoverinfo = "text", colorbar = list(title = sign_lab(sign), thickness = 12),
          height = max(280, 24 * nrow(Z) + 200)) |>
    base_layout(xaxis = list(tickangle = -40), yaxis = list(autorange = "reversed", tickfont = list(family = "IBM Plex Mono, monospace", size = 11)),
                margin = list(l = 10, r = 10, t = 10, b = 140))
}

volcano_plot <- function(m, sign) {
  if (!nrow(m)) return(empty_plot("No genes for this grouping"))
  m <- copy(m); m[, `:=`(x = sign * pooled_log2fc, y = -log10(pmax(padj, 1e-300)))]
  m[, col := fifelse(padj < 0.05 & abs(pooled_log2fc) >= 1, fifelse(x > 0, "#B2182B", "#2166AC"), "#B8C0C7")]
  plot_ly(m, x = ~x, y = ~y, type = "scattergl", mode = "markers", key = ~gene_id,
          marker = list(color = ~col, size = 5, opacity = 0.75),
          hovertext = ~paste0(label, "<br>", fmt(x), " · padj ", fmtp(padj), "<br>", k_studies, " studies · I² ",
                              ifelse(is.na(i2), "–", paste0(round(100 * i2), "%"))),
          hoverinfo = "text", source = "volcano", height = 420) |>
    event_register("plotly_click") |>
    base_layout(xaxis = list(title = paste("pooled", sign_lab(sign)), zeroline = TRUE, zerolinecolor = MUTED, gridcolor = GRID),
                yaxis = list(title = "−log10 adjusted p", gridcolor = GRID), margin = list(l = 55, r = 10, t = 10, b = 45))
}

class_legend <- function(classes = CLASSES) {
  htmltools::span(class = "legend", lapply(classes, function(cl)
    htmltools::span(class = "lg", htmltools::span(class = "sw", style = paste0("background:", CLASS_COL[[cl]])), CLASS_LABEL[[cl]])),
    htmltools::span(class = "lg", htmltools::span(class = "sw open"), "unreplicated"))
}
