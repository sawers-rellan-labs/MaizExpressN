# MaizExpressN: maize gene expression across nitrogen RNA-seq studies.
#   shiny::runApp("app")        (from the project root; data in app/data, built by app/prepare_data.R)
# Set MAIZEXPRESSN_DATA to use a data folder elsewhere.

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(plotly); library(DT); library(data.table); library(arrow)
  library(markdown)   # used by includeMarkdown(); loaded explicitly so deployment tools see the dependency
})
# helpers in R/ (data.R, plots.R, pool.R) are sourced by Shiny before this file

theme <- bs_theme(
  version = 5, bg = "#FFFFFF", fg = "#22303A", primary = "#1F6F5C", secondary = "#6B7884",
  base_font = font_google("IBM Plex Sans", local = FALSE), heading_font = font_google("Spectral", wght = c(500, 600), local = FALSE),
  code_font = font_google("IBM Plex Mono", local = FALSE), "border-radius" = "8px", "font-size-base" = "0.93rem"
)

check <- function(id, label, choices, selected = choices)
  checkboxGroupInput(id, label, choices = choices, selected = selected)

TIME_MAX <- ceiling(max(CONTRASTS$time_h, na.rm = TRUE))
STUDY_CHOICES <- sort(unique(CONTRASTS[status %in% c("tested", "display only"), study_id]))

filters <- sidebar(
  width = 270, open = "desktop",
  radioButtons("sign", "Effect shown as", inline = FALSE,
               choices = c("log2(low N / high N)" = "1", "log2(high N / low N)" = "-1"), selected = "1"),
  tags$small(class = "text-muted", "Positive = higher under N limitation (first option) or induced by N (second)."),
  hr(),
  check("f_tissue", "Tissue group", TISSUES),
  check("f_class", "Contrast class", setNames(CLASSES, CLASS_LABEL[CLASSES])),
  check("f_system", "Growth system", SYSTEMS),
  check("f_stage", "Stage", STAGES),
  input_switch("f_b73", "B73 only", FALSE),
  input_switch("f_rep", "Replicated contrasts only", FALSE),
  sliderInput("f_time", "Time since the N change (h)", min = 0, max = TIME_MAX, value = c(0, TIME_MAX), step = 0.25),
  tags$small(class = "text-muted", "Steady-state contrasts (no time) are always kept."),
  selectizeInput("f_exclude", "Exclude studies", choices = STUDY_CHOICES, multiple = TRUE,
                 options = list(placeholder = "none")),
  actionButton("reset", "Reset filters", class = "btn-sm btn-outline-secondary"),
  textOutput("n_contrasts", inline = TRUE)
)

lookup_page <- nav_panel(
  "Gene lookup", value = "lookup",
  layout_columns(col_widths = 12,
    card(class = "searchcard",
      selectizeInput("gene", NULL, choices = NULL, width = "100%",
                     options = list(placeholder = "Search a gene: symbol (nii1, gln2), v5/v4/v3 ID, Arabidopsis name (NIR1, NRT2.1) or description",
                                    maxOptions = 40, searchField = c("label", "search"), loadThrottle = 250)),
      uiOutput("genecard"))),
  layout_columns(col_widths = c(7, 5),
    card(full_screen = TRUE, card_header("Every contrast", span(class = "hint", "point = log2FC, bar = 95% CI")),
         div(class = "legendrow", class_legend()), plotlyOutput("forest", height = "auto")),
    card(full_screen = TRUE, card_header("Pooled by tissue and contrast class", span(class = "hint", "random effects across studies")),
         plotlyOutput("pooled", height = "auto"), uiOutput("pooled_note"))),
  layout_columns(col_widths = c(6, 6),
    card(full_screen = TRUE, card_header("Expression by condition"),
         layout_columns(col_widths = c(8, 4), selectInput("expr_study", NULL, choices = NULL), input_switch("expr_log", "log scale", TRUE)),
         plotlyOutput("expr", height = "auto"), uiOutput("expr_note")),
    card(full_screen = TRUE, card_header("Time courses"),
         selectInput("tc_study", NULL, choices = NULL), plotlyOutput("tc", height = "auto"))),
  card(card_header("Results table", downloadButton("dl_gene", "Download CSV", class = "btn-sm")), DTOutput("restable"))
)

genelist_page <- nav_panel(
  "Gene list", value = "genelist",
  card(card_header("Gene list"),
       p(class = "lede", "Paste symbols or IDs (v5, v4, v3 or Arabidopsis), separated by spaces, commas or new lines. An identifier that maps to several genes brings in all of them."),
       textAreaInput("gl_text", NULL, rows = 4, width = "100%",
                     value = "nnr1 nnr2 nii1 gln1 gln2 gln4 gln5 gln6 asn1 asn4 nrt2 hb2 cle6 NRT2.1 NLP7 Zm00001d018161"),
       layout_columns(col_widths = c(3, 3, 6), actionButton("gl_go", "Show heatmap", class = "btn-primary"),
                      selectInput("gl_order", NULL, c("as entered" = "input", "by mean effect" = "mean", "by similarity" = "cluster")),
                      uiOutput("gl_matched"))),
  card(full_screen = TRUE, card_header("Genes × pooled tissue and contrast class",
                                       downloadButton("dl_list", "Download CSV", class = "btn-sm")),
       plotlyOutput("heat", height = "auto"))
)

coexpr_page <- nav_panel(
  "Co-expression", value = "coexpr",
  card(class = "searchcard",
    layout_columns(col_widths = c(5, 5, 2),
      selectizeInput("cg_a", "Gene A", choices = NULL, width = "100%",
                     options = list(placeholder = "Search gene A", maxOptions = 40, searchField = c("label", "search"), loadThrottle = 250)),
      selectizeInput("cg_b", "Gene B", choices = NULL, width = "100%",
                     options = list(placeholder = "Search gene B", maxOptions = 40, searchField = c("label", "search"), loadThrottle = 250)),
      div(input_switch("cg_within", "Within conditions", FALSE),
          numericInput("cg_min", "Min mean TPM", 1, min = 0, step = 0.5, width = "100%"))),
    p(class = "hint", "Across samples: correlation of log2(TPM + 1) over the libraries of each study and tissue, pooled across studies. ",
      "\u201cWithin conditions\u201d removes each condition group\u2019s mean first, so a shared response to the N treatment does not count. ",
      "Across contrasts: the two genes\u2019 log2FC in every contrast where both were tested.")),
  card(card_header(textOutput("pt_title", inline = TRUE), span(class = "hint", "precomputed \u00b7 click a row to make it gene B"),
                   downloadButton("dl_pt", "Download CSV", class = "btn-sm")),
    radioButtons("pt_scope", NULL, c("All tissues" = "All", "Root" = "Root", "Leaf" = "Leaf"), inline = TRUE),
    navset_underline(
      nav_panel("Positively co-expressed", div(class = "pt-tables", DTOutput("pt_pos"))),
      nav_panel("Negatively co-expressed", div(class = "pt-tables", DTOutput("pt_neg")))),
    p(class = "hint", "Pooled r across study \u00d7 tissue units (k: units where both genes are expressed), computed as below with min mean TPM 1 and without the sidebar filters; ",
      "follows the \u201cWithin conditions\u201d switch. A pair is ranked only when at least half of the units in the scope express both genes. ",
      "Very long transcripts, and neighbouring paralogs that share reads, can co-vary for technical reasons: \u201cnearby\u201d marks partners within 100 kb.")),
  layout_columns(col_widths = c(6, 6),
    card(full_screen = TRUE, card_header("Co-expression across samples, by study", span(class = "hint", "r with 95% CI · diamonds: pooled (random effects)")),
         plotlyOutput("cx_forest", height = "auto"), uiOutput("cx_note")),
    card(full_screen = TRUE, card_header("Libraries of one study"),
         selectInput("cx_study", NULL, choices = NULL), plotlyOutput("cx_scatter", height = "auto"),
         p(class = "hint", "Colour: red = lowest N level in the study, blue = highest. Shape: tissue."))),
  layout_columns(col_widths = c(7, 5),
    card(full_screen = TRUE, card_header("Co-response across contrasts", span(class = "hint", "each point is a contrast")),
         input_switch("cr_all", "Include unreplicated contrasts", FALSE),
         plotlyOutput("cr_plot", height = "auto"),
         p(class = "hint", "Contrasts from one study share libraries and controls, so the points are not independent: read the correlation as descriptive.")),
    card(card_header("Summary", downloadButton("dl_cx", "Download CSV", class = "btn-sm")),
         tags$h6("Across contrasts"), DTOutput("cr_table"), tags$h6("Across samples, by study"), DTOutput("cx_table")))
)

rank_page <- nav_panel(
  "Rank genes", value = "rank",
  card(card_header("Rank genes", span(class = "hint", "precomputed pooling: only the B73-only filter applies here")),
       layout_columns(col_widths = c(4, 3, 3, 2),
         selectInput("rk_group", "Grouping", choices = NULL),
         sliderInput("rk_k", "Minimum studies", min = 1, max = 20, value = 2, step = 1),
         selectInput("rk_sort", "Rank by", c("adjusted p" = "padj", "pooled log2FC, up" = "up", "pooled log2FC, down" = "down")),
         numericInput("rk_lfc", "|log2FC| ≥", 1, min = 0, step = 0.5)),
       plotlyOutput("volcano", height = "auto"),
       p(class = "hint", "Click a point or a table row to open the gene in Gene lookup.")),
  card(DTOutput("rk_table"))
)

studies_page <- nav_panel(
  "Studies", value = "studies",
  card(card_header("Studies"), uiOutput("study_kpis"), DTOutput("study_table")),
  card(card_header(textOutput("study_title", inline = TRUE)),
       navset_underline(nav_panel("Contrasts", DTOutput("study_contrasts")), nav_panel("Libraries", DTOutput("study_libs"))))
)

methods_page <- nav_panel("Methods", value = "methods", card(class = "methods", includeMarkdown("methods.md")))

ui <- page_navbar(
  id = "nav", title = span(class = "brand", "MaizExpressN"), theme = theme, sidebar = filters, fillable = FALSE,
  header = tags$head(tags$link(rel = "stylesheet", href = "styles.css")),
  lookup_page, genelist_page, coexpr_page, rank_page, studies_page, methods_page,
  nav_spacer(), nav_item(span(class = "built", paste("data", BUILT)))
)

server <- function(input, output, session) {
  sign <- reactive(as.numeric(input$sign))

  # ---- gene search choices (server side)
  alias_txt <- ALIASES[, .(search = paste(unique(alias), collapse = " ")), by = gene_id]
  choices <- GENES[alias_txt, on = "gene_id"][, .(value = gene_id,
    label = paste0(ifelse(nzchar(symbol), paste0(symbol, " · "), ""), gene_id,
                   ifelse(nzchar(description), paste0(" · ", substr(description, 1, 60)), "")),
    search = paste(search, description, full_name))]
  start_gene <- "Zm00001eb255880"
  updateSelectizeInput(session, "gene", choices = choices, server = TRUE, selected = start_gene)
  gene <- reactive({ g <- input$gene; if (is.null(g) || !nzchar(g)) start_gene else g })

  observeEvent(input$reset, {
    updateCheckboxGroupInput(session, "f_tissue", selected = TISSUES)
    updateCheckboxGroupInput(session, "f_class", selected = CLASSES)
    updateCheckboxGroupInput(session, "f_system", selected = SYSTEMS)
    updateCheckboxGroupInput(session, "f_stage", selected = STAGES)
    update_switch("f_b73", value = FALSE); update_switch("f_rep", value = FALSE)
    updateSliderInput(session, "f_time", value = c(0, TIME_MAX))
    updateSelectizeInput(session, "f_exclude", selected = character())
  })

  allowed <- reactive({
    cm <- CONTRASTS[status %in% c("tested", "display only") & tissue_group %in% input$f_tissue &
                    contrast_class %in% input$f_class & system_group %in% input$f_system & stage_class %in% input$f_stage]
    if (isTRUE(input$f_b73)) cm <- cm[genotype_class == "B73"]
    if (isTRUE(input$f_rep)) cm <- cm[replicated == TRUE]
    if (length(input$f_exclude)) cm <- cm[!(study_id %in% input$f_exclude)]
    tr <- input$f_time
    cm[is.na(time_h) | (time_h >= tr[1] & time_h <= tr[2])]
  })
  output$n_contrasts <- renderText(sprintf("%d of %d contrasts", nrow(allowed()), CONTRASTS[status %in% c("tested", "display only"), .N]))

  est <- reactive({
    d <- get_de(gene())
    d[, study_id := NULL]
    allowed()[d, on = "contrast_id", nomatch = 0]
  })

  # ---- gene card
  output$genecard <- renderUI({
    g <- GENES[.(gene())]
    e <- est()
    P <- pool_dl(e[pool == TRUE & replicated == TRUE], c("tissue_group", "contrast_class"))
    best <- if (nrow(P[k >= 2])) P[k >= 2][which.max(abs(est))] else NULL
    tp <- get_tpm(gene())
    top <- if (nrow(tp)) { x <- SAMPLES[tp, on = "sample_id", nomatch = 0][, .(m = mean(tpm)), by = .(study_id, tissue_group, group)]; x[which.max(m)] } else NULL
    gos <- if (nzchar(g$go_terms)) strsplit(g$go_terms, ",")[[1]] else character()
    pred <- grepl("^EnTAP", g$description_source)
    stat <- function(v, l) div(class = "stat", div(class = "v", v), div(class = "l", l))
    div(class = "genecard",
      div(class = "gc-main",
        div(class = "sym", if (nzchar(g$symbol)) g$symbol else g$gene_id),
        div(class = "ids", span(class = "mono", g$gene_id),
            a(href = paste0("https://www.maizegdb.org/gene_center/gene/", g$gene_id), target = "_blank", "MaizeGDB gene page ↗")),
        tags$dl(
          tags$dt("Description"), tags$dd(if (nzchar(g$description)) g$description else span(class = "muted", "none"),
                                          if (pred) span(class = "tag", "predicted (EnTAP)")),
          if (nzchar(g$full_name)) tagList(tags$dt("MaizeGDB name"), tags$dd(g$full_name)),
          tags$dt("Location (v5)"), tags$dd(class = "mono", sprintf("%s:%s–%s (%s)", g$chrom, format(g$start, big.mark = ","), format(g$end, big.mark = ","), g$strand)),
          tags$dt("v4 / v3"), tags$dd(class = "mono", paste(c(if (nzchar(g$v4_ids)) g$v4_ids else "–", if (nzchar(g$v3_ids)) g$v3_ids else "–"), collapse = " / ")),
          tags$dt("Arabidopsis"), tags$dd(if (nzchar(g$ath_best_hit)) tagList(span(class = "mono", g$ath_best_hit), if (nzchar(g$ath_symbol)) paste0(" ", g$ath_symbol),
                                          if (lengths(strsplit(g$ath_orthologs, ",")) > 1) span(class = "muted", paste0("  (orthologs: ", gsub(",", ", ", g$ath_orthologs), ")"))) else span(class = "muted", "no Ensembl ortholog")),
          if (length(gos)) tagList(tags$dt("GO"), tags$dd(class = "small", paste(head(paste0(go_name[gos], " (", gos, ")"), 6), collapse = "; "),
                                                         if (length(gos) > 6) span(class = "muted", sprintf(" … %d more", length(gos) - 6)))))),
      div(class = "stats",
        stat(if (!is.null(best)) fmt(sign() * best$est) else "–",
             if (!is.null(best)) sprintf("largest pooled: %s %s (%d studies)", best$tissue_group, tolower(CLASS_LABEL[best$contrast_class]), best$k) else "pooled (needs 2 studies)"),
        stat(nrow(e), sprintf("contrasts in view (%d studies)", uniqueN(e$study_id))),
        stat(if (!is.null(top)) signif(top$m, 3) else "–", if (!is.null(top)) sprintf("highest mean TPM: %s %s", top$study_id, tolower(top$tissue_group)) else "highest mean TPM"),
        stat(e[replicated == TRUE & !is.na(padj) & padj < 0.05, .N], "contrasts with padj < 0.05")))
  })

  output$forest <- renderPlotly(forest_plot(est(), sign()))

  pooled <- reactive({
    e <- est()[pool == TRUE & replicated == TRUE]
    rbind(pool_dl(e, c("tissue_group", "contrast_class")),
          pool_dl(e[contrast_class %in% LEVEL_CLASSES][, contrast_class := "all_levels"], c("tissue_group", "contrast_class")))
  })
  output$pooled <- renderPlotly(pooled_plot(pooled(), sign()))
  output$pooled_note <- renderUI(p(class = "hint", "Contrasts are combined within each study first. Unreplicated and display-only contrasts are left out. Faint diamonds rest on one study. “All N levels” pools every class except N form."))

  # ---- expression
  observeEvent(gene(), {
    tp <- get_tpm(gene())
    sts <- if (nrow(tp)) tp[, .(m = max(tpm)), by = study_id][order(-m), study_id] else character()
    e <- est(); pick <- if (nrow(e[replicated == TRUE & is.finite(stat)])) e[replicated == TRUE & is.finite(stat)][which.max(abs(stat)), study_id] else sts[1]
    updateSelectInput(session, "expr_study", choices = sort(sts), selected = if (length(pick) && pick %in% sts) pick else sts[1])
    tcs <- CONTRASTS[!is.na(time_h) & status %in% c("tested", "display only"), .(n = uniqueN(time_h)), by = study_id][n >= 3, study_id]
    updateSelectInput(session, "tc_study", choices = tcs, selected = if (!is.null(input$tc_study) && input$tc_study %in% tcs) input$tc_study else "Huang2025")
  })
  expr_data <- reactive({
    req(input$expr_study)
    x <- SAMPLES[study_id == input$expr_study][get_tpm(gene())[study_id == input$expr_study], on = "sample_id", nomatch = 0]
    if (!nrow(x)) return(x)
    x[, cond := paste(group, sep = "")]
    x
  })
  output$expr <- renderPlotly(expr_plot(expr_data(), isTRUE(input$expr_log)))
  output$expr_note <- renderUI({
    s <- STUDIES[study_id == input$expr_study]
    if (!nrow(s)) return(NULL)
    p(class = "hint", sprintf("%s: %s · %s. Salmon TPM per sample; red = lowest N level in the study, blue = highest.",
                              s$study_id, s$tissue, s$growth_system))
  })
  output$tc <- renderPlotly({ req(input$tc_study); tc_plot(est()[study_id == input$tc_study], sign()) })

  # ---- results table
  res_tab <- reactive({
    e <- est()
    if (!nrow(e)) return(data.table())
    e[order(match(tissue_group, TISSUES), match(contrast_class, CLASSES), study_id, time_h),
      .(contrast = contrast_id, study = study_id, class = CLASS_LABEL[contrast_class], tissue = tissue_group, genotype,
        `lower N` = N_numerator, `higher N` = N_denominator, time = hrs(time_h), n = paste0(n_num_used, "/", n_den_used),
        log2FC = round(sign() * log2fc, 3), SE = round(lfc_se, 3), padj = signif(padj, 3),
        `TPM lower N` = round(tpm_numerator_mean, 2), `TPM higher N` = round(tpm_denominator_mean, 2),
        replicated, pooled = pool & replicated, doi = doi_url)]
  })
  output$restable <- renderDT(datatable(res_tab()[, !"doi"], rownames = FALSE, class = "compact stripe",
                                        options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")))
  output$dl_gene <- downloadHandler(filename = function() paste0("MaizExpressN_", gene(), ".csv"),
                                    content = function(f) fwrite(res_tab()[, `log2FC sign` := sign_lab(sign())], f))

  # ---- co-expression
  cx_default <- c(a = "Zm00001eb255880", b = "Zm00001eb176470")   # nii1, nnr1
  updateSelectizeInput(session, "cg_a", choices = choices, server = TRUE, selected = cx_default[["a"]])
  updateSelectizeInput(session, "cg_b", choices = choices, server = TRUE, selected = cx_default[["b"]])
  ga <- reactive({ g <- input$cg_a; if (is.null(g) || !nzchar(g)) cx_default[["a"]] else g })
  gb <- reactive({ g <- input$cg_b; if (is.null(g) || !nzchar(g)) cx_default[["b"]] else g })
  short_name <- function(g) { x <- GENES[.(g)]; if (nrow(x) && nzchar(x$symbol)) x$symbol else g }
  # libraries under the sidebar filters that apply to samples (contrast class and time apply to contrasts only)
  cx_samples <- reactive({
    s <- SAMPLES[tissue_group %in% input$f_tissue & study_id %in% STUDIES[system_group %in% input$f_system, study_id]]
    if (isTRUE(input$f_b73)) s <- s[genotype_class == "B73"]
    if (length(input$f_exclude)) s <- s[!(study_id %in% input$f_exclude)]
    s
  })
  # strongest partners of gene A (precomputed); a click on a row makes the partner gene B
  pt_all <- reactive(get_partners(ga()))
  pt <- reactive({
    p <- pt_all(); if (!nrow(p)) return(p)
    p <- p[scope == input$pt_scope & view == if (isTRUE(input$cg_within)) "within" else "across"][order(direction, rank)]
    g <- GENES[.(p$partner_id)]; a <- GENES[.(ga())]
    d <- ifelse(g$chrom == a$chrom, pmax(0, pmax(g$start, a$start) - pmin(g$end, a$end)) / 1e3, Inf)
    p[, `:=`(symbol = fifelse(is.na(g$symbol), "", g$symbol), description = fifelse(is.na(g$description), "", g$description),
             nearby = fifelse(!is.na(d) & d < 100, sprintf("%.0f kb", d), ""))]
  })
  output$pt_title <- renderText(paste("Strongest partners of", short_name(ga())))
  pt_dt <- function(d) {
    if (!nrow(d)) return(datatable(data.table(note = sprintf(
      "No list for %s here: it needs mean TPM \u2265 1 in at least half of the study \u00d7 tissue units in this scope.", short_name(ga()))),
      rownames = FALSE, selection = "none", options = list(dom = "t")))
    esc <- htmltools::htmlEscape
    short <- fifelse(nchar(d$description) > 60, paste0(substr(d$description, 1, 59), "\u2026"), d$description)
    datatable(d[, .(`#` = rank,
                    gene = esc(fifelse(nzchar(symbol), symbol, "\u2013")), `gene ID` = partner_id,
                    description = sprintf('<span title="%s">%s</span>', esc(description), esc(short)),
                    `r [95% CI]` = sprintf("%.2f [%.2f, %.2f]", r, lo, hi), k = units, nearby)],
              rownames = FALSE, escape = FALSE, selection = "single", class = "compact stripe hover nowrap",
              options = list(dom = "tp", pageLength = 10, ordering = FALSE, scrollX = TRUE))
  }
  output$pt_pos <- renderDT(pt_dt(pt()[direction == "+"]))
  output$pt_neg <- renderDT(pt_dt(pt()[direction == "-"]))
  pick_partner <- function(dir, i) {
    d <- pt()[direction == dir]
    if (length(i) && i <= nrow(d)) updateSelectizeInput(session, "cg_b", choices = choices, server = TRUE, selected = d$partner_id[i])
  }
  observeEvent(input$pt_pos_rows_selected, pick_partner("+", input$pt_pos_rows_selected))
  observeEvent(input$pt_neg_rows_selected, pick_partner("-", input$pt_neg_rows_selected))
  output$dl_pt <- downloadHandler(
    filename = function() sprintf("MaizExpressN_partners_%s_%s_%s.csv", short_name(ga()), input$pt_scope,
                                  if (isTRUE(input$cg_within)) "within" else "across"),
    content = function(f) fwrite(pt()[, .(gene_id, scope, view = fifelse(view == "within", "within conditions", "across all libraries"),
                                          direction = fifelse(direction == "+", "positive", "negative"), rank, partner_id,
                                          partner_symbol = symbol, units, r, ci_low = lo, ci_high = hi, i2, nearby, partner_description = description)], f))

  tpa <- reactive(get_tpm(ga())); tpb <- reactive(get_tpm(gb()))
  cx_u <- reactive(coexpr_units(tpa(), tpb(), cx_samples(), isTRUE(input$cg_within), max(0, input$cg_min %||% 1)))
  cx_P <- reactive(pool_r(cx_u()))
  output$cx_forest <- renderPlotly(coexpr_forest(cx_u(), cx_P(), short_name(ga()), short_name(gb())))
  output$cx_note <- renderUI({
    u <- cx_u(); if (!nrow(u)) return(NULL)
    out <- u[!is.finite(r)]
    why <- if (nrow(out)) paste0(" Left out: ", paste(sprintf("%s %s (%s)", out$study_id, tolower(out$tissue_group), out$reason), collapse = "; "), ".") else ""
    p(class = "hint", sprintf("%d study \u00d7 tissue units shown.%s", u[is.finite(r), .N], why))
  })
  observe({
    u <- cx_u()[is.finite(r)]
    ch <- unique(u$study_id)
    cur <- isolate(input$cx_study)
    updateSelectInput(session, "cx_study", choices = ch,
                      selected = if (!is.null(cur) && cur %in% ch) cur else if (nrow(u)) u[which.max(n), study_id] else character())
  })
  output$cx_scatter <- renderPlotly({
    req(input$cx_study)
    x <- cx_samples()[study_id == input$cx_study]
    x[tpa(), a := i.tpm, on = "sample_id"]; x[tpb(), b := i.tpm, on = "sample_id"]
    x[is.na(a), a := 0]; x[is.na(b), b := 0]
    x[, `:=`(la_ = log2(a + 1), lb_ = log2(b + 1))]
    coexpr_scatter(x, short_name(ga()), short_name(gb()), isTRUE(input$cg_within))
  })
  cr <- reactive(coresponse(get_de(ga()), get_de(gb()), allowed(), !isTRUE(input$cr_all)))
  output$cr_plot <- renderPlotly(coresponse_plot(cr(), short_name(ga()), short_name(gb()), sign()))
  cr_sum <- reactive(coresponse_summary(cr()))
  output$cr_table <- renderDT({
    s <- cr_sum(); if (!nrow(s)) return(datatable(data.table(note = "fewer than 3 shared contrasts"), rownames = FALSE, options = list(dom = "t")))
    datatable(s[, .(tissue = tissue_group, contrasts, studies, Pearson = round(pearson, 2), Spearman = round(spearman, 2),
                    `same direction` = paste0(round(100 * same_sign), "%"))], rownames = FALSE, class = "compact", options = list(dom = "t", scrollX = TRUE))
  })
  output$cx_table <- renderDT({
    u <- cx_u()[is.finite(r)]; P <- cx_P()
    if (!nrow(u)) return(datatable(data.table(note = "no study with both genes expressed"), rownames = FALSE, options = list(dom = "t")))
    t <- rbind(P[, .(unit = paste0("Pooled · ", tissue_group), libraries = n, studies = k, r = round(r, 2),
                     CI = sprintf("%.2f to %.2f", lo, hi), p = signif(p, 2))],
               u[, .(unit = paste(study_id, tissue_group, sep = " · "), libraries = n, studies = 1L, r = round(r, 2),
                     CI = sprintf("%.2f to %.2f", lo, hi), p = signif(p, 2))])
    datatable(t, rownames = FALSE, class = "compact stripe", options = list(dom = "tp", pageLength = 12, scrollX = TRUE))
  })
  output$dl_cx <- downloadHandler(
    filename = function() sprintf("MaizExpressN_coexpression_%s_%s.csv", short_name(ga()), short_name(gb())),
    content = function(f) {
      u <- cx_u(); P <- cx_P(); s <- cr_sum()
      out <- rbind(
        u[, .(view = "across samples", unit = paste(study_id, tissue_group, sep = " · "), n_libraries = n, n_contrasts = NA_integer_, k_studies = 1L,
              r = r, ci_low = lo, ci_high = hi, p = p, spearman = rho, note = reason)],
        P[, .(view = "across samples, pooled", unit = tissue_group, n_libraries = n, n_contrasts = NA_integer_, k_studies = k,
              r = r, ci_low = lo, ci_high = hi, p = p, spearman = NA_real_, note = ifelse(is.na(i2), "", paste0("I2 = ", round(i2, 2))))],
        if (nrow(s)) s[, .(view = "across contrasts", unit = tissue_group, n_libraries = NA_integer_, n_contrasts = contrasts, k_studies = studies,
                           r = pearson, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, spearman = spearman,
                           note = paste0(round(100 * same_sign), "% same direction"))])
      out[, `:=`(gene_a = ga(), gene_b = gb(), within_conditions = isTRUE(input$cg_within))]
      fwrite(out, f)
    })

  # ---- gene list
  gl <- eventReactive(input$gl_go, ignoreNULL = FALSE, {
    r <- resolve_tokens(strsplit(input$gl_text, "[\\s,;]+", perl = TRUE)[[1]])
    if (!length(r$genes)) return(list(r = r, Z = data.table(), hov = NULL))
    d <- get_de_many(head(r$genes, 300)); d[, study_id := NULL]
    d <- allowed()[d, on = "contrast_id", nomatch = 0][pool == TRUE & replicated == TRUE]
    P <- pool_dl(d, c("gene_id", "tissue_group", "contrast_class"))
    if (!nrow(P)) return(list(r = r, Z = data.table(), hov = NULL))
    P[, col := paste(tissue_group, CLASS_LABEL[contrast_class], sep = " · ")]
    cols <- unique(P[order(match(tissue_group, TISSUES), match(contrast_class, CLASSES)), col])
    Z <- dcast(P, gene_id ~ col, value.var = "est")
    H <- dcast(P, gene_id ~ col, value.var = "k")
    Z <- Z[match(intersect(r$genes, Z$gene_id), gene_id)]; H <- H[match(Z$gene_id, gene_id)]
    if (input$gl_order == "mean") { o <- order(-rowMeans(as.matrix(Z[, -1]), na.rm = TRUE)); Z <- Z[o]; H <- H[o] }
    if (input$gl_order == "cluster" && nrow(Z) > 2) { m <- as.matrix(Z[, -1]); m[is.na(m)] <- 0; o <- hclust(dist(m))$order; Z <- Z[o]; H <- H[o] }
    Z <- Z[, c("gene_id", cols), with = FALSE]; H <- H[, c("gene_id", cols), with = FALSE]
    hov <- sapply(cols, function(cc) sprintf("%s<br>%s<br>%s · %s studies", gene_label(Z$gene_id), cc,
                                             fmt(sign() * Z[[cc]]), ifelse(is.na(H[[cc]]), "–", H[[cc]])))
    Z[, gene := gene_label(gene_id)]
    list(r = r, Z = Z, hov = matrix(hov, nrow = nrow(Z)))
  })
  output$gl_matched <- renderUI({ r <- gl()$r
    p(class = "hint", sprintf("%d genes matched", length(r$genes)), if (length(r$missing)) span(class = "warn", paste0(" · not found: ", paste(r$missing, collapse = ", ")))) })
  output$heat <- renderPlotly({ x <- gl(); if (!nrow(x$Z)) return(empty_plot("No pooled estimates for these genes"))
    heat_plot(x$Z[, c("gene", setdiff(names(x$Z), c("gene", "gene_id"))), with = FALSE], x$hov, sign()) })
  output$dl_list <- downloadHandler(filename = function() "MaizExpressN_gene_list.csv",
                                    content = function(f) { x <- gl()$Z; fwrite(x[, c("gene_id", "gene", setdiff(names(x), c("gene_id", "gene"))), with = FALSE], f) })

  # ---- rank
  groups <- rbind(unique(META_TC[, .(tissue_group, contrast_class)]), unique(META_TL[, .(tissue_group, contrast_class)]))
  groups <- groups[order(match(tissue_group, TISSUES), match(contrast_class, c(CLASSES, "all_levels")))]
  gch <- setNames(paste(groups$tissue_group, groups$contrast_class, sep = "|"), paste(groups$tissue_group, CLASS_LABEL[groups$contrast_class], sep = " · "))
  updateSelectInput(session, "rk_group", choices = gch, selected = "Root|resupply")
  rank_data <- reactive({
    req(input$rk_group); tg <- strsplit(input$rk_group, "|", fixed = TRUE)[[1]]
    src <- if (tg[2] == "all_levels") META_TL else META_TC
    m <- src[tissue_group == tg[1] & contrast_class == tg[2] & filter_set == (if (isTRUE(input$f_b73)) "B73" else "all") & k_studies >= input$rk_k]
    m <- GENES[, .(gene_id, symbol, description)][m, on = "gene_id"]
    m[, label := ifelse(nzchar(symbol), paste0(symbol, " · ", gene_id), gene_id)]
    m
  })
  output$volcano <- renderPlotly(volcano_plot(rank_data(), sign()))
  rank_tab <- reactive({
    m <- rank_data()[padj < 0.05 & abs(pooled_log2fc) >= input$rk_lfc]
    m <- switch(input$rk_sort, padj = m[order(padj)], up = m[order(-sign() * pooled_log2fc)], down = m[order(sign() * pooled_log2fc)])
    m[, .(gene_id, symbol, description = substr(description, 1, 60), log2FC = round(sign() * pooled_log2fc, 2),
          CI = paste0("[", fmt(sign() * ci_low), ", ", fmt(sign() * ci_high), "]"), padj = signif(padj, 3), studies = k_studies,
          `up/down` = paste0(k_studies_up, "/", k_studies_down), I2 = round(i2, 2))]
  })
  output$rk_table <- renderDT(datatable(rank_tab(), rownames = FALSE, selection = "single", class = "compact stripe hover",
                                        options = list(pageLength = 20, scrollX = TRUE, dom = "tip")))
  go_gene <- function(g) { updateSelectizeInput(session, "gene", choices = choices, server = TRUE, selected = g); nav_select("nav", "lookup") }
  observeEvent(input$rk_table_rows_selected, go_gene(rank_tab()$gene_id[input$rk_table_rows_selected]))
  observeEvent(event_data("plotly_click", source = "volcano"), { k <- event_data("plotly_click", source = "volcano")$key; if (length(k)) go_gene(k[[1]]) })

  # ---- studies
  output$study_kpis <- renderUI({
    k <- function(v, l) div(class = "stat", div(class = "v", v), div(class = "l", l))
    div(class = "kpis", k(nrow(STUDIES), "studies in the analysis"), k(nrow(SAMPLES), "libraries"),
        k(CONTRASTS[status == "tested", .N], "replicated contrasts tested"),
        k(CONTRASTS[status == "display only", .N], "display-only contrasts"))
  })
  study_tab <- STUDIES[order(study_id),
                       .(study = study_id, year, title = substr(title, 1, 90), tissue = tissue_group, system = system_group,
                         libraries, `median % mapped` = round(median_mapped_pct, 1), contrasts = n_contrasts, doi = doi_url)]
  output$study_table <- renderDT(datatable(study_tab[, !"doi"], rownames = FALSE, selection = list(mode = "single", selected = 1),
                                           class = "compact stripe hover", options = list(pageLength = 15, scrollX = TRUE)))
  sel_study <- reactive({ i <- input$study_table_rows_selected; if (length(i)) study_tab$study[i] else study_tab$study[1] })
  output$study_title <- renderText({ s <- STUDIES[study_id == sel_study()]; sprintf("%s · %s (%s)", s$study_id, s$title, s$year) })
  output$study_contrasts <- renderDT(datatable(CONTRASTS[study_id == sel_study(),
    .(contrast = contrast_id, class = CLASS_LABEL[contrast_class], tissue = tissue_detail, genotype, `lower N` = N_numerator,
      `higher N` = N_denominator, time = hrs(time_h), status, n = paste0(n_num_used, "/", n_den_used), up = n_up, down = n_down)], rownames = FALSE, class = "compact stripe", options = list(dom = "tip", scrollX = TRUE)))
  output$study_libs <- renderDT(datatable(SAMPLES[study_id == sel_study(), .(sample = sample_id, condition = N_label_paper, tissue = tissue_detail,
    genotype, `mapped reads (M)` = round(reads_mapped / 1e6, 1), `% mapped` = round(percent_mapped, 1))],
    rownames = FALSE, class = "compact stripe", options = list(dom = "tip", scrollX = TRUE, pageLength = 25)))
}

shinyApp(ui, server)
