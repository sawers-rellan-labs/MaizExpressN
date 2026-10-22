# Data access for the MaizExpressN app. Shiny sources app/R/*.R before app.R (alphabetically),
# so the packages these helpers need are loaded here.
suppressPackageStartupMessages({ library(data.table); library(arrow); library(plotly) })

# Data access for the MaizExpressN app. The data folder is app/data by default; set
# MAIZEXPRESSN_DATA to point elsewhere (built by app/prepare_data.R).

DATA_DIR <- Sys.getenv("MAIZEXPRESSN_DATA", "data")
rp <- function(f) as.data.table(arrow::read_parquet(file.path(DATA_DIR, f)))

GENES     <- rp("genes.parquet")
ALIASES   <- rp("gene_aliases.parquet")
GO_TERMS  <- rp("go_terms.parquet")
CONTRASTS <- rp("contrasts.parquet")
SAMPLES   <- rp("samples.parquet")
STUDIES   <- rp("studies.parquet")
META_TC   <- rp("meta_tissue_class.parquet")
META_TL   <- rp("meta_tissue_level.parquet")
BUILT     <- tryCatch(readLines(file.path(DATA_DIR, "BUILT"))[1], error = function(e) "")
setkey(GENES, gene_id)
setkey(ALIASES, alias_lower)
go_name <- setNames(GO_TERMS$go_name, GO_TERMS$go_id)

CLASSES <- c("deficiency", "starvation", "resupply", "excess", "field_rate", "local", "form")
CLASS_LABEL <- c(deficiency = "Deficiency", starvation = "Starvation", resupply = "Resupply",
                 excess = "Excess", field_rate = "Field rate", local = "Local (split-root)", form = "N form",
                 all_levels = "All N levels")
# Okabe-Ito based, colour-blind safe
CLASS_COL <- c(deficiency = "#0072B2", starvation = "#56B4E9", resupply = "#009E73", excess = "#E69F00",
               field_rate = "#CC79A7", local = "#D55E00", form = "#7A5C3A", all_levels = "#3C4650")
TISSUES <- c("Root", "Leaf", "Shoot", "Ear/kernel", "Whole seedling")
SYSTEMS <- c("Hydroponic", "Solid substrate", "Soil/pots", "Field", "Not reported")
STAGES  <- c("seedling", "vegetative", "reproductive", "Not reported")
LEVEL_CLASSES <- c("deficiency", "starvation", "resupply", "excess", "field_rate", "local")

bucket <- function(g) substr(g, 10, 12)

get_de <- function(gid) {
  f <- file.path(DATA_DIR, "de", paste0(bucket(gid), ".parquet"))
  if (!file.exists(f)) return(data.table())
  as.data.table(arrow::read_parquet(f))[gene_id == gid]
}
get_de_many <- function(ids) {
  rbindlist(lapply(split(ids, bucket(ids)), function(g) {
    f <- file.path(DATA_DIR, "de", paste0(bucket(g[1]), ".parquet"))
    if (!file.exists(f)) return(NULL)
    as.data.table(arrow::read_parquet(f))[gene_id %in% g]
  }))
}
get_tpm <- function(gid) {
  f <- file.path(DATA_DIR, "tpm", paste0(bucket(gid), ".parquet"))
  if (!file.exists(f)) return(data.table())
  as.data.table(arrow::read_parquet(f))[gene_id == gid]
}

# strongest co-expression partners (app/prepare_partners.R); partner IDs are stored as their number
get_partners <- function(gid) {
  f <- file.path(DATA_DIR, "partners", paste0(bucket(gid), ".parquet"))
  if (!file.exists(f)) return(data.table())
  x <- as.data.table(arrow::read_parquet(f))[gene_id == gid]
  x[, partner_id := sprintf("Zm00001eb%06d", partner)][, partner := NULL][]
}

# gene search: exact alias match first, then prefix and description matches
search_genes <- function(q, limit = 50) {
  q <- trimws(tolower(q)); if (nchar(q) < 2) return(character())
  hit <- unique(ALIASES[.(q), gene_id, nomatch = 0])
  if (length(hit) < limit) hit <- unique(c(hit, ALIASES[startsWith(alias_lower, q), gene_id]))
  if (length(hit) < limit && nchar(q) >= 3)
    hit <- unique(c(hit, GENES[grepl(q, tolower(description), fixed = TRUE) |
                                 grepl(q, tolower(full_name), fixed = TRUE), gene_id]))
  head(hit, limit)
}
resolve_tokens <- function(tokens) {
  tokens <- unique(trimws(tokens[nzchar(trimws(tokens))]))
  hits <- lapply(tokens, function(t) unique(ALIASES[.(tolower(t)), gene_id, nomatch = 0]))
  list(genes = unique(unlist(hits)), missing = tokens[lengths(hits) == 0])
}
gene_label <- function(ids) {
  g <- GENES[.(ids)]
  ifelse(nzchar(g$symbol), paste0(g$symbol, " (", g$gene_id, ")"), g$gene_id)
}
