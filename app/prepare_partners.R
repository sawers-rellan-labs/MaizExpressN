#!/usr/bin/env Rscript
# Strongest co-expression partners of every gene, for the Co-expression page. Run after
# app/prepare_data.R, from the project root (about 15 min on 5 cores, ~4 GB):
#   Rscript app/prepare_partners.R
#
# Same statistic as the page's two-gene view (app/R/coexpr.R), with its defaults: per study x tissue
# unit, Pearson r of log2(TPM + 1) over the unit's libraries (zeros counted), across all libraries or
# within conditions (partial r given the condition group); a unit counts for a pair only when both
# genes have mean TPM >= 1 there and it has >= 3 residual df; Fisher z with SE 1/sqrt(df - 1);
# DerSimonian-Laird pooling. Three tissue scopes: all units, root units, leaf units.
# A gene gets a list in a scope when it has mean TPM >= 1 in at least half of the scope's units, and a
# pair is ranked only when at least half of the units have both genes. Pairs are ranked by pooled r.
#   app/data/partners/<bucket>.parquet   gene_id, scope, view, direction (+/-), rank, partner (number of its ID),
#                                        units, r, lo, hi, i2 (top K per gene, scope, view, direction)
#   results/coexpr_strongest_pairs.tsv   the strongest pairs overall (all units), both views
suppressPackageStartupMessages({ library(data.table); library(arrow); library(parallel) })
MIN_TPM <- 1; MIN_DF <- 3; BLOCK <- 3000; K <- 25; CORES <- 5
SCOPES <- strsplit(Sys.getenv("PARTNER_SCOPES", "All,Root,Leaf"), ",")[[1]]   # subset only for testing
out <- "app/data"
bucket <- function(g) substr(g, 10, 12)

S <- as.data.table(read_parquet(file.path(out, "samples.parquet")))[, .(sample_id, study_id, tissue_group, group)]
S[, unit := paste(study_id, tissue_group, sep = " | ")]
T <- rbindlist(lapply(Sys.glob(file.path(out, "tpm", "*.parquet")), function(f)
  as.data.table(read_parquet(f, col_select = c("gene_id", "sample_id", "tpm")))))[sample_id %in% S$sample_id]
units <- S[, .(n = .N, g = uniqueN(group), tissue_group = tissue_group[1]), by = unit]
m <- T[S, on = "sample_id", nomatch = 0][, .(s = sum(tpm)), by = .(gene_id, unit)][units, on = "unit"][, mean_tpm := s / n]

# gene x library matrices for every gene expressed (mean TPM >= 1) in at least one unit
G_all <- sort(unique(m[mean_tpm >= MIN_TPM, gene_id]))
T <- T[gene_id %in% G_all]
TP <- matrix(0, length(G_all), nrow(S))
TP[cbind(match(T$gene_id, G_all), match(T$sample_id, S$sample_id))] <- T$tpm
rm(T); invisible(gc())
L <- log2(TP + 1)

# per unit: rows scaled to unit norm after centring (overall or by condition group), so r = Z_a . Z_b
prep <- function(G, us, within) {
  gi <- match(G, G_all)
  lapply(us, function(uu) {
    j <- which(S$unit == uu); grp <- S$group[j]
    df <- if (within) length(j) - uniqueN(grp) - 1 else length(j) - 2
    X <- L[gi, j, drop = FALSE]
    if (within) {
      H <- model.matrix(~ 0 + factor(grp)); X <- X - X %*% (H %*% (t(H) / colSums(H)))   # minus group means
    } else X <- X - rowMeans(X)
    ss <- sqrt(rowSums(X^2))
    ok <- rowMeans(TP[gi, j, drop = FALSE]) >= MIN_TPM & ss > 1e-12
    X[!ok, ] <- 0; X[ok, ] <- X[ok, ] / ss[ok]
    list(Z = X, ok = ok, w = df - 1)
  })
}

# pooled r for genes ia x genes ib, then each gene's top K partners in both directions
block_pair <- function(U, ia, ib, kmin) {
  same <- identical(ia, ib)
  zfun <- function(u) atanh(pmin(pmax(tcrossprod(u$Z[ia, , drop = FALSE], u$Z[ib, , drop = FALSE]), -0.9999), 0.9999))
  mk <- function(u) outer(u$ok[ia], u$ok[ib])
  k <- SW <- SW2 <- SWZ <- SWZ2 <- 0
  for (u in U) { M <- mk(u); z <- zfun(u) * M
    k <- k + M; SW <- SW + u$w * M; SW2 <- SW2 + u$w^2 * M; SWZ <- SWZ + u$w * z; SWZ2 <- SWZ2 + u$w * z^2 }
  Q <- SWZ2 - SWZ^2 / SW
  tau2 <- pmax(0, (Q - (k - 1)) / (SW - SW2 / SW)); tau2[k < 2] <- 0
  A <- B <- 0
  for (u in U) { M <- mk(u); ws <- M / (1 / u$w + tau2); A <- A + ws * zfun(u); B <- B + ws }
  zr <- A / B; s <- sqrt(1 / B)
  r <- tanh(zr); r[k < kmin] <- NA
  if (same) diag(r) <- NA
  lo <- tanh(zr - 1.96 * s); hi <- tanh(zr + 1.96 * s)
  i2 <- ifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q), 0)
  tops <- function(R, gi, pj, tr) {
    rbindlist(lapply(c(TRUE, FALSE), function(dec) {
      idx <- lapply(seq_len(nrow(R)), function(i) head(order(R[i, ], decreasing = dec, na.last = NA), K))
      a <- rep(seq_along(idx), lengths(idx)); b <- unlist(idx)
      ij <- if (tr) cbind(b, a) else cbind(a, b)
      data.table(g = gi[a], p = pj[b], direction = if (dec) "+" else "-", units = as.integer(k[ij]),
                 r = r[ij], lo = lo[ij], hi = hi[ij], i2 = i2[ij])
    }))
  }
  rbind(tops(r, ia, ib, FALSE), if (!same) tops(t(r), ib, ia, TRUE))
}

res <- list()
for (sc in SCOPES) for (within in c(FALSE, TRUE)) {
  t0 <- Sys.time()
  uu <- units[(sc == "All" | tissue_group == sc) & (if (within) n - g - 1 else n - 2) >= MIN_DF, unit]
  G <- m[unit %in% uu & mean_tpm >= MIN_TPM, .N, by = gene_id][N >= length(uu) / 2, sort(gene_id)]
  kmin <- max(3L, as.integer(ceiling(length(uu) / 2)))
  U <- prep(G, uu, within)
  blocks <- split(seq_along(G), ceiling(seq_along(G) / BLOCK))
  bp <- CJ(a = seq_along(blocks), b = seq_along(blocks))[a <= b]
  x <- mclapply(seq_len(nrow(bp)), function(i) block_pair(U, blocks[[bp$a[i]]], blocks[[bp$b[i]]], kmin),
                mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(x, inherits, TRUE, "try-error"); if (any(bad)) stop(x[[which(bad)[1]]])
  x <- rbindlist(x)
  x[, rank := frank(if (direction[1] == "+") -r else r, ties.method = "first"), by = .(g, direction)]
  x <- x[rank <= K]
  x[, `:=`(gene_id = G[g], partner_id = G[p], scope = sc, view = if (within) "within" else "across")]
  res[[length(res) + 1]] <- x[, .(gene_id, scope, view, direction, rank = as.integer(rank), partner_id, units,
                                  r = round(r, 3), lo = round(lo, 3), hi = round(hi, 3), i2 = round(i2, 2))]
  message(sprintf("%s %s: %d units, %d genes, pairs need >= %d units, %.0f s", sc, if (within) "within" else "across",
                  length(uu), length(G), kmin, as.numeric(Sys.time() - t0, units = "secs")))
  rm(U, x); invisible(gc())
}
P <- rbindlist(res)

dir.create(file.path(out, "partners"), showWarnings = FALSE)
unlink(Sys.glob(file.path(out, "partners", "*.parquet")))
P[, b := bucket(gene_id)]
# partner IDs are stored as their number (Zm00001eb + 6 digits): a string dictionary of ~20k IDs in every
# file would otherwise dominate the size; get_partners() rebuilds the IDs
stopifnot(all(grepl("^Zm00001eb[0-9]{6}$", P$partner_id)))
for (bb in unique(P$b))
  write_parquet(P[b == bb, .(gene_id, scope, view, direction, rank, partner = as.integer(substr(partner_id, 10, 15)),
                             units, r, lo, hi, i2)], file.path(out, "partners", paste0(bb, ".parquet")), compression = "zstd")
P[, b := NULL]
message(sprintf("partners: %d rows for %d genes in %d buckets", nrow(P), uniqueN(P$gene_id), uniqueN(bucket(P$gene_id))))

# strongest pairs overall (all units): each pair once
if (!"All" %in% SCOPES) quit(save = "no")
A <- as.data.table(read_parquet(file.path(out, "genes.parquet")))[, .(gene_id, symbol, description, chrom, start, end)]
top <- P[scope == "All"][, `:=`(g1 = pmin(gene_id, partner_id), g2 = pmax(gene_id, partner_id))]
top <- unique(top, by = c("view", "g1", "g2"))[, `:=`(gene_id = g1, partner_id = g2)]
top <- rbindlist(lapply(split(top, top$view), function(d) rbind(d[order(-r)][1:1000], d[order(r)][1:200])))
a <- A[match(top$gene_id, gene_id)]; p <- A[match(top$partner_id, gene_id)]
top[, `:=`(symbol_a = a$symbol, symbol_b = p$symbol, description_a = a$description, description_b = p$description,
           distance_kb = ifelse(a$chrom == p$chrom, round(pmax(0, pmax(a$start, p$start) - pmin(a$end, p$end)) / 1e3, 1), NA_real_))]
fwrite(top[, .(view = ifelse(view == "across", "across all libraries", "within conditions"), gene_a = gene_id, symbol_a,
               gene_b = partner_id, symbol_b, units, pooled_r = r, lo, hi, i2, distance_kb, description_a, description_b)],
       "results/coexpr_strongest_pairs.tsv", sep = "\t")
message("strongest pairs written to results/coexpr_strongest_pairs.tsv")
