# App data (app/prepare_data.R): bucketed files hold exactly the project results.

test_that("the bucketed DE files hold every row of results/de.parquet", {
  skip_if_not(dir.exists(file.path(ROOT, "app", "data", "de")), "app data not built")
  n <- sum(vapply(Sys.glob(file.path(ROOT, "app", "data", "de", "*.parquet")),
                  function(f) nrow(read_parquet(f, col_select = "gene_id")), numeric(1)))
  used <- pq("app", "data", "contrasts.parquet")$contrast_id
  expect_equal(n, nrow(pq("results", "de.parquet")[contrast_id %in% used]))
})

test_that("the app lists only studies used in the analysis, with no reasons for dropping any", {
  skip_if_not(file.exists(file.path(ROOT, "app", "data", "studies.parquet")), "app data not built")
  st <- pq("app", "data", "studies.parquet"); cm <- pq("app", "data", "contrasts.parquet")
  expect_setequal(st$study_id, unique(cm$study_id))
  expect_true(all(cm$status %in% c("tested", "display only")))
  expect_false(any(c("status", "discard_reason") %in% names(st)))
  not_used <- c("Rosman2024", "Li2023a", "He2016", "Ning2023", "Pan2015", "Ravazzolo2021b")
  expect_length(intersect(not_used, c(st$study_id, cm$study_id, pq("app", "data", "samples.parquet")$study_id)), 0)
  # library QC exclusions are not reported per study or contrast; library counts describe libraries in use
  expect_false("excluded" %in% names(st)); expect_false("excluded_libraries" %in% names(cm))
  smp <- pq("app", "data", "samples.parquet")
  n_used <- smp[, .N, by = study_id][st[, .(study_id, libraries)], on = "study_id"]
  expect_equal(n_used$N, as.integer(n_used$libraries))
  md <- paste(readLines(file.path(ROOT, "app", "methods.md")), collapse = " ")
  for (s in not_used) expect_false(grepl(s, md, fixed = TRUE))
})

test_that("a gene lookup returns the same rows as the full table, quickly", {
  skip_if_not(dir.exists(file.path(ROOT, "app", "data", "de")), "app data not built")
  f <- file.path(ROOT, "app", "data", "de", paste0(substr(GENE[["nii1"]], 10, 12), ".parquet"))
  t0 <- Sys.time(); a <- as.data.table(read_parquet(f))[gene_id == GENE[["nii1"]]]
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 1)
  b <- pq("results", "de.parquet")[gene_id == GENE[["nii1"]]]
  expect_equal(nrow(a), nrow(b))
  expect_setequal(a$contrast_id, b$contrast_id)
})

test_that("excluded libraries never reach the app", {
  skip_if_not(file.exists(file.path(ROOT, "app", "data", "samples.parquet")), "app data not built")
  ex <- tsv("metadata", "library_exclusions.tsv")
  s <- pq("app", "data", "samples.parquet")
  expect_length(intersect(ex$sample_id, s$sample_id), 0)
})

test_that("MaizeMine gene-card tables are present and consistent", {
  d <- file.path(ROOT, "app", "data")
  skip_if_not(file.exists(file.path(d, "pangenes.parquet")), "MaizeMine tables not built")
  rd <- function(f) data.table::as.data.table(arrow::read_parquet(file.path(d, f)))
  pw <- rd("pathways.parquet"); gh <- rd("gwas_hits.parquet"); v <- rd("variants.parquet"); pg <- rd("pangenes.parquet")
  expect_true(all(grepl("^Zm00001eb[0-9]{6}$", c(pw$gene_id, gh$gene_id, v$gene_id, pg$gene_id))))
  expect_true("PWY-381" %in% pw[gene_id == GENE[["nii1"]], pathway_id])           # nitrate reduction II
  expect_true(all(pg$genomes_present >= 1 & pg$genomes_present <= pg$n_genomes))
  expect_equal(anyDuplicated(pg$gene_id), 0L)
  expect_true(all(gh$distance_bp >= 0 & gh$distance_bp <= 10000))
  expect_true(all(v$missense_deleterious <= v$missense))
  expect_true(file.exists(file.path(d, "maizemine_release.json")))
})

test_that("B73 atlas and literature tables are present and consistent", {
  d <- file.path(ROOT, "app", "data")
  skip_if_not(file.exists(file.path(d, "atlas_samples.parquet")), "atlas not built")
  smp <- data.table::as.data.table(arrow::read_parquet(file.path(d, "atlas_samples.parquet")))
  expect_true(all(startsWith(smp$sample, "B73")))
  expect_equal(sort(smp$order), seq_len(nrow(smp)))
  a <- data.table::as.data.table(arrow::read_parquet(file.path(d, "atlas", paste0(substr(GENE[["nii1"]], 10, 12), ".parquet"))))
  a <- a[gene_id == GENE[["nii1"]]]
  expect_equal(nrow(a), nrow(smp))                                 # one value per atlas sample
  expect_true(all(a$sample %in% smp$sample) && all(a$fpkm >= 0, na.rm = TRUE))
  lit <- data.table::as.data.table(arrow::read_parquet(file.path(d, "literature.parquet")))
  expect_true(all(lit$n_genes <= 10) && all(grepl("^[0-9]+$", lit$pmid)))
  expect_gte(lit[gene_id == GENE[["nnr1"]], .N], 3)                  # corn nitrate reductase papers
  expect_equal(anyDuplicated(lit[, .(gene_id, pmid)]), 0L)
})
