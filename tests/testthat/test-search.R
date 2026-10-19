# Gene ID search (HANDOFF M5): any identifier or symbol, case-insensitive, every match returned.

al <- pq("annotation", "gene_aliases.parquet")
hits <- function(q) unique(al[alias_lower == tolower(q), gene_id])

test_that("v5, v4 and v3 IDs and the symbol all find nii1", {
  for (q in c("Zm00001eb255880", "Zm00001d018161", "GRMZM2G102959", "nii1")) expect_true(GENE[["nii1"]] %in% hits(q), info = q)
})

test_that("search is case-insensitive", {
  expect_identical(hits("NII1"), hits("nii1"))
  expect_identical(hits("zm00001eb255880"), hits("Zm00001eb255880"))
})

test_that("Arabidopsis names find their maize orthologs, including dotted forms of colon names", {
  expect_true(GENE[["nii1"]] %in% hits("NIR1"))
  expect_true(GENE[["nii1"]] %in% hits("AT2G15620"))
  expect_true(GENE[["nrt2"]] %in% hits("NRT2.1"))
})

test_that("an identifier that maps to several genes returns all of them", {
  expect_gt(length(hits("NIA1")), 1)
  expect_gt(length(hits("gln1.3")), 1)
})

test_that("every quantified gene is in the annotation", {
  g <- pq("annotation", "genes.parquet")
  counts <- fread(Sys.glob(file.path(ROOT, "results", "*", "salmon.merged.gene_counts.tsv"))[1], select = "gene_id")
  expect_true(all(counts$gene_id %in% g$gene_id))
})
