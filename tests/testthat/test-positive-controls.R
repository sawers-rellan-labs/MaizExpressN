# Positive controls (HANDOFF M8): nitrate-assimilation genes come out strongly nitrate-induced in
# roots, i.e. negative log2(lower N / higher N), in resupply studies; nitrate reductase and
# nitrite reductase are also lower under steady N deficiency.

meta <- pq("results", "meta_tissue_class.parquet")[filter_set == "all"]
get <- function(g, tissue, cls) meta[gene_id == GENE[[g]] & tissue_group == tissue & contrast_class == cls]

test_that("nitrite reductase (nii1) is strongly nitrate-induced in root resupply studies", {
  m <- get("nii1", "Root", "resupply")
  expect_gte(m$k_studies, 5)
  expect_lt(m$pooled_log2fc, -2)
  expect_lt(m$padj, 0.05)
})

test_that("nitrate reductase (nnr1) and the NRT2 transporter are nitrate-induced in roots", {
  expect_lt(get("nnr1", "Root", "resupply")$pooled_log2fc, -1)
  expect_lt(get("nrt2", "Root", "resupply")$pooled_log2fc, -0.5)
})

test_that("N deficiency lowers nii1 and nnr1 in roots and leaves", {
  for (t in c("Root", "Leaf")) {
    expect_lt(get("nii1", t, "deficiency")$pooled_log2fc, 0, label = paste("nii1", t))
    expect_lt(get("nnr1", t, "deficiency")$pooled_log2fc, 0, label = paste("nnr1", t))
  }
})

test_that("no root contrast shows nii1 significantly higher under low N", {
  de <- pq("results", "de.parquet")[gene_id == GENE[["nii1"]]]
  cons <- tsv("metadata", "contrasts.tsv")
  d <- de[cons[tissue_group == "Root" & contrast_class != "form", .(contrast_id)], on = "contrast_id", nomatch = 0]
  expect_equal(nrow(d[padj < 0.05 & log2fc > 0]), 0)
})
