# Contrast construction and the sign convention (HANDOFF M4): log2FC = log2(lower N / higher N),
# form contrasts log2(alternative form / nitrate).

cons <- tsv("metadata", "contrasts.tsv")
smp <- tsv("metadata", "samples.tsv")
setkey(smp, sample_id)

test_that("every contrast refers to included samples of one tissue group, with disjoint sides", {
  for (i in seq_len(nrow(cons))) {
    c1 <- cons[i]; num <- split_ids(c1$numerator_samples); den <- split_ids(c1$denominator_samples)
    expect_true(length(num) > 0 && length(den) > 0, info = c1$contrast_id)
    expect_true(all(c(num, den) %in% smp$sample_id), info = c1$contrast_id)
    expect_true(all(smp[c(num, den), include] == "TRUE"), info = c1$contrast_id)
    expect_length(intersect(num, den), 0)
    expect_equal(unique(smp[c(num, den), tissue_group]), c1$tissue_group, info = c1$contrast_id)
  }
})

test_that("the lower-N side is the numerator in every poolable N-level contrast", {
  lvl <- cons[contrast_class != "form"]
  for (i in seq_len(nrow(lvl))) {
    c1 <- lvl[i]
    rn <- as.integer(smp[split_ids(c1$numerator_samples), N_rank_in_study])
    rd <- as.integer(smp[split_ids(c1$denominator_samples), N_rank_in_study])
    if (c1$pool == "TRUE") expect_lt(max(rn), min(rd), label = c1$contrast_id)
    else expect_lte(max(rn), min(rd), label = c1$contrast_id)   # display-only systemic split-root sides share an N level
  }
})

test_that("form contrasts put nitrate in the denominator and an alternative form in the numerator", {
  for (i in which(cons$contrast_class == "form")) {
    c1 <- cons[i]
    expect_equal(unique(smp[split_ids(c1$denominator_samples), N_form]), "nitrate", info = c1$contrast_id)
    expect_false("nitrate" %in% smp[split_ids(c1$numerator_samples), N_form], info = c1$contrast_id)
  }
})

test_that("replicated contrasts have at least two samples per side", {
  r <- cons[replicated == "TRUE"]
  expect_true(all(lengths(lapply(r$numerator_samples, split_ids)) >= 2))
  expect_true(all(lengths(lapply(r$denominator_samples, split_ids)) >= 2))
})

test_that("time-resolved contrasts compare matched times unless tagged otherwise", {
  for (i in which(nzchar(cons$time_h))) {
    c1 <- cons[i]
    if (grepl("not_time_matched", c1$tags)) next
    tn <- unique(smp[split_ids(c1$numerator_samples), time_h]); td <- unique(smp[split_ids(c1$denominator_samples), time_h])
    expect_equal(tn, td, info = c1$contrast_id)
  }
})

test_that("DESeq2 log2FC follows the numerator / denominator orientation", {
  de <- pq("results", "de.parquet")
  for (cid in c("Mager2018_deficiency_B73_Root_ControlP", "Cao2024_resupply_ND101_Root_2h", "Wang2019_form_Root_ammonium")) {
    d <- de[contrast_id == cid & base_mean > 50]
    r <- cor(d$log2fc, log2((d$tpm_numerator_mean + 1) / (d$tpm_denominator_mean + 1)))
    expect_gt(r, 0.8, label = cid)
  }
})
