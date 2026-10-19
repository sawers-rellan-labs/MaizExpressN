# M4 validation (19 Oct 2026)

`analysis/de_per_study.R` on 33 studies (Huang2025 still quantifying; Pan2015 deferred), with the
strict library exclusions in `metadata/library_exclusions.tsv` (54 of 645 libraries).
Output: `results/de.parquet` (3.23 M gene × contrast rows) and `results/de_summary.tsv`.

## Contrasts

| | Contrasts |
|---|---|
| Replicated (Wald test, SE, padj) | 89 |
| Unreplicated (log2 ratio of normalised counts, display only) | 33 |
| Skipped: no samples left after library exclusions | 12 |
| Not run yet (Huang2025) or deferred (Pan2015) | 18 + 2 |

Replicated contrasts: median 1,309 DEGs (padj < 0.05 and |log2FC| ≥ 1), range 19 to 10,866.
The largest, Li2026a systemic chloride side (10,866), compares split-root chloride halves with
homogeneous chloride roots; it is display-only (`pool = FALSE`) and looks like a sampling difference
between the two root systems rather than an N effect.

## Positive controls (root contrasts; expected negative log2(lower N / higher N))

| Gene | Class | Contrasts | Negative | Significantly negative | Significantly positive | Median log2FC |
|---|---|---|---|---|---|---|
| nii1 | deficiency | 15 | 15 | 10 | 0 | −2.93 |
| nii1 | starvation | 7 | 7 | 6 | 0 | −2.71 |
| nii1 | resupply | 37 | 36 | 14 | 0 | −2.79 |
| nii1 | local | 8 | 5 | 3 | 0 | −0.66 |
| nnr1 | deficiency | 15 | 12 | 6 | 0 | −0.65 |
| nnr1 | starvation | 7 | 6 | 6 | 0 | −0.84 |
| nnr1 | resupply | 37 | 30 | 14 | 3 | −1.49 |
| nnr1 | local | 8 | 4 | 3 | 1 | −0.30 |
| nrt2 | resupply | 37 | 31 | 9 | 1 | −1.53 |
| nrt2 | deficiency | 15 | 6 | 1 | 1 | +0.24 |
| nrt2 | starvation | 7 | 3 | 1 | 2 | +0.49 |

nrt2 (high-affinity nitrate transporter) is nitrate-inducible (negative in resupply) and mildly
induced by steady N limitation (positive in deficiency and starvation), as expected.

## DEG counts vs the papers (padj < 0.05, |log2FC| ≥ 1)

| Contrast | Paper | Here | Note |
|---|---|---|---|
| Cao2024 wild-type roots, 2 h nitrate | 891 induced, 840 repressed | 730 induced, 865 repressed | matches |
| Ravazzolo2020 NO3 vs −N, 24 h | 1,349 | 229 (471 at padj < 0.05 alone) | ~3–6× lower |
| Ravazzolo2020 NH4 vs −N, 24 h | 1,264 | 55 | NH4_R2 excluded (looked like a nitrate sample); 2 replicates left |
| Li2026a split nitrate vs control nitrate | 661 | 19 | SpNO3_R3 excluded; 2 vs 2 replicates |

Lower counts are expected with an independent pipeline, fewer replicates after strict exclusions,
and a fold-change cut-off; the papers' own thresholds were not always stated.
