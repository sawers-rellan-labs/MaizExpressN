# Vertical slice: quantification QC (17 Oct 2026)

nf-core/rnaseq 3.27.0, fastp + salmon (decoy-aware index, B73 RefGen_v5, Zm00001eb.1; 39,090 genes).
Full outputs are on sara under `MaizeExpressN/results/<study>`. Each study folder here has the
MultiQC report, the DESeq2 QC plots and the merged gene count and TPM tables.

The marker columns use the project sign convention, log2(lower N / higher N), computed from mean CPM.
These are raw ratios as a sanity check, not DESeq2 estimates.

| Study | Samples | Reads mapped | Library | nii1 | nnr1 | nrt2 | Verdict |
|---|---|---|---|---|---|---|---|
| Mu2017 | 4 | 74–86% | IU | down | down | — | OK: gln2 down; gln4, gln5 and asn1 up |
| Ravazzolo2020 | 9 | 79–87% | SR | −4.1 | −2.0 | — | OK |
| Chen2015 | 4 | 63–78% | IU | weak | weak | — | Usable. Field-grown and unreplicated, so the marker response is weak; contrasts are already `replicated = FALSE` (display only, not pooled) |
| Mager2018 | 9 | 84–88% | IU | −6.9 | −2.0 | +1.8 | OK: gln4 +1.4, gln5 +2.5. The −P samples show P-starvation markers (pht2 +6.5, spx1 +3.8), so the labels are not swapped |
| Wang2022 | 21 | 83–90% | IU | −0.8 / −4.9 / −4.2 | −0.4 / −1.9 / −2.9 | −0.2 / −2.8 / −3.4 | OK at 48 h and 72 h; weak at 24 h (values given for 24 / 48 / 72 h) |

## Notes
- All libraries are 100% compatible with the detected library type, and every sample has about 24,000 genes with 10 or more counts.
- Wang2022:
  - Samples cluster by day (S0–S3) before they cluster by N. Within a day, CK and LN replicates are about as similar to each other as to their own group.
  - At 24 h the CK roots had barely induced the nitrate-assimilation genes (nii1 at 6.9 CPM vs 2.4 at S0). The 24 h contrast is therefore weak but real data, and it stays in the pool as approved.
  - There is no consistent replicate pairing between CK and LN, so no `batch_term` is needed.
- Mager2018: replicates cluster cleanly by treatment.
