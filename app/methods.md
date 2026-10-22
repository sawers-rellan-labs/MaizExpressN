## Methods

### Studies
Published maize RNA-seq studies that vary nitrogen supply, with raw reads in public repositories (NCBI SRA/GEO, ENA, CNCB GSA). Accessions were resolved from papers and repository records, and every run was mapped to a biological sample with its genotype, tissue, N supply, phase and time. The app shows the studies used in the analysis: those whose data passed quantification and library QC and that give at least one N contrast.

### Quantification
FASTQ files were downloaded on NC State's Hazel cluster with MD5 checks, trimmed with fastp and quantified with salmon (selective alignment, decoy-aware index) against B73 RefGen_v5 gene models Zm00001eb.1 using nf-core/rnaseq 3.27.0, then summarised to genes with tximport. Huang2025 (QuantSeq 3′ tags with UMIs) was aligned with STAR, deduplicated with UMI-tools and quantified with salmon without length correction. Other genotypes are mapped to B73.

### Library QC (strict)
A library is left out if it has fewer than 5 M mapped reads (2 M for 3′ tags), maps below 60% or more than 10 points below other libraries of the same study, tissue and genotype, has a different library type from the rest of its study, or disagrees with its replicates while they agree with each other (when only two replicates disagree, both go). Contrasts that lose replication become display-only.

### Contrasts and sign
Each study is analysed on its own with DESeq2 (`~ 0 + group`), one Wald test per contrast with BH-adjusted p within the contrast. log2FC = log2(lower N / higher N), so positive means higher under N limitation; N-form contrasts are log2(alternative form / nitrate). The sidebar toggle flips the sign to read results as a response to N. Shrunken (ashr) estimates are stored for display; unshrunken log2FC and SE are used for pooling. Unreplicated contrasts show the log2 ratio of mean normalised counts without SE and are never pooled.

A fold change carried by a single library is not tested. In a replicated contrast, a gene is marked “one library” when one library holds more than 90% of its normalised counts across the contrast's libraries, or when it reaches padj < 0.05 but its higher condition, without its top library, is no higher than the lower condition. Such genes (for example 0, 0, 96 reads against 0, 0, 0) can otherwise reach log2FC ±20 with small p-values when DESeq2's Cook's filter misses them. As with Cook's outliers, they get no p-value; their log2FC is still shown (× in the plots) but is not pooled.

### Pooling
Counts are never pooled across studies, only effect sizes. Contrasts that share study, tissue group and class are first combined by inverse variance, so one study counts once; studies are then combined with DerSimonian–Laird random effects (pooled log2FC, 95% CI, p, τ², I²). Estimates without a p-value (one library, or a DESeq2 Cook's outlier) are left out. “All N levels” pools every class except N form. Gene lookup and Gene list pool on the fly under the current filters; Rank genes uses precomputed tables (all genotypes or B73 only).

### Co-expression
Two genes are compared in two ways. **Across samples**: within each study and tissue group, the Pearson correlation of log2(TPM + 1) over its libraries (a library missing from the TPM tables counts as 0). With “within conditions”, each condition group’s mean is removed first, which gives the partial correlation given the condition, so a shared response to the N treatment does not count. Fisher-z standard errors are 1/√(n − 3) across all libraries and 1/√(n − g − 2) within conditions (g condition groups). Units where either gene has a mean TPM below the threshold, or with fewer than three residual degrees of freedom, are left out. Studies are pooled with DerSimonian–Laird random effects on Fisher z, by tissue group and overall. **Across contrasts**: each point is a contrast where both genes were tested, with its two log2FC values; contrasts from one study share libraries and controls, so the Pearson and Spearman correlations there are descriptive. The sidebar’s tissue, growth-system, B73-only and excluded-study filters apply to the sample view; all filters apply to the contrast view.

**Strongest partners** of gene A are precomputed for every gene with the same across-samples statistic (all included libraries, minimum mean TPM 1, no sidebar filters), in three scopes: all study × tissue units, root units and leaf units, each across all libraries and within conditions. A gene has a list in a scope when its mean TPM is at least 1 in at least half of the scope's units; a pair is ranked only when at least half of the units express both genes, and pairs are ranked by pooled r (top 25 in each direction). Two technical sources of co-variation are worth keeping in mind when reading these lists: very long transcripts (about 7 kb and up, against a typical 1.7 kb) tend to rise and fall together with library-to-library differences in RNA integrity, and neighbouring paralogs can share reads; partners within 100 kb are marked “nearby”.

### Gene annotation
MaizeGDB Zm00001eb.1 symbols, names, gene products and GO terms; v3/v4 cross-references from MaizeGDB; Arabidopsis orthologs and gene features from Ensembl Plants. Descriptions marked “predicted” come from MaizeGDB's EnTAP annotation.

### Caveats
- N levels are relative within a study; actual amounts and units are listed with every contrast.
- Non-B73 genotypes are mapped to B73, which biases expression for divergent genes. Use “B73 only” to check.
- Effects from 3′ tag-seq, laser-capture or field designs behave differently from standard mRNA-seq.
