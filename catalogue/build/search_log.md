# Search log (16 Oct 2026)

Provenance for `catalogue/maize_N_rnaseq_studies.tsv`. Searches and screening were run by Claude (Cowork) for Rubén Rellán-Álvarez.

## PubMed queries (Title/Abstract unless noted)

| # | Query | Hits |
|---|---|---|
| 1 | (maize OR "Zea mays") AND (nitrogen OR nitrate OR ammonium) AND ("RNA-seq" OR "RNA sequencing" OR "RNA-sequencing" OR transcriptome OR transcriptomic OR transcriptomics) | 193 |
| 2 | ("Zea mays"[MeSH] OR maize) AND ("Nitrogen"[MeSH] OR "Nitrates"[MeSH] OR "Ammonium Compounds"[MeSH]) AND ("Gene Expression Profiling"[MeSH] OR "Transcriptome"[MeSH] OR "Sequence Analysis, RNA"[MeSH] OR "RNA-Seq"[MeSH]) | 82 |
| 3 | (maize OR "Zea mays") AND ("low nitrogen" OR "nitrogen deficiency" OR "nitrogen starvation" OR "nitrogen stress" OR "N deficiency" OR "low N" OR "N starvation" OR "nitrogen limitation" OR "N limitation" OR "nitrogen-deficient") AND (transcriptome OR transcriptomic OR transcriptomics OR "RNA-seq" OR "differentially expressed" OR "gene expression" OR "expression profiling") | 73 |
| 4 | (maize OR "Zea mays" OR corn) AND (nitrogen OR nitrate OR ammonium) AND ("differentially expressed genes" OR DEGs OR WGCNA OR "co-expression" OR lncRNA OR "long non-coding" OR "single-cell" OR "single cell" OR eQTL OR TWAS OR "alternative splicing") | 87 |
| 5 | (maize OR "Zea mays" OR corn) AND ("nitrogen use efficiency" OR "nitrogen fertilizer" OR "nitrogen application" OR "nitrogen rate" OR "nitrogen fertilization" OR "nitrogen supply" OR "nitrogen levels" OR "nitrogen availability" OR "nitrate supply" OR "nitrate availability") AND (transcriptome OR transcriptomic OR transcriptomics OR "RNA-seq" OR "RNA sequencing" OR "differentially expressed" OR "gene expression profiles") | 48 |
| 6 | (maize OR "Zea mays") AND (nitrogen OR nitrate) AND (microRNA OR miRNA OR "small RNA" OR degradome OR microarray) | 32 |
| 7 | (maize OR "Zea mays") AND ("nitrogen form" OR "nitrogen forms" OR "ammonium nutrition" OR "ammonium supply" OR "nitrate and ammonium" OR "ammonium and nitrate" OR urea OR "nitrification inhibitor" OR "nitrification inhibition") AND ("RNA-seq" OR transcriptome OR transcriptomic OR transcriptomics OR "differentially expressed") | 12 |
| 8 | (maize OR "Zea mays") AND ("N-efficient" OR "N efficiency" OR "N supply" OR "N levels" OR "N rates" OR "high nitrogen" OR "nitrogen-efficient" OR "nitrogen efficiency" OR "nitrogen tolerance" OR "low-nitrogen") AND ("RNA-seq" OR transcriptome OR transcriptomic OR "differentially expressed" OR "RNA sequencing") | 33 |
| 9 | (maize OR "Zea mays") AND ("nitrogen remobilization" OR "nitrogen resupply" OR "nitrate induction" OR "nitrate response" OR "nitrate signaling" OR "nitrate signalling" OR "nitrate treatment" OR "nitrate-responsive" OR "N-responsive" OR "nitrogen-responsive") AND ("RNA-seq" OR transcriptome OR transcriptomic OR "differentially expressed" OR "RNA sequencing") | 25 |

Union: 276 unique PMIDs (`pubmed_candidates_276.txt`). The last two queries added only 2 new records, so PubMed recall was close to saturation.

## Screening and extraction

- The 276 records were split into 8 batches and screened by abstract with the rules in `agent_instructions.md`; included papers were read in full text from PubMed Central where available (raw outputs in `agent_out/pubmed_batch_*.tsv`).
- Web searches covered journals PubMed does not index, gene-function papers whose abstracts don't mention RNA-seq, and preprints (`agent_out/web_literature.tsv`), plus repository searches for deposited N datasets (`agent_out/web_datasets.tsv`).
- `build_catalogue.py` merges the batches, applies manual curation (study IDs, accession status, tissue/growth groups) and writes the final tables. Run as `python3 build_catalogue.py <this build folder>`; outputs go to `<build folder>/out/`.
- 60 studies passed screening; 8 were dropped as not core on 16 Oct 2026 (see the catalogue README), leaving 52.

## Known limits

- The PubMed Central full-text service drops hyperlinked accession numbers and some data-availability sections.
- NCBI E-utilities, the ENA Portal API and the Europe PMC annotations API were not reachable from the cloud workspace, and the web-search budget ran out during the first pass, so 27 accessions remain to be resolved.
