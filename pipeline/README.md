# pipeline: download and quantification on Hazel (M3)

Hazel's scheduler is **Slurm** (LSF retires end of Fall 2026). Quantification uses
**nf-core/rnaseq** in salmon pseudo-alignment mode with the official nf-core
`hazel` profile (Slurm + Singularity), plus `conf/maizexpressn.config`.
FASTQ downloads use our own MD5-checked script, fed from `metadata/runs.tsv`
(nf-core/fetchngs can't fetch from CNCB GSA, and we already hold URLs and MD5s).

## Locations

| What | Where |
|---|---|
| Pipeline files (copied from the project) | `/rsstu/users/r/rrellan/sara/MaizeExpressN/pipeline` (Mac: `/Volumes/sara/MaizeExpressN/pipeline`) |
| Final outputs (salmon, gene tables, MultiQC) | `/rsstu/users/r/rrellan/sara/MaizeExpressN/results/<study>` |
| FASTQ (temporary), Nextflow work, containers, reference, logs | `/share/maize/rrellan/MaizExpressN/{fastq,work,singularity,reference,logs}` |

## Files

| File | Purpose |
|---|---|
| `make_manifests.py` | (Mac) Writes `manifests/<study>.downloads.tsv` and `<study>.samplesheet.csv` from runs.tsv and samples.tsv. `--slice` = the five vertical-slice studies. |
| `download_fastq.sh` | Downloads one study's FASTQ: resumable, parallel, MD5-checked, re-runnable; log in the study's FASTQ folder. |
| `download_study.sbatch` | Runs `download_fastq.sh` as a job on Hazel's data-transfer node (`xfer` partition, QoS `xfer`). |
| `get_reference.sh` | Fetches the B73 v5 genome and Zm00001eb.1 GFF3 from MaizeGDB and checks MD5s (~660 MB, once). |
| `conf/maizexpressn.config` | nf-core/rnaseq parameters: reference, salmon only, fastp, no bbsplit, keep the built index. |
| `run_study.sh` | Nextflow head job for one study (`sbatch --job-name=<study> run_study.sh <study>`), on the `xfer` node because only it has internet; pipeline steps run as `compute` jobs. Pins Nextflow 25.10.8 (`NXF_VER`) for nf-core/rnaseq 3.27.0; reuses the saved salmon index after the first run. Restores any reference file that the scratch purge deleted, using the copy on sara. |
| `backup_reference.sbatch` | Copies the reference bundle (MaizeGDB files, filtered GTF, transcript FASTA, salmon index) from scratch to `sara/MaizeExpressN/reference` and verifies it by checksum (`sbatch backup_reference.sbatch`). Scratch purges files that haven't been accessed for 30 days. |
| `cleanup_study.sh` | Deletes a study's FASTQ and work files after its gene table exists (asks first). |
| `hazel_check.sh`, `slurm_check.sh`, `net_check.sh` | One-off environment checks (output saved in the share). |

## Running the vertical slice (after approval)

On a Hazel login node:

```bash
cd /rsstu/users/r/rrellan/sara/MaizeExpressN/pipeline
bash get_reference.sh                               # ~660 MB, a minute or two
for s in Mu2017 Ravazzolo2020 Mager2018 Wang2022 Chen2015; do sbatch --job-name=$s download_study.sbatch $s; done
sbatch --job-name=Mu2017 run_study.sh Mu2017      # after Mu2017's download; first run builds the salmon index
# then the other four (they reuse the index)
```

Compute nodes have **no internet**; the `xfer` node does (`net_check/`, 16 Oct 2026). QoS `xfer` caps each user at 4 CPUs / 24 GB on that node, so at most two of our 2-CPU jobs (downloads or head jobs) run at once; Slurm queues the rest.

Slurm on Hazel: account `maize_cpu`; `compute` partition with QoS `normal` (4 days, 256 CPUs/job) for pipeline jobs; `xfer` partition with QoS `xfer` (4 CPUs, 24 GB, 4 days) for downloads; interactive jobs need QoS `short` (2 h, `compute_partners`).

## Known special cases (later studies)

- Ma2021a: R1 and R2 of each library were uploaded as two separate runs. `build_runs.py` pairs them and gives the read-1 run both files, so the normal downloader fetches them from ENA.
- Huang2025: QuantSeq 3' tags with UMIs; needs UMI-aware handling and `--noLengthCorrection` (HANDOFF M3).
- Pan2015: 21-nt DGE tags; deferred.
- Rosman2024: FASTQ derived from BAM submissions (HS_ runs only).
