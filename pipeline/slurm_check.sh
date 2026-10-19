#!/usr/bin/env bash
# M3 step 1b: Slurm details on Hazel (LSF is being retired; Slurm is production).
# Run on a Hazel login node:
#   bash slurm_check.sh > /rsstu/users/r/rrellan/sara/MaizeExpressN/slurm_check.txt 2>&1
# Creates the scratch folders under /share/maize/rrellan/MaizExpressN.
# Submits one 5-minute, 1-core test job to see whether compute nodes reach the internet.

SCRATCH=/share/maize/rrellan/MaizExpressN
section() { printf '\n== %s\n' "$1"; }

section "scratch folders"
mkdir -p "$SCRATCH"/{fastq,work,logs,reference,singularity} && ls -ld "$SCRATCH" "$SCRATCH"/*
df -h "$SCRATCH" | tail -1

section "Slurm partitions (sinfo -s)"
sinfo -s 2>&1 | head -30

section "my associations (account / partition / qos)"
sacctmgr -nP show assoc user="$USER" format=account,partition,qos 2>&1 | head -20

section "sqos (Hazel helper), if present"
command -v sqos >/dev/null && sqos 2>&1 | head -30 || echo "no sqos"

section "nextflow modules (default path and BRC path)"
module avail nextflow 2>&1 | tr -s ' ' '\n' | grep -i nextflow
module load nextflow/26.04.3 2>&1 && nextflow -version 2>&1 | grep -i version || echo "nextflow/26.04.3 not loadable from the default module path"
module avail singularity apptainer 2>&1 | tr -s ' ' '\n' | grep -iE 'singularity|apptainer'

section "BRC shared container library used by the nf-core hazel profile"
ls -ld /gpfs_backup/bioinfo_data/nf-core_ood/container_cache 2>&1
ls /gpfs_backup/bioinfo_data/nf-core_ood/container_cache 2>/dev/null | grep -iE 'salmon|fastp|fastqc|multiqc|tximeta|gffread' | head -20

section "compute-node internet (srun, 5 min, 1 core)"
srun --partition=compute --time=00:05:00 --ntasks=1 --mem=1G bash -c '
  hostname
  for u in "https://ftp.sra.ebi.ac.uk/vol1/fastq/" "https://download.cncb.ac.cn/gsa2/CRA009546/md5sum.txt" "https://github.com"; do
    printf "%s  %s\n" "$(curl -s -o /dev/null -m 20 -w "%{http_code} %{time_total}s" "$u")" "$u"
  done' 2>&1

section "done"
