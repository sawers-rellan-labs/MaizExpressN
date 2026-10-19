#!/usr/bin/env bash
# M3 step 1: check what Hazel offers before choosing the workflow engine.
# Run on a Hazel login node:  bash hazel_check.sh [/path/to/MaizeExpressN] > hazel_check.txt 2>&1
# Read-only, except one empty test file written to and removed from the data folder.
# Nothing is downloaded (network checks fetch a few bytes of metadata only).

DATA="${1:-/rsstu/users/r/rrellan/sara/MaizeExpressN}"
SOFT=/usr/local/usrapps/maize

section() { printf '\n== %s\n' "$1"; }

section "host and user"
hostname; whoami; date

section "data folder: $DATA"
ls -ld "$DATA" 2>&1
df -h "$DATA" 2>&1 | tail -1
if t=$(mktemp "$DATA/.write_test.XXXXXX" 2>/dev/null); then
  echo "writable: yes"; rm -f "$t"
else
  echo "writable: NO"
fi

section "software folder: $SOFT"
ls -ld "$SOFT" 2>&1
df -h "$SOFT" 2>&1 | tail -1
ls "$SOFT" 2>&1 | head -20

section "scratch candidates"
for d in /share/maize "/share/$USER" /gpfs_common/share02/maize; do
  [ -d "$d" ] && { ls -ld "$d"; df -h "$d" | tail -1; }
done

section "modules (nextflow, apptainer/singularity, java, conda, salmon, sra-tools)"
if command -v module >/dev/null 2>&1; then
  module avail 2>&1 | tr -s ' ' '\n' | grep -iE 'nextflow|apptainer|singularity|java|jdk|conda|mamba|salmon|sra|fastp|snakemake' | sort -u
else
  echo "no 'module' command"
fi

section "commands on PATH"
for c in nextflow apptainer singularity java conda mamba micromamba snakemake salmon prefetch fasterq-dump curl wget; do
  printf '%-14s %s\n' "$c" "$(command -v "$c" 2>/dev/null || echo -)"
done
command -v java >/dev/null 2>&1 && java -version 2>&1 | head -1
command -v apptainer >/dev/null 2>&1 && apptainer --version
command -v singularity >/dev/null 2>&1 && singularity --version

section "LSF"
command -v bsub >/dev/null 2>&1 && echo "bsub: $(command -v bsub)" || echo "bsub: not found"
bqueues 2>/dev/null | head -15

section "network from this login node (HTTP code and seconds)"
for u in \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR6829713&result=read_run&fields=fastq_ftp&format=tsv" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=sra" \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/" \
  "https://sra-pub-run-odp.s3.amazonaws.com/" \
  "https://download.cncb.ac.cn/gsa2/CRA009546/md5sum.txt" \
  "ftp://download.big.ac.cn/gsa2/CRA009546/md5sum.txt" \
  "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.README" \
  "https://github.com" "https://registry.hub.docker.com/v2/" "https://depot.galaxyproject.org/singularity/" ; do
  printf '%s  %s\n' "$(curl -s -o /dev/null -m 25 -w '%{http_code} %{time_total}s' "$u" 2>&1)" "$u"
done

section "done"
echo "Also useful, if you can: run the network block above inside a short batch job"
echo "(e.g. bsub -Is -n 1 -W 10 bash) to see whether compute nodes reach the internet."
