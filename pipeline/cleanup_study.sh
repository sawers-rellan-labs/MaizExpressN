#!/usr/bin/env bash
# Delete one study's FASTQ and Nextflow work files from scratch after its
# results are confirmed. Refuses to run unless the salmon gene table exists.
# The run list and MD5s stay in pipeline/manifests and in the download log
# (copied to the results folder), so everything can be re-downloaded.
#   bash cleanup_study.sh <study>
set -euo pipefail

STUDY="${1:?study_id}"
SCRATCH=/share/maize/rrellan/MaizExpressN
OUT=/rsstu/users/r/rrellan/sara/MaizeExpressN/results/$STUDY
GENES="$OUT/salmon/salmon.merged.gene_counts.tsv"

[[ -s "$GENES" ]] || { echo "no $GENES: not cleaning $STUDY"; exit 1; }
cp -n "$SCRATCH/fastq/$STUDY/download.log" "$OUT/download.log" 2>/dev/null || true
du -sh "$SCRATCH/fastq/$STUDY" "$SCRATCH/work/$STUDY" 2>/dev/null || true
read -r -p "Delete FASTQ and work files for $STUDY? [y/N] " ok
[[ "$ok" == "y" ]] || { echo "kept"; exit 0; }
rm -rf "$SCRATCH/fastq/$STUDY" "$SCRATCH/work/$STUDY"
echo "removed scratch files for $STUDY"
