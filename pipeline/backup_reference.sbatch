#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=06:00:00
#SBATCH --job-name=mxn_ref_backup
#SBATCH --output=/share/maize/rrellan/MaizExpressN/logs/ref_backup_%j.out
# Copy the reference bundle (MaizeGDB genome and gene models, plus the nf-core
# reference built from them: filtered GTF, transcript FASTA, salmon index) from
# scratch, which purges files untouched for 30 days, to research storage (sara).
# Verifies the copy by checksum. Safe to re-run.
#   sbatch backup_reference.sbatch
set -euo pipefail
SRC=/share/maize/rrellan/MaizExpressN/reference/
DST=/rsstu/users/r/rrellan/sara/MaizeExpressN/reference/
mkdir -p "$DST"
rsync -a "$SRC" "$DST"
# second pass compares checksums; any file listed differs and is re-copied
diff=$(rsync -anc --out-format='%n' "$SRC" "$DST" | grep -v '/$' || true)
if [[ -n "$diff" ]]; then
  echo "checksum differences, re-copying:"; echo "$diff"
  rsync -ac "$SRC" "$DST"
fi
echo "files: $(find "$DST" -type f | wc -l)  size: $(du -sh "$DST" | cut -f1)"
echo "reference bundle backed up to $DST (checksums verified)"
