#!/usr/bin/env bash
# Fetch the B73 RefGen_v5 genome and gene models from MaizeGDB into Hazel scratch
# and check them against MaizeGDB's published MD5s. About 660 MB. Run once.
#   bash get_reference.sh [dest_dir]
set -euo pipefail

DEST="${1:-/share/maize/rrellan/MaizExpressN/reference}"
BASE=https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0
FILES=(Zm-B73-REFERENCE-NAM-5.0.fa.gz
       Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
       Zm-B73-REFERENCE-NAM-5.0.README)

mkdir -p "$DEST"
cd "$DEST"
curl -fsSL -o Zm-B73-REFERENCE-NAM-5.0.md5 "$BASE/Zm-B73-REFERENCE-NAM-5.0.md5"
for f in "${FILES[@]}"; do
  echo "fetching $f"
  curl -fsSL --retry 5 -C - -o "$f" "$BASE/$f"
done
grep -E " ($(IFS='|'; echo "${FILES[*]}"))\$" Zm-B73-REFERENCE-NAM-5.0.md5 | md5sum -c -
echo "reference ready in $DEST"
