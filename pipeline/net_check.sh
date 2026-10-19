#!/usr/bin/env bash
# M3 step 1c: exact Nextflow version, and internet access from batch jobs on the
# compute and xfer partitions (interactive srun needs QoS short; batch jobs don't).
# Run on a Hazel login node:  bash net_check.sh
# Submits two 5-minute, 1-core jobs; results land in $OUT within a few minutes.

OUT=/rsstu/users/r/rrellan/sara/MaizeExpressN/net_check
mkdir -p "$OUT"

{
  echo "== nextflow/25 module"
  module load nextflow/25 2>&1 && nextflow -version 2>&1 | grep -iE 'version|build'
  echo "== java"
  java -version 2>&1 | head -1
} > "$OUT/nextflow_version.txt" 2>&1

cat > "$OUT/net_test.sh" <<'EOF'
#!/usr/bin/env bash
echo "node: $(hostname)  partition: $SLURM_JOB_PARTITION  qos: $SLURM_JOB_QOS"
for u in "https://ftp.sra.ebi.ac.uk/vol1/fastq/" \
         "https://download.cncb.ac.cn/gsa2/CRA009546/md5sum.txt" \
         "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.README" \
         "https://github.com" "https://community-cr-prod.seqera.io/v2/" "https://depot.galaxyproject.org/singularity/"; do
  printf '%s  %s\n' "$(curl -s -o /dev/null -m 20 -w '%{http_code} %{time_total}s' "$u")" "$u"
done
command -v sbatch >/dev/null && echo "sbatch available on this node" || echo "no sbatch on this node"
EOF
chmod +x "$OUT/net_test.sh"

sbatch --account=maize_cpu --partition=compute --qos=normal --time=00:05:00 --ntasks=1 --mem=1G \
       --job-name=mxn_net_compute --output="$OUT/compute_%j.txt" "$OUT/net_test.sh"
sbatch --account=maize_cpu --partition=xfer --qos=xfer --time=00:05:00 --ntasks=1 --mem=1G \
       --job-name=mxn_net_xfer --output="$OUT/xfer_%j.txt" "$OUT/net_test.sh"
squeue -u "$USER" -n mxn_net_compute,mxn_net_xfer
echo "results will appear in $OUT"
