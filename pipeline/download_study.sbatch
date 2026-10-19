#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/share/maize/rrellan/MaizExpressN/logs/download_%x_%j.out
# Download one study's FASTQ on Hazel's data-transfer node (xfer partition).
#   sbatch --job-name=<study> download_study.sbatch <study> [connections, default 6]
# Re-submitting is safe: verified files are skipped, partial files resume.
# QoS xfer allows 4 CPUs / 24 GB per user in total: 1 CPU per download (curl is I/O-bound)
# leaves room for two downloads plus a 2-CPU Nextflow head job.
set -uo pipefail
STUDY="${1:?study_id}"
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressN/pipeline
bash "$PIPE/download_fastq.sh" "$PIPE/manifests/$STUDY.downloads.tsv" "${2:-6}"
