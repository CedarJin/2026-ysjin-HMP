#!/usr/bin/env bash
#SBATCH --job-name=hmp_gus_download
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/hmp_gus_download.%j.out
#SBATCH --error=logs/hmp_gus_download.%j.err

set -euo pipefail

cd /home/jys0914/2026-ysjin-HMP

echo "[$(date -Is)] Job started on $(hostname)"
echo "[$(date -Is)] Working directory: $(pwd)"

scripts/download_hmp_gus_data.sh

echo "[$(date -Is)] Job finished"
