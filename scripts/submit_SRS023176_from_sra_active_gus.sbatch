#!/bin/bash
#SBATCH --job-name=SRS023176_sra_qc
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/SRS023176_sra_qc.%j.out
#SBATCH --error=logs/SRS023176_sra_qc.%j.err

# Phase 1 only: reconstruct both WGS runs from SRA and validate every FASTQ pair.
# This job deliberately does not run HMP preprocessing or active GUS.

set -euo pipefail

PROJECT=/home/jys0914/2026-ysjin-HMP
GUS_ENV=/home/jys0914/2026-ysjin-gus/envs/gus-map
cd "$PROJECT"

mkdir -p logs .cache "data/HMP_GUS/SRA/SRS023176/raw"
export PATH="$GUS_ENV/bin:$PATH"
export XDG_CACHE_HOME="$PROJECT/.cache"
export TMPDIR="/tmp/${USER}/SRS023176_sra_${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"

for cmd in snakemake fasterq-dump python3; do
    command -v "$cmd" >/dev/null || {
        echo "ERROR: $cmd is not available" >&2
        exit 2
    }
done

targets=(
    results/hmp_active_gus/sra_rebuild/SRS023176/SRR059900.fastq_validation.tsv
    results/hmp_active_gus/sra_rebuild/SRS023176/SRR059901.fastq_validation.tsv
)

echo "[$(date -Is)] Downloading and validating SRS023176 SRA runs on $(hostname)"
snakemake -s workflow/hmp_active_gus.smk \
    --configfile config/hmp_active_gus.yaml \
    --cores "${SLURM_CPUS_PER_TASK:-8}" \
    --rerun-incomplete --printshellcmds \
    "${targets[@]}"

tracked_count=$(awk -F '\t' '$1 ~ /^tracked_read_found:/ {sum += $2} END {print sum + 0}' "${targets[@]}")
if [[ "$tracked_count" -ne 1 ]]; then
    echo "ERROR: expected the target instrument read in exactly one run; observed $tracked_count" >&2
    exit 3
fi

echo "[$(date -Is)] Phase 1 complete; target read found exactly once"
awk 'FNR == 1 {print "# " FILENAME} {print}' "${targets[@]}"
echo "[$(date -Is)] Active GUS was not started; inspect these reports before phase 2"
