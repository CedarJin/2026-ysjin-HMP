#!/bin/bash
#SBATCH --job-name=SRS023176_hmp_gus
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/SRS023176_hmp_gus.%j.out
#SBATCH --error=logs/SRS023176_hmp_gus.%j.err

#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

# Phase 2: submit only after both phase-1 FASTQ validation reports are reviewed.

set -euo pipefail

PROJECT=/home/jys0914/2026-ysjin-HMP
GUS_ENV=/home/jys0914/2026-ysjin-gus/envs/gus-map
cd "$PROJECT"

mkdir -p logs .cache
export PATH="$GUS_ENV/bin:$PATH"
export XDG_CACHE_HOME="$PROJECT/.cache"

reports=(
    results/hmp_active_gus/sra_rebuild/SRS023176/SRR059900.fastq_validation.tsv
    results/hmp_active_gus/sra_rebuild/SRS023176/SRR059901.fastq_validation.tsv
)
for report in "${reports[@]}"; do
    [[ -s "$report" ]] || {
        echo "ERROR: missing phase-1 validation report: $report" >&2
        exit 2
    }
    awk -F '\t' '$1 == "error_pairs" && $2 != 0 {exit 1}
                    $1 == "pair_count_matches_expected" && $2 != 1 {exit 1}' "$report" || {
        echo "ERROR: failed phase-1 validation report: $report" >&2
        exit 3
    }
done

tracked_count=$(awk -F '\t' '$1 ~ /^tracked_read_found:/ {sum += $2} END {print sum + 0}' "${reports[@]}")
[[ "$tracked_count" -eq 1 ]] || {
    echo "ERROR: target instrument read was not found exactly once" >&2
    exit 4
}

targets=(
    results/hmp_active_gus/sra_rebuild/SRS023176/SRS023176.processed_fastq_validation.tsv
    results/hmp_active_gus/summary/SRS023176.active_gus_abundance.tsv
    results/hmp_active_gus/summary/active_gus_abundance_all_samples.tsv
    results/hmp_active_gus/summary/active_gus_tpm_genome_wide_matrix.tsv
    results/hmp_active_gus/summary/active_gus_tpm_within_active_db_matrix.tsv
)

echo "[$(date -Is)] Starting HMP preprocessing and active GUS on $(hostname)"
snakemake -s workflow/hmp_active_gus.smk \
    --configfile config/hmp_active_gus.yaml \
    --config threads=32 threads_bowtie2=32 threads_featurecounts=32 threads_diamond=32 \
    --cores "${SLURM_CPUS_PER_TASK:-32}" \
    --rerun-incomplete --printshellcmds \
    "${targets[@]}"
echo "[$(date -Is)] Phase 2 complete"
