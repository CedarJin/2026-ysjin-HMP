#!/bin/bash
#SBATCH --job-name=hmp_asm_gus
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-12:00:00
#SBATCH --array=1-139%10
#SBATCH --output=logs/hmp_assembly_gus_%A_%a.out
#SBATCH --error=logs/hmp_assembly_gus_%A_%a.err

#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=FAIL


# Per-sample MEGAHIT + GUS abundance array (mirrors 2026-ysjin-gus).
# 139 tasks, at most 10 running at once (%10); each task = one sample.
#
# Setup once:
#   bash scripts/setup_hmp_assembly_gus_refs.sh
#
# Submit all:
#   sbatch scripts/submit_hmp_assembly_gus_array.sbatch
#
# Review progress (no need to open 139 logs):
#   bash scripts/review_hmp_assembly_gus_array.sh
#
# Pilot first 5 samples:
#   sbatch --array=1-5%5 scripts/submit_hmp_assembly_gus_array.sbatch

set -euo pipefail

PROJECT="${PROJECT:-/home/jys0914/2026-ysjin-HMP}"
GUS_ROOT="${GUS_ROOT:-/home/jys0914/2026-ysjin-gus}"
GUS_ENV="${GUS_ENV:-${GUS_ROOT}/envs/gus-map}"
ASSEMBLE_ENV="${ASSEMBLE_ENV:-/home/jys0914/.conda/envs/assemble}"
SAMPLE_LIST="${SAMPLE_LIST:-${PROJECT}/data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt}"
SMK_FILE="${PROJECT}/workflow/hmp_assembly_gus.smk"
CONFIG="${PROJECT}/config/hmp_assembly_gus.yaml"
STATUS_TSV="${PROJECT}/logs/hmp_assembly_gus_status.tsv"

cd "$PROJECT"
mkdir -p logs/hmp_assembly_gus_samples \
    results/hmp_assembly_gus/{reads,assembly/megahit,mapping,prodigal,counts,search,summary,reference,logs/megahit}
mkdir -p results/hmp_assembly_gus/summary/pident95_cov90/{abundance,hits}

export XDG_CACHE_HOME="${PROJECT}/.cache"
export TMPDIR="${TMPDIR:-/tmp/${USER}/hmp_asm_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}}"
mkdir -p "$XDG_CACHE_HOME" "$TMPDIR"

# megahit from assemble env; everything else from gus-map
export PATH="${ASSEMBLE_ENV}/bin:${GUS_ENV}/bin:${PATH}"

if [[ -n "${SAMPLE:-}" ]]; then
    :
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")"
else
    echo "ERROR: set SAMPLE=... or submit as a Slurm array job" >&2
    exit 2
fi

[[ -n "${SAMPLE}" ]] || {
    echo "ERROR: empty SAMPLE for array task ${SLURM_ARRAY_TASK_ID:-NA}" >&2
    exit 2
}

ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-NA}"
SAMPLE_LOG="logs/hmp_assembly_gus_samples/${SAMPLE}.log"

record_status() {
    local status="$1"
    local exit_code="$2"
    # Write the newline directly; do not capture via $(...) which strips trailing newlines.
    mkdir -p logs
    (
        flock 9
        if [[ ! -s "${STATUS_TSV}" ]]; then
            printf 'timestamp\tarray_job_id\ttask_id\tsample\tstatus\texit_code\thost\n' > "${STATUS_TSV}"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$(date -Is)" "${ARRAY_JOB_ID}" "${TASK_ID}" "${SAMPLE}" "${status}" "${exit_code}" "$(hostname)" \
            >> "${STATUS_TSV}"
    ) 9>>"${STATUS_TSV}.lock"
}

on_exit() {
    local rc=$?
    if [[ "${rc}" -eq 0 ]]; then
        record_status "DONE" "0"
        echo "[$(date -Is)] DONE sample=${SAMPLE}"
    else
        record_status "FAILED" "${rc}"
        echo "[$(date -Is)] FAILED sample=${SAMPLE} exit=${rc}" >&2
        echo "[$(date -Is)] see ${SAMPLE_LOG} and logs/hmp_assembly_gus_${ARRAY_JOB_ID}_${TASK_ID}.err" >&2
    fi
    exit "${rc}"
}
trap on_exit EXIT

echo "[$(date -Is)] host=$(hostname) job=${SLURM_JOB_ID:-local} task=${TASK_ID} sample=${SAMPLE}"
echo "[$(date -Is)] detailed log: ${SAMPLE_LOG}"
echo "[$(date -Is)] cores=${SLURM_CPUS_PER_TASK:-16} mem=${SLURM_MEM_PER_NODE:-NA}"

for cmd in megahit snakemake prodigal diamond bowtie2 bowtie2-build samtools featureCounts python3; do
    command -v "$cmd" >/dev/null || {
        echo "ERROR: ${cmd} not on PATH" >&2
        exit 1
    }
done

# Shared DIAMOND DBs (flock so concurrent array tasks do not race)
LOCK="${PROJECT}/results/hmp_assembly_gus/reference/.makedb.lock"
mkdir -p "$(dirname "$LOCK")"
(
    flock 9
    bash scripts/setup_hmp_assembly_gus_refs.sh
) 9>"$LOCK"

# Prefer repaired reads for SRS023176; otherwise HMASM trimmed reads.
R1="results/hmp_active_gus/reads/${SAMPLE}/${SAMPLE}.denovo_duplicates_marked.trimmed.1.fastq"
R2="results/hmp_active_gus/reads/${SAMPLE}/${SAMPLE}.denovo_duplicates_marked.trimmed.2.fastq"
if [[ "${SAMPLE}" == "SRS023176" ]]; then
    R1="results/hmp_active_gus/repaired_reads/${SAMPLE}/${SAMPLE}.denovo_duplicates_marked.trimmed.repaired.1.fastq"
    R2="results/hmp_active_gus/repaired_reads/${SAMPLE}/${SAMPLE}.denovo_duplicates_marked.trimmed.repaired.2.fastq"
fi
for f in "$R1" "$R2" "$SMK_FILE" "$CONFIG"; do
    [[ -e "$f" ]] || {
        echo "ERROR: missing required file: $f" >&2
        exit 1
    }
done

THREAD_CFG="${THREADS_PER_TOOL:-8}"

# Keep Slurm .out/.err short; put full snakemake transcript in per-sample log.
snakemake \
    -s "$SMK_FILE" \
    --configfile "$CONFIG" \
    --config "sample=${SAMPLE}" \
        "threads=${THREAD_CFG}" \
        "threads_megahit=${THREAD_CFG}" \
        "threads_bowtie2=${THREAD_CFG}" \
        "threads_featurecounts=${THREAD_CFG}" \
        "threads_prodigal=${THREAD_CFG}" \
        "threads_diamond=${THREAD_CFG}" \
    --cores "${SLURM_CPUS_PER_TASK:-16}" \
    --rerun-incomplete \
    --nolock \
    --printshellcmds \
    ${SNAKEMAKE_EXTRA:-} \
    >"${SAMPLE_LOG}" 2>&1
