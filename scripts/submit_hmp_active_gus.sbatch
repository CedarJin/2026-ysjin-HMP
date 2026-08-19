#!/bin/bash
#SBATCH --job-name=hmp-active-gus
#SBATCH --output=logs/hmp_active_gus-%j.out
#SBATCH --error=logs/hmp_active_gus-%j.err
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=low
#SBATCH --account=publicgrp

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

# Run the HMP 139-sample active GUS Snakemake workflow.
#
# From /home/jys0914/2026-ysjin-HMP:
#   sbatch scripts/submit_hmp_active_gus.sbatch
#
# Dry run:
#   sbatch --export=ALL,SNAKEMAKE_EXTRA='-n' scripts/submit_hmp_active_gus.sbatch
#
# Smoke-test-style subset:
#   sbatch --export=ALL,SAMPLE_LIST=config/hmp_active_gus_smoketest_samples.txt,READ_PAIR_LIMIT=100000 scripts/submit_hmp_active_gus.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-/home/jys0914/2026-ysjin-HMP}"
ROOT_DIR="$(pwd)"
mkdir -p logs .cache

SMK_FILE="${ROOT_DIR}/workflow/hmp_active_gus.smk"
CONFIG="${ROOT_DIR}/config/hmp_active_gus.yaml"
GUS_ENV="${GUS_ENV:-/home/jys0914/2026-ysjin-gus/envs/gus-map}"
SAMPLE_LIST="${SAMPLE_LIST:-data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt}"
READ_PAIR_LIMIT="${READ_PAIR_LIMIT:-0}"

export XDG_CACHE_HOME="${ROOT_DIR}/.cache"
export PATH="${GUS_ENV}/bin:${PATH}"

echo "[$(date -Is)] Job started on $(hostname)"
echo "[$(date -Is)] ROOT_DIR=${ROOT_DIR}"
echo "[$(date -Is)] SMK_FILE=${SMK_FILE}"
echo "[$(date -Is)] CONFIG=${CONFIG}"
echo "[$(date -Is)] GUS_ENV=${GUS_ENV}"
echo "[$(date -Is)] SAMPLE_LIST=${SAMPLE_LIST}"
echo "[$(date -Is)] READ_PAIR_LIMIT=${READ_PAIR_LIMIT}"
echo "[$(date -Is)] CORES=${SLURM_CPUS_PER_TASK:-8}"

for f in "${SMK_FILE}" "${CONFIG}" "${SAMPLE_LIST}" \
  "data/HMP_GUS/HMGI/HMGI-stool-gene-seqs.tar.gz" \
  "results/reference/active_gus.faa" \
  "results/reference/active_gus_metadata.tsv"
do
  [[ -f "${f}" ]] || {
    echo "ERROR: missing required file: ${f}" >&2
    exit 1
  }
done

for cmd in snakemake diamond bowtie2-build bowtie2 samtools python3; do
  command -v "${cmd}" >/dev/null || {
    echo "ERROR: ${cmd} not on PATH; check GUS_ENV=${GUS_ENV}" >&2
    exit 1
  }
done

if [[ "${SNAKEMAKE_SKIP_UNLOCK:-0}" != "1" ]]; then
  snakemake --directory "${ROOT_DIR}" -s "${SMK_FILE}" --configfile "${CONFIG}" --unlock
fi

exec snakemake \
  --directory "${ROOT_DIR}" \
  -s "${SMK_FILE}" \
  -c "${SLURM_CPUS_PER_TASK:-8}" \
  --configfile "${CONFIG}" \
  --config "sample_list=${SAMPLE_LIST}" "read_pair_limit=${READ_PAIR_LIMIT}" \
  --rerun-incomplete \
  --printshellcmds \
  ${SNAKEMAKE_EXTRA:-}
