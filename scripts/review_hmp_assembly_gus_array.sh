#!/bin/bash
# Concise review for the HMP assembly GUS array — no need to open 139 logs.
#
# Usage:
#   bash scripts/review_hmp_assembly_gus_array.sh
#   bash scripts/review_hmp_assembly_gus_array.sh --failed-tail 40
#   bash scripts/review_hmp_assembly_gus_array.sh --job 19901234

set -euo pipefail

PROJECT="${PROJECT:-/home/jys0914/2026-ysjin-HMP}"
SAMPLE_LIST="${SAMPLE_LIST:-${PROJECT}/data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt}"
STATUS_TSV="${PROJECT}/logs/hmp_assembly_gus_status.tsv"
RESULTS="${PROJECT}/results/hmp_assembly_gus"
FAILED_TAIL=20
JOB_FILTER=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --failed-tail)
            FAILED_TAIL="${2:?}"
            shift 2
            ;;
        --job)
            JOB_FILTER="${2:?}"
            shift 2
            ;;
        -h|--help)
            sed -n '2,10p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

cd "$PROJECT"
mapfile -t ALL_SAMPLES < <(grep -v '^[[:space:]]*$' "${SAMPLE_LIST}")
N="${#ALL_SAMPLES[@]}"

done_out() {
    local s="$1"
    [[ -s "${RESULTS}/summary/pident95_cov90/abundance/${s}.min1000.active_gus_abundance.tsv" ]] \
        && [[ -s "${RESULTS}/summary/${s}.min1000.total_gus_abundance.tsv" ]] \
        && [[ -s "${RESULTS}/summary/${s}.min1000.candidate_gus_abundance.tsv" ]]
}

have_status=0
[[ -s "${STATUS_TSV}" ]] && have_status=1

declare -A LAST_STATUS
declare -A LAST_EXIT
declare -A LAST_JOB
declare -A LAST_TASK

if [[ "${have_status}" -eq 1 ]]; then
    # Robust to older status files that accidentally omitted newlines between rows.
    while IFS=$'\t' read -r ts job task sample status exit_code host; do
        [[ "${ts}" == "timestamp" || -z "${ts}" ]] && continue
        [[ -n "${JOB_FILTER}" && "${job}" != "${JOB_FILTER}" ]] && continue
        [[ -z "${sample}" || -z "${status}" ]] && continue
        LAST_STATUS["${sample}"]="${status}"
        LAST_EXIT["${sample}"]="${exit_code}"
        LAST_JOB["${sample}"]="${job}"
        LAST_TASK["${sample}"]="${task}"
    done < <(STATUS_TSV_PATH="${STATUS_TSV}" python3 - <<'PY'
from pathlib import Path
import os
import re

text = Path(os.environ["STATUS_TSV_PATH"]).read_text()
parts = re.split(r"(?=20\d{2}-\d{2}-\d{2}T)", text)
print("timestamp\tarray_job_id\ttask_id\tsample\tstatus\texit_code\thost")
for part in parts:
    part = part.strip()
    if not part or part.startswith("timestamp"):
        continue
    cols = part.split("\t")
    if len(cols) < 6:
        continue
    print("\t".join(cols[:7]))
PY
)
fi

n_done=0
n_failed=0
n_running_or_unknown=0
n_missing_out=0
FAILED_SAMPLES=()
MISSING_OUT=()

for s in "${ALL_SAMPLES[@]}"; do
    if done_out "$s"; then
        n_done=$((n_done + 1))
    else
        n_missing_out=$((n_missing_out + 1))
        MISSING_OUT+=("$s")
    fi

    st="${LAST_STATUS[$s]:-}"
    if [[ "${st}" == "FAILED" ]]; then
        n_failed=$((n_failed + 1))
        FAILED_SAMPLES+=("$s")
    elif [[ "${st}" == "DONE" ]]; then
        :
    else
        n_running_or_unknown=$((n_running_or_unknown + 1))
    fi
done

echo "=== HMP assembly GUS array review ==="
echo "samples total:           ${N}"
echo "final outputs present:   ${n_done}"
echo "final outputs missing:   ${n_missing_out}"
if [[ "${have_status}" -eq 1 ]]; then
    echo "status DONE (last):      $(grep -c $'\tDONE\t' "${STATUS_TSV}" 2>/dev/null || true) rows in ${STATUS_TSV}"
    echo "status FAILED (unique):  ${n_failed}"
    echo "no status yet/unknown:   ${n_running_or_unknown}"
else
    echo "status file:             missing (${STATUS_TSV})"
fi
echo

# Live Slurm view if available
if command -v squeue >/dev/null 2>&1; then
    echo "=== currently running / pending (squeue) ==="
    if [[ -n "${JOB_FILTER}" ]]; then
        squeue -j "${JOB_FILTER}" -o '%.18i %.9P %.10T %.10M %.6D %R' 2>/dev/null || true
    else
        squeue -u "${USER}" -n hmp_asm_gus -o '%.18i %.9P %.10T %.10M %.6D %R' 2>/dev/null || true
    fi
    echo
fi

if [[ "${#FAILED_SAMPLES[@]}" -gt 0 ]]; then
    echo "=== FAILED samples (from status TSV) ==="
    printf '%-12s %-12s %-6s %-8s %s\n' SAMPLE JOB TASK EXIT LOG
    for s in "${FAILED_SAMPLES[@]}"; do
        job="${LAST_JOB[$s]}"
        task="${LAST_TASK[$s]}"
        slog="logs/hmp_assembly_gus_samples/${s}.log"
        elogg="logs/hmp_assembly_gus_${job}_${task}.err"
        printf '%-12s %-12s %-6s %-8s %s\n' "$s" "$job" "$task" "${LAST_EXIT[$s]}" "$slog"
        if [[ "${FAILED_TAIL}" -gt 0 && -f "${slog}" ]]; then
            echo "  --- tail ${slog} ---"
            tail -n "${FAILED_TAIL}" "${slog}" | sed 's/^/  /'
        elif [[ "${FAILED_TAIL}" -gt 0 && -f "${elogg}" ]]; then
            echo "  --- tail ${elogg} ---"
            tail -n "${FAILED_TAIL}" "${elogg}" | sed 's/^/  /'
        fi
    done
    echo
fi

if [[ "${n_missing_out}" -gt 0 && "${n_missing_out}" -le 30 ]]; then
    echo "=== samples missing final outputs ==="
    printf '%s\n' "${MISSING_OUT[@]}"
elif [[ "${n_missing_out}" -gt 30 ]]; then
    echo "=== samples missing final outputs (first 30 / ${n_missing_out}) ==="
    printf '%s\n' "${MISSING_OUT[@]:0:30}"
    echo "... (${n_missing_out} total; see status TSV / rerun review later)"
fi

echo
echo "Tips:"
echo "  status table:  ${STATUS_TSV}"
echo "  per-sample:    logs/hmp_assembly_gus_samples/<SAMPLE>.log"
echo "  only failures: bash scripts/review_hmp_assembly_gus_array.sh --failed-tail 40"
echo "  resubmit one:  sbatch --array=<task_id> scripts/submit_hmp_assembly_gus_array.sbatch"
