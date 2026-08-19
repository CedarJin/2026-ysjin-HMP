#!/bin/bash
# Build shared DIAMOND databases once before submitting the sample array.
set -euo pipefail

PROJECT="${PROJECT:-/home/jys0914/2026-ysjin-HMP}"
GUS_ROOT="${GUS_ROOT:-/home/jys0914/2026-ysjin-gus}"
GUS_ENV="${GUS_ENV:-${GUS_ROOT}/envs/gus-map}"
cd "$PROJECT"

export PATH="${GUS_ENV}/bin:${PATH}"
command -v diamond >/dev/null || {
    echo "ERROR: diamond not on PATH; check GUS_ENV=${GUS_ENV}" >&2
    exit 1
}

REF_OUT=results/hmp_assembly_gus/reference
REF_IN="${GUS_ROOT}/results/reference"
mkdir -p "$REF_OUT" logs

makedb_if_needed() {
    local faa="$1"
    local prefix="$2"
    if [[ -s "${prefix}.dmnd" ]]; then
        echo "[skip] ${prefix}.dmnd"
        return 0
    fi
    echo "[makedb] ${faa} -> ${prefix}"
    diamond makedb --in "$faa" --db "$prefix"
}

makedb_if_needed "${REF_IN}/gus_279_plus_fmn.faa" "${REF_OUT}/gus_279_plus_fmn"
makedb_if_needed "${REF_IN}/candidate_gus.faa" "${REF_OUT}/candidate_gus"
makedb_if_needed "${REF_IN}/active_gus.faa" "${REF_OUT}/active_gus"

echo "[$(date -Is)] Reference databases ready under ${REF_OUT}"
