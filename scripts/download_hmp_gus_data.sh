#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HMGI_DIR="$ROOT_DIR/data/HMP_GUS/HMGI"
HMASM_DIR="$ROOT_DIR/data/HMP_GUS/HMASM/WGS/stool"
MANIFEST_DIR="$ROOT_DIR/data/HMP_GUS/manifests"

mkdir -p "$HMGI_DIR" "$HMASM_DIR" "$MANIFEST_DIR" "$ROOT_DIR/logs"

cat > "$MANIFEST_DIR/HMGI_stool_urls.txt" <<'URLS'
https://human-microbiome-project.s3.amazonaws.com/HHS/HMGI/stool/gene-seqs/HMGI-stool-gene-seqs.tar.gz
https://human-microbiome-project.s3.amazonaws.com/HHS/HMGI/stool/gene-seqs/stool-fastq-seq-length-summary.txt
URLS

echo "[$(date -Is)] Downloading HMGI stool gene/protein sequences"
wget \
  --continue \
  --progress=dot:giga \
  --tries=0 \
  --retry-connrefused \
  --waitretry=30 \
  --timeout=60 \
  --directory-prefix="$HMGI_DIR" \
  --input-file="$MANIFEST_DIR/HMGI_stool_urls.txt"

echo "[$(date -Is)] Downloading HMASM WGS stool sample packages"
wget \
  --continue \
  --progress=dot:giga \
  --tries=0 \
  --retry-connrefused \
  --waitretry=30 \
  --timeout=60 \
  --directory-prefix="$HMASM_DIR" \
  --input-file="$MANIFEST_DIR/HMASM_WGS_stool_urls.txt"

echo "[$(date -Is)] HMP GUS data download complete"
