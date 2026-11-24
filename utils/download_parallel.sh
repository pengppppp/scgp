#!/usr/bin/env bash
# Parallel downloads for scRNA datasets (NCBI FTP is per-connection throttled ~145KB/s,
# so we use several parallel connections).
set -u
OUT="/d/Journal of Translational Medicine/data/raw"
LOG="$OUT/download_parallel.log"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series"
mkdir -p "$OUT"

dl() {
    local url="$1" out="$2"
    echo "[$(date '+%F %T')] START $(basename "$out")" >> "$LOG"
    curl -sSL --retry 5 --retry-delay 10 -C - -o "$out" "$url" >> "$LOG" 2>&1
    echo "[$(date '+%F %T')] DONE rc=$? $(basename "$out") ($(du -h "$out" 2>/dev/null | cut -f1))" >> "$LOG"
}

# GSE242425: 12 filtered h5 matrices (hypoxia vs normoxia kidney organoids, d20/d25)
for f in h1_d20 h1_d25 h2_d20 h2_d25 h3_d20 h3_d25 n1_d20 n1_d25 n2_d20 n2_d25 n3_d20 n3_d25; do
    dl "$BASE/GSE242nnn/GSE242425/suppl/GSE242425_${f}_sample_filtered_feature_bc_matrix.h5" \
       "$OUT/GSE242425_${f}_sample_filtered_feature_bc_matrix.h5" &
done
dl "$BASE/GSE242nnn/GSE242425/suppl/GSE242425_barcode_sample.csv.gz" "$OUT/GSE242425_barcode_sample.csv.gz" &

# GSE209781: DKD vs normal human kidney snRNA (single big tar, own connection)
dl "$BASE/GSE209nnn/GSE209781/suppl/GSE209781_RAW.tar" "$OUT/GSE209781_RAW.tar" &

wait
echo "[$(date '+%F %T')] PARALLEL BATCH FINISHED" >> "$LOG"
