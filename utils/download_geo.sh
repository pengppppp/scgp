#!/usr/bin/env bash
# Download real public data for the JTM nephrology secretome project.
# All files go to D:\Journal of Translational Medicine\data\raw
set -u
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series"
OUT="/d/Journal of Translational Medicine/data/raw"
mkdir -p "$OUT"
LOG="$OUT/download.log"

dl() {  # dl <url> <outfile>
    local url="$1" out="$2"
    echo "[$(date '+%F %T')] START $out" >> "$LOG"
    curl -sSL --retry 5 --retry-delay 10 -C - -o "$out" "$url" >> "$LOG" 2>&1
    local rc=$?
    echo "[$(date '+%F %T')] DONE rc=$rc $out ($(du -h "$out" 2>/dev/null | cut -f1))" >> "$LOG"
    return $rc
}

# filelists first (tiny, tells us archive contents)
dl "$BASE/GSE209nnn/GSE209781/suppl/filelist.txt" "$OUT/GSE209781_filelist.txt"
dl "$BASE/GSE133nnn/GSE133288/suppl/filelist.txt" "$OUT/GSE133288_filelist.txt"
dl "$BASE/GSE180nnn/GSE180395/suppl/filelist.txt" "$OUT/GSE180395_filelist.txt"
dl "$BASE/GSE286nnn/GSE286184/suppl/filelist.txt" "$OUT/GSE286184_filelist.txt"

# series matrices (smaller, quick clinical-annotation inspection)
dl "$BASE/GSE133nnn/GSE133288/matrix/GSE133288_series_matrix.txt.gz" "$OUT/GSE133288_series_matrix.txt.gz"
dl "$BASE/GSE96nnn/GSE96804/matrix/GSE96804_series_matrix.txt.gz" "$OUT/GSE96804_series_matrix.txt.gz"
dl "$BASE/GSE180nnn/GSE180395/matrix/GSE180395_series_matrix.txt.gz" "$OUT/GSE180395_series_matrix.txt.gz"

# big raw tars (run sequentially to be gentle; each resumable)
dl "$BASE/GSE209nnn/GSE209781/suppl/GSE209781_RAW.tar" "$OUT/GSE209781_RAW.tar"
dl "$BASE/GSE180nnn/GSE180395/suppl/GSE180395_RAW.tar" "$OUT/GSE180395_RAW.tar"
dl "$BASE/GSE133nnn/GSE133288/suppl/GSE133288_RAW.tar" "$OUT/GSE133288_RAW.tar"
dl "$BASE/GSE286nnn/GSE286184/suppl/GSE286184_RAW.tar" "$OUT/GSE286184_RAW.tar"

echo "[$(date '+%F %T')] ALL DOWNLOADS FINISHED" >> "$LOG"
