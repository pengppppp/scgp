#!/usr/bin/env bash
# Download eQTLGen cis-eQTL file + FinnGen R12 kidney endpoints for Modules B & D.
set -u
OUT="/d/Journal of Translational Medicine/data/raw"
LOG="$OUT/download_mr.log"
mkdir -p "$OUT"

dl() {
    local url="$1" out="$2"
    # integrity check: gzip files must pass gzip -t; skip only if complete
    if [ -s "$out" ] && gzip -t "$out" 2>/dev/null; then
        echo "[$(date '+%F %T')] SKIP $out (exists and valid)" >> "$LOG"
        return 0
    fi
    if [ -s "$out" ]; then
        echo "[$(date '+%F %T')] REDOWNLOAD $out (existing file incomplete)" >> "$LOG"
        rm -f "$out"
    fi
    echo "[$(date '+%F %T')] START $out" >> "$LOG"
    curl -sSL --retry 5 --retry-delay 10 -C - -o "$out" "$url" >> "$LOG" 2>&1
    local rc=$?
    if [ $rc -eq 0 ] && gzip -t "$out" 2>/dev/null; then
        echo "[$(date '+%F %T')] DONE rc=0 $out ($(du -h "$out" 2>/dev/null | cut -f1))" >> "$LOG"
    else
        echo "[$(date '+%F %T')] INCOMPLETE rc=$rc $out ($(du -h "$out" 2>/dev/null | cut -f1))" >> "$LOG"
    fi
}

dl "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
   "$OUT/eqtlgen_cis_FDR05.txt.gz"

for ep in N14_CHRONKIDNEYDIS N14_RENFAIL N14_RENALTUB DM_NEPHROPATHY N14_IGA_NEPHROPATHY I9_HYPTENSRD; do
    dl "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_${ep}.gz" \
       "$OUT/finngen_R12_${ep}.gz"
done

echo "[$(date '+%F %T')] MR DOWNLOADS FINISHED" >> "$LOG"
