#!/usr/bin/env bash
# SCGP: Single-Cell anchored Genetic Prioritisation
# Companion code for "From kidney single-cell states to blood-accessible
# secreted candidates: a reproducible prioritisation resource for chronic
# kidney disease" (Journal of Translational Medicine submission).
# Public data only. Each step logs to results/logs/.
set -euo pipefail
# Path-independent: derive the package root from this script's location.
# After unzipping the archive, run:  bash scgp.sh
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python3}"   # set PYTHON=/path/to/python to override

echo "[1/13] Module A: single-cell cell states (GSE242425 organoids)"
$PY moduleA_scRNA/01_scRNA_cellstates.py --dataset GSE242425
echo "[2/13] Module A: single-cell cell states (GSE209781 DKD tissue)"
$PY moduleA_scRNA/01_scRNA_cellstates.py --dataset GSE209781
echo "[3/13] Module B: state-restricted secretome"
$PY moduleB_secretome/02_secretome.py
echo "[4/13] Module C: bulk cohort validation"
$PY moduleC_bulk/03_bulk_cohorts.py
echo "[5/13] Module D: cis-eQTL MR (six kidney endpoints, lead-SNP primary)"
$PY moduleD_MR/04_MR_causal.py
echo "[6/13] Module D: endpoint genetic correlation (D4)"
$PY moduleD_MR/04b_endpoint_correlation.py
echo "[7/13] Module D: colocalisation (D5)"
$PY moduleD_MR/04c_coloc.py
echo "[8/13] Module D: focused multi-instrument cis-pQTL MR (D6)"
$PY moduleD_MR/04d_pQTL_MR.py
echo "[9/13] Module E: cross-endpoint replication under FDR control (D3)"
$PY moduleE_population/05_population.py
echo "[10/13] Module F: druggability mapping (F1-F3)"
$PY moduleF_drug/06_drug_mapping.py
echo "[11/13] Module D: panel-wide multi-instrument cis-pQTL MR (D8) +"
echo "        cross-layer evidence summary (T0)"
$PY moduleD_MR/04f_panel_pQTL_MR.py
$PY moduleD_MR/make_cross_layer_table.py
echo "[12/13] Module D: independent-outcome replication (CKDGen, Biobank Japan)"
$PY moduleD_MR/04e_independent_replication.py
echo "[13/14] Figures: module scripts write their own; Fig 9 here"
$PY utils/make_fig9.py
echo "[14/14] Consolidate result tables into S1-S7 Excel workbooks (GigaDB)"
$PY utils/consolidate_tables.py
echo "ALL DONE. Tables in results_tables/, GigaDB workbooks in 03_tables_for_GigaDB/, logs in results/logs/."
echo "Supplementary analyses (need large downloads into data/raw/, see README):"
echo "  make_kuppe_table.py (D11), nid2_ckd_binary.py (D12), uacr_replication.py (D13)"
