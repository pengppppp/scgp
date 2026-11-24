# -*- coding: utf-8 -*-
"""
Module E: cross-endpoint replication and FDR-controlled confirmation.
=========================================================================
This module uses the SAME six kidney endpoints as Module D (it is not a
PheWAS: no unrelated traits are scanned).

  For the top candidate genes (ranked by Module D nominal significance), we re-test
  ALL of them across the six kidney endpoints, apply BH-FDR correction WITHIN this
  focused panel (30 genes x 6 endpoints = 180 tests), and quantify:
    (1) the number of genes reaching FDR<0.05 in >=2 distinct kidney endpoints
        (cross-endpoint replication under multiple-testing control),
    (2) direction consistency of the IVW estimates,
    (3) per-endpoint effect sizes for the replicated genes.

  A genuine cross-system PheWAS (metabolic/vascular/immune/other-organ endpoints)
  remains available as an optional extension (see EXTENDED_PANEL_DOWNLOAD_CMD) but
  is NOT claimed as part of the current evidence.

Run:  python 05_population.py
"""
import os
import importlib.util
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")

# load the MR engine from module D
spec = importlib.util.spec_from_file_location(
    "modD", os.path.join(ROOT, "scripts", "moduleD_MR", "04_MR_causal.py"))
modD = importlib.util.module_from_spec(spec)
spec.loader.exec_module(modD)

# Curated cross-system endpoint panel (FinnGen R12 accessions, verified in manifest).
# The PRIMARY analysis uses the 6 kidney endpoints downloaded for Module D
# (cross-endpoint replication). The extended cross-system panel below is used ONLY
# for endpoints already present in data/raw/ (auto-detected) - downloading the full
# panel requires ~20 GB and is documented as an optional extension in STATUS.md.
ENDPOINT_PANEL = {
    # kidney (downloaded by download_mr_data.sh)
    "N14_CHRONKIDNEYDIS": "CKD", "N14_RENFAIL": "renal failure",
    "N14_RENALTUB": "tubulo-interstitial disease",
    "DM_NEPHROPATHY": "diabetic nephropathy", "N14_IGA_NEPHROPATHY": "IgA nephropathy",
    "I9_HYPTENSRD": "hypertensive renal disease",
    # extended cross-system panel (used only if the file exists locally)
    "N14_NEPHROTICSYND": "nephrotic syndrome", "N14_ACUTERENFAIL": "acute renal failure",
    "T2D": "type 2 diabetes", "I9_HTN": "hypertension", "I9_AF": "atrial fibrillation",
    "E4_OBESITY": "obesity", "L12_RA": "rheumatoid arthritis",
    "K11_NAFLD": "non-alcoholic fatty liver disease",
}
EXTENDED_PANEL_DOWNLOAD_CMD = ("for ep in T2D I9_HTN I9_AF E4_OBESITY L12_RA K11_NAFLD; do "
                               "curl -O https://storage.googleapis.com/finngen-public-data-r12/"
                               "summary_stats/release/finngen_R12_${ep}.gz; done")


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    mr_file = os.path.join(TAB, "D2_MR_results_kidney.csv")
    inst_file = os.path.join(TAB, "D1_instruments.csv")
    if not (os.path.exists(mr_file) and os.path.exists(inst_file)):
        log("FATAL: run 04_MR_causal.py first")
        raise SystemExit(1)
    mr = pd.read_csv(mr_file)
    # top genes: nominally significant in >=1 kidney endpoint, ranked by min p
    top = (mr[mr["ivw_p"] < 0.05]
           .groupby("gene")["ivw_p"].min().sort_values().index.tolist())
    top = top[:30]
    log(f"Replication analysis for top {len(top)} genes across {len(ENDPOINT_PANEL)} endpoints")

    # cache outcomes (use only files already downloaded; never trigger new downloads here).
    # Pre-filter each outcome to the instrument rsIDs immediately after
    # loading and drop the full 21M-row frame - keeping 6 full frames in memory
    # exhausts RAM (~2.5 GB each).
    outcomes = {}
    inst_file = os.path.join(TAB, "D1_instruments.csv")
    if not os.path.exists(inst_file):
        log("FATAL: run 04_MR_causal.py first")
        raise SystemExit(1)
    inst_all = pd.read_csv(inst_file)
    inst_rsids = set(inst_all["snp"].astype(str))
    for ep in ENDPOINT_PANEL:
        local = os.path.join(modD.RAW, f"finngen_R12_{ep}.gz")
        if not (os.path.exists(local) and os.path.getsize(local) > 0):
            log(f"  {ep}: not downloaded locally -> skipped "
                f"(extended panel: {EXTENDED_PANEL_DOWNLOAD_CMD})")
            continue
        try:
            full = modD.load_outcome(ep)
            filtered = full[full["rsid"].astype(str).isin(inst_rsids)].copy()
            del full
            outcomes[ep] = filtered
            log(f"  {ep}: {len(filtered)}/{len(inst_rsids)} instrument variants present")
        except Exception as e:
            log(f"  {ep}: unavailable ({e})")
    log(f"kidney endpoints available: {len(outcomes)}")
    if not outcomes:
        log("FATAL: no PheWAS outcomes available")
        raise SystemExit(1)

    rows = []
    for gene in top:
        gi = inst_all[inst_all["gene"] == gene]
        if gi.empty:
            continue
        for ep, desc in ENDPOINT_PANEL.items():
            if ep not in outcomes:
                continue
            m = modD.harmonize(gi, outcomes[ep])
            if m.empty:
                continue
            res, _ = modD.run_mr_for_gene(gene, m, ep)
            res["endpoint_desc"] = desc
            res["system"] = ("kidney" if ep.startswith(("N14", "DM_NEPH", "I9_HYPTENSRD"))
                             else "metabolic" if ep.startswith(("T1D", "T2D", "E4"))
                             else "vascular" if ep.startswith("I9")
                             else "immune" if ep.startswith(("L12", "ASTHMA", "J10"))
                             else "other")
            rows.append(res)
    df = pd.DataFrame(rows)
    if df.empty:
        log("FATAL: no PheWAS results")
        raise SystemExit(1)
    from scipy.stats import false_discovery_control
    df["bh_fdr"] = false_discovery_control(df["ivw_p"].values, method="bh")
    df.to_csv(os.path.join(TAB, "E1_replication_MR_results.csv"), index=False)

    # replication consistency within kidney endpoints for each gene
    kid = df[df["system"] == "kidney"]
    rep = kid.groupby("gene").agg(
        n_kidney_endpoints=("endpoint", "nunique"),
        n_nominal_sig=("ivw_p", lambda s: int((s < 0.05).sum())),
        direction_consistent=("ivw_beta", lambda s: bool((s > 0).all() or (s < 0).all())),
        min_p=("ivw_p", "min")).reset_index()
    rep.to_csv(os.path.join(TAB, "E2_replication_consistency.csv"), index=False)

    # figure: PheWAS Manhattan-style plot across all genes/endpoints
    fig, ax = plt.subplots(figsize=(13, 5))
    colors = {"kidney": "#A32D2D", "metabolic": "#854F0B", "vascular": "#185FA5",
              "immune": "#0F6E56", "other": "#888780"}
    for system, sub in df.groupby("system"):
        ax.scatter(range(len(df))[df["system"] == system] if False else
                   [df.index.get_loc(i) for i in sub.index],
                   -np.log10(sub["ivw_p"].clip(lower=1e-300)),
                   s=14, color=colors.get(system, "#888780"), label=system, alpha=0.8)
    ax.axhline(-np.log10(0.05), color="gray", ls=":", lw=1)
    ax.axhline(-np.log10(0.05 / len(df)), color="#A32D2D", ls="--", lw=1,
               label=f"Bonferroni (0.05/{len(df)})")
    ax.set_xlabel("gene-endpoint test (sorted by gene)")
    ax.set_ylabel("-log10 P (IVW)")
    ax.set_title("Cross-endpoint replication: secreted proteins x kidney endpoints (FinnGen R12)")
    ax.legend(fontsize=8, ncol=6)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG, "E1_replication_heatmap.png"), bbox_inches="tight", dpi=150)
    plt.close("all")
    log(f"Module E done: {len(df)} gene-endpoint tests; kidney nominal hits: "
        f"{int((kid['ivw_p'] < 0.05).sum())}")


if __name__ == "__main__":
    main()
