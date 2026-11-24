# -*- coding: utf-8 -*-
"""Figure 9 (pQTL-MR): LD-clumped panel-wide results.
a: forest of the 15 nominal tests, 6 FDR survivors ringed;
b: NID2 cross-endpoint replication (3 FDR endpoints) - the headline;
c: LD-clumped instrument coverage (238 instruments, 24/25 genes >=3).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results_tables")
OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
C = {"blue": "#4878CF", "red": "#D65F5F", "green": "#6ACC65",
     "orange": "#EE854A", "purple": "#956CB4", "teal": "#82C6E2",
     "gray": "#8C8C8C", "cb_blue": "#0072B2", "cb_orange": "#E69F00"}
EP_SHORT = {"N14_CHRONKIDNEYDIS": "CKD", "N14_RENFAIL": "RenFail",
            "N14_RENALTUB": "RenTub", "DM_NEPHROPATHY": "DM-Neph",
            "N14_IGA_NEPHROPATHY": "IgAN", "I9_HYPTENSRD": "HyptRen"}


def fig9_pqtl():
    d8 = pd.read_csv(os.path.join(TAB, "D8_panel_pQTL_MR.csv"))
    W, H = 7.2, 9.8
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.6, 1.05, 2.2], hspace=0.50)

    # ---- a: forest of nominal tests (15) ----
    axa = fig.add_subplot(gs[0])
    nom = d8[d8["ivw_p"] < 0.05].copy()
    nom["_ep"] = nom["endpoint"].map({e: i for i, e in enumerate(
        ["N14_CHRONKIDNEYDIS", "N14_RENFAIL", "N14_RENALTUB",
         "DM_NEPHROPATHY", "N14_IGA_NEPHROPATHY", "I9_HYPTENSRD"])})
    nom = nom.sort_values(["gene", "_ep"]).reset_index(drop=True)
    z = norm.isf(nom["ivw_p"].clip(lower=1e-300) / 2)
    se = (nom["ivw_beta"] / z).abs()
    y = np.arange(len(nom))[::-1]
    or_l10 = nom["ivw_beta"].values / np.log(10)
    se_l10 = se.values / np.log(10)
    fdr = nom["bh_fdr"].values < 0.05
    axa.errorbar(or_l10, y, xerr=1.96 * se_l10, fmt="none",
                 ecolor="#9A9A9A", elinewidth=1.4, zorder=1)
    axa.scatter(or_l10[~fdr], y[~fdr], s=42, color=C["cb_orange"], zorder=2)
    axa.scatter(or_l10[fdr], y[fdr], s=140, facecolors="none",
                edgecolors=C["red"], linewidths=2.0, zorder=3)
    axa.axvline(0, color="#555555", lw=1.0, ls="--", zorder=0)
    axa.set_yticks(y)
    axa.set_yticklabels([f"{r.gene} - {EP_SHORT[r.endpoint]} "
                         f"(n={int(r.n_instruments)})"
                         for r in nom.itertuples()], fontsize=7.0)
    axa.set_ylim(-0.8, len(nom) - 0.2)
    axa.set_xlabel("IVW effect on kidney outcome (log10 OR per SD of "
                   "genetically predicted protein level)", fontsize=7.5)
    axa.set_title(f"a | panel-wide cis-pQTL MR: 26 proteins x 6 endpoints = "
                  f"150 tests; {len(nom)} nominal, 6 FDR-controlled (ringed)",
                  fontsize=8.6, fontweight="bold", loc="left", pad=9)
    axa.scatter([], [], s=42, color=C["cb_orange"], label="nominal (p<0.05)")
    axa.scatter([], [], s=120, facecolors="none", edgecolors=C["red"],
                label="FDR<0.05")
    axa.legend(fontsize=6.8, loc="lower left", frameon=True,
               framealpha=0.95, edgecolor="#CCCCCC")
    for sp in ["top", "right"]:
        axa.spines[sp].set_visible(False)

    # ---- b: NID2 cross-endpoint replication ----
    axb = fig.add_subplot(gs[1])
    nid2 = d8[d8.gene == "NID2"].copy()
    nid2["_ep"] = nid2["endpoint"].map({e: i for i, e in enumerate(
        ["N14_CHRONKIDNEYDIS", "N14_RENFAIL", "N14_RENALTUB",
         "DM_NEPHROPATHY", "N14_IGA_NEPHROPATHY", "I9_HYPTENSRD"])})
    nid2 = nid2.sort_values("_ep").reset_index(drop=True)
    xb = np.arange(len(nid2))
    vals = -np.log10(nid2["ivw_p"].values)
    colors = [C["red"] if q < 0.05 else C["gray"] for q in nid2["bh_fdr"]]
    axb.bar(xb, vals, color=colors, alpha=0.85, width=0.62)
    axb.axhline(-np.log10(0.05), color="#555555", lw=0.9, ls="--")
    axb.text(5.45, -np.log10(0.05) + 0.15, "p = 0.05", fontsize=6.2,
             color="#555555", ha="right")
    for i, (v, q) in enumerate(zip(vals, nid2["bh_fdr"])):
        axb.text(i, v + 0.12, f"q={q:.4f}" if q < 0.05 else f"q={q:.2f}",
                 fontsize=6.2, ha="center",
                 color=C["red"] if q < 0.05 else "#666666")
    axb.set_xticks(xb)
    axb.set_xticklabels([EP_SHORT[e] for e in nid2["endpoint"]], fontsize=7)
    axb.set_ylabel("-log10 p", fontsize=7.5)
    axb.set_title(f"b | NID2 replicates across THREE kidney endpoints under "
                  f"FDR control (11 LD-clumped cis-pQTL instruments)",
                  fontsize=8.6, fontweight="bold", loc="left", pad=9)
    for sp in ["top", "right"]:
        axb.spines[sp].set_visible(False)

    # ---- c: instrument coverage ----
    axc = fig.add_subplot(gs[2])
    inst = pd.read_csv(os.path.join(TAB, "D8_pQTL_instruments.csv"))
    per_gene = inst.groupby("gene").size().sort_values()
    colors = [C["cb_blue"] if v >= 3 else C["gray"] for v in per_gene.values]
    axc.barh(np.arange(len(per_gene)), per_gene.values, color=colors,
             alpha=0.85, height=0.68)
    axc.set_yticks(np.arange(len(per_gene)))
    axc.set_yticklabels(per_gene.index, fontsize=6.0)
    axc.axvline(3, color=C["red"], lw=1.1, ls="--")
    axc.text(3.3, 0.2, "3-instrument threshold", fontsize=6.4,
             color=C["red"], va="bottom")
    axc.set_xlabel("LD-clumped cis-pQTL instruments (1000G EUR, r2<0.1)",
                   fontsize=7.5)
    axc.set_title(f"c | instrument coverage: {len(inst)} independent "
                  f"instruments; 24 of 25 genes reach >=3 (full MR battery)",
                  fontsize=8.6, fontweight="bold", loc="left", pad=9)
    for sp in ["top", "right"]:
        axc.spines[sp].set_visible(False)

    fig.savefig(os.path.join(OUT, "Fig9_pQTL_MR.png"), dpi=300,
                bbox_inches="tight", pad_inches=0.15)
    plt.close("all")
    print("Fig9 v3 drawn (LD-based results)")


if __name__ == "__main__":
    fig9_pqtl()
