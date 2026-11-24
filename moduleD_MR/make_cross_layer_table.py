# -*- coding: utf-8 -*-
"""Cross-layer evidence table (T0): one row per nominally replicated gene,
joining the eQTL-MR (D2), colocalisation (D5), pQTL-MR (D6/D8),
CKDGen replication (D9) and annotation (B2, F1) layers.

Usage: python make_cross_layer_table.py
Writes: results_tables/T0_cross_layer_evidence.csv
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, "results_tables")

EP_ABBR = {
    "CHRONKIDNEYDIS": "CKD",
    "RENFAIL": "renal failure",
    "RENALTUB": "tubulo-interstitial",
    "N14_NEPHROPAT": "nephropathy",
    "IGA_NEPHROPATHY": "IgA nephropathy",
    "HYPTENSRD": "hypertensive renal",
}


def load(name):
    p = os.path.join(T, name)
    if not os.path.exists(p):
        print(f"missing {name}; run the module scripts first", file=sys.stderr)
        sys.exit(1)
    return pd.read_csv(p)


def main():
    d2 = load("D2_MR_results_kidney.csv")
    d5 = load("D5_coloc_results.csv")
    d6 = load("D6_pQTL_MR_replication.csv")
    d8 = load("D8_panel_pQTL_MR.csv")
    d9 = load("D9_D10_independent_outcome_summary.csv")
    b2 = load("B2_panel_plasma_detectability.csv")
    f1 = load("F1_drug_mapping.csv")

    rows = []
    for gene, g in d2.groupby("gene"):
        g = g.copy()
        g["nominal"] = g.ivw_p < 0.05 if "ivw_p" in g else g.p < 0.05
        pcol = "ivw_p" if "ivw_p" in g else "p"
        nom = g[g[pcol] < 0.05]
        if len(nom) < 2:
            continue  # T0 covers genes nominally replicated in >= 2 endpoints

        eps = [EP_ABBR.get(e, e) for e in nom.endpoint]
        ors = []
        if "OR" in nom:
            ors = [float(x) for x in nom.OR if pd.notna(x)]
        direction = ("risk" if np.nanmedian(ors) > 1 else "protective"
                     if ors else "-")

        # FDR 端点（D2 或 D5 无 fdr 列时用 D3 的结论：此处按 D2 p 与
        # focused-panel 阈值复算不可行，直接采用 D3 表）
        try:
            d3 = load("D3_replication_enrichment.csv")
        except SystemExit:
            d3 = None

        # pQTL 层
        pq = d6[d6.gene == gene] if "gene" in d6 else pd.DataFrame()
        pq8 = d8[d8.gene == gene] if "gene" in d8 else pd.DataFrame()
        if len(pq8) and "fdr" in pq8:
            n_fdr_pq = int(pq8.fdr.sum())
        elif len(pq):
            n_fdr_pq = 0
        else:
            n_fdr_pq = "not covered"
        if len(pq8) and "ivw_p" in pq8 and pq8.ivw_p.notna().any():
            best = pq8.loc[pq8.ivw_p.idxmin()]
            pq_strong = f"min p={best.ivw_p:.3g}, {best.endpoint}"
        else:
            pq_strong = "-"

        # coloc
        c5 = d5[d5.gene == gene] if "gene" in d5 else pd.DataFrame()
        pp4 = (float(c5.PP_H4.max()) if len(c5) and "PP_H4" in c5
               else (float(c5.PP.H4.max()) if len(c5) and "PP.H4" in c5
                     else np.nan))

        # CKDGen eGFRcrea (D9/D10 summary 无逐基因 p；用 D9 trans 表)
        try:
            d9g = load("D9_eGFRcrea_trans_eQTL_wald.csv")
            w = d9g[d9g.gene == gene]
            ck = f"{w.wald_p.min():.3g}" if len(w) else "not tested"
        except SystemExit:
            ck = "not tested"

        # 注释层
        soma = ("yes" if len(b2[b2.gene == gene]) else "no")
        drugs = f1[f1.gene == gene] if "gene" in f1 else pd.DataFrame()
        n_drug = ("-" if not len(drugs) else
                  str(int((drugs.get("n_approved", pd.Series([0]))).sum())))

        rows.append({
            "Gene": gene,
            "Tier": "A",  # tier 复算见 D3；此处只标名义复现基因
            "Nominal endpoints (n/6)": len(nom),
            "Endpoints (abbreviated)": ", ".join(eps),
            "Direction": direction,
            "OR range": (f"{min(ors):.2f}-{max(ors):.2f}" if ors else "-"),
            "FDR endpoints (n/6)": "-",  # 见 D3_replication_enrichment.csv
            "pQTL-MR FDR endpoints (n/6)": n_fdr_pq,
            "pQTL-MR strongest": pq_strong,
            "Colocalisation max PP.H4": round(pp4, 2) if pp4 == pp4 else "-",
            "CKDGen eGFRcrea nominal p": ck,
            "SomaScan plasma-detectable": soma,
            "Approved drugs": n_drug,
        })

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(T, "T0_cross_layer_evidence.csv"), index=False)
    print(f"wrote T0_cross_layer_evidence.csv ({len(out)} genes)")


if __name__ == "__main__":
    main()
