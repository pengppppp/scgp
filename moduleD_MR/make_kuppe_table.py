# -*- coding: utf-8 -*-
"""Cell-type expression of Tier candidates in the Kuppe 2021 kidney
fibrosis atlas (Supplementary Data 2, human/mouse CD10+/CD10-/PDGFRB+).
Requires data/raw/kuppe/MOESM4_gene_expression.xlsx (see README).

Cell-type expression of Tier candidates in the human kidney single-cell
atlas (independent validation of cellular provenance).
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, "results_tables")
import warnings
import pandas as pd

warnings.filterwarnings("ignore")

XLSX = os.path.join(ROOT, "data", "raw", "kuppe",
                     "MOESM4_gene_expression.xlsx")
OUT = os.path.join(T, "D11_kuppe_celltype_expression.csv")

xl = pd.ExcelFile(XLSX)

GENES = ["NID2", "LOXL4", "TUFT1", "C2", "TGFB3", "CASP7", "LGALS3",
         "TNFRSF1A", "CCL2", "LCN2", "ITIH4", "MMP7", "HAVCR1", "UMOD",
         "TNF", "IL6R", "ACE", "DPP4", "HGF", "VEGFA", "KIM1"]

rows = []
for sn in xl.sheet_names:
    df = pd.read_excel(xl, sn)
    if df.shape[1] < 3:
        continue
    gene_col = df.columns[0]
    # ENSG00000227232.5;WASH7P -> WASH7P
    df["_sym"] = (df[gene_col].astype(str).str.split(";").str[-1]
                  .str.split(".").str[0].str.strip().str.upper())
    hit = df[df["_sym"].isin(GENES)]
    is_human = sn.startswith("Human")
    is_level2 = "Level2" in sn
    is_fraction = "Fraction" in sn
    print(f"{sn}: {len(df)} 基因 | 命中 {len(hit)}")
    for _, r in hit.iterrows():
        g = r["_sym"]
        vals = {}
        for c in df.columns[1:]:
            try:
                v = float(r[c])
            except (TypeError, ValueError):
                continue
            if pd.notna(v):
                vals[str(c).strip()] = v
        if not vals:
            continue
        top = sorted(vals.items(), key=lambda x: -x[1])[:3]
        rows.append({
            "dataset": sn,
            "species": "human" if is_human else "mouse",
            "level": 2 if is_level2 else 3,
            "metric": "fraction_detected" if is_fraction else "specificity",
            "gene": g,
            "n_clusters_detected": sum(1 for v in vals.values() if v > 0),
            "top1_celltype": top[0][0], "top1_value": round(top[0][1], 4),
            "top2_celltype": top[1][0] if len(top) > 1 else "",
            "top2_value": round(top[1][1], 4) if len(top) > 1 else None,
            "top3_celltype": top[2][0] if len(top) > 2 else "",
            "top3_value": round(top[2][1], 4) if len(top) > 2 else None,
        })

if rows:
    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)
    # 人类 Level2/3 Fraction 汇总（主证据）
    hum = out[(out.species == "human") & (out.metric == "fraction_detected")]
    print(f"\n=== 人类数据集（{len(hum)}行）===")
    for ds in sorted(hum.dataset.unique()):
        sub = hum[hum.dataset == ds]
        print(f"\n--- {ds} ---")
        for _, r in sub.iterrows():
            print(f"  {r.gene:10s} top1: {r.top1_celltype:42s} {r.top1_value:.3f} "
                  f"| top2: {r.top2_celltype[:30]:30s} {r.top2_value:.3f}")
else:
    print("无命中")
