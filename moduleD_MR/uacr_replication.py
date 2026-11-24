# -*- coding: utf-8 -*-
"""Independent-outcome replication against UACR (Teumer 2019,
trans-ethnic): protein layer (pQTL instruments, multi-instrument IVW with
Egger/median/Q) and gene layer (cis-eQTL Wald ratios, direction
concordance, panel-background calibration).
Requires data/raw/ckdgen/UACR_overall_ALL.tbl.rsid.gz (see README).

P0-4: UACR (Teumer 2019, n=564,257 trans-ethnic) local-file analysis.

Same design as D9/D10: protein layer (238 LD-clumped deCODE pQTL
instruments, multi-instrument IVW) and gene layer (453 cis-eQTL lead
instruments, Wald ratio), direction concordance with binomial test,
background-rate comparison, BH-FDR within each family.
"""
import gzip
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, "results_tables")
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, norm, chi2

T = r"D:\133\results\tables"
UACR = os.path.join(ROOT, "data", "raw", "ckdgen",
                    "UACR_overall_ALL.tbl.rsid.gz")
OUT_P = os.path.join(T, "D13_UACR_pQTL_MR.csv")
OUT_E = os.path.join(T, "D13_UACR_eQTL_wald.csv")
OUT_S = os.path.join(T, "D13_UACR_summary.csv")

print("读取UACR summary stats...", flush=True)
rows = {}
with gzip.open(UACR, "rt") as f:
    f.readline()
    for line in f:
        p = line.split()
        if len(p) < 10:
            continue
        try:
            rows[p[2]] = {
                "a1_out": p[3],
                "beta_out": float(p[6]),
                "se_out": float(p[7]),
                "p_out": float(p[8]),
                "eaf": float(p[5]),
            }
        except ValueError:
            continue
print(f"  {len(rows)} 个变异", flush=True)


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    running = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(running, 1.0)
    return out


def mr_battery(bx, bxse, by, byse):
    """IVW (multi), Egger, weighted median, Cochran Q."""
    res = {}
    w = 1.0 / byse**2
    beta = np.sum(w * by) / np.sum(w)
    se = np.sqrt(1.0 / np.sum(w))
    res["ivw_beta"] = beta
    res["ivw_p"] = 2 * norm.sf(abs(beta / se))
    # Egger
    if len(bx) >= 3:
        X = np.column_stack([np.ones(len(bx)), bx])
        sw = np.sqrt(w)
        WX, Wy = X * sw[:, None], by * sw
        try:
            coef, *_ = np.linalg.lstsq(WX, Wy, rcond=None)
            resid = Wy - WX @ coef
            sigma2 = np.sum(resid**2) / (len(bx) - 2)
            cov = sigma2 * np.linalg.inv(WX.T @ WX)
            res["egger_p"] = 2 * norm.sf(abs(coef[0] / np.sqrt(cov[0, 0])))
        except Exception:
            res["egger_p"] = np.nan
    else:
        res["egger_p"] = np.nan
    # weighted median
    ratio = by / bx
    wm = 1.0 / (byse / bx)**2
    order = np.argsort(ratio)
    r_s, w_s = ratio[order], wm[order]
    cw = np.cumsum(w_s) - 0.5 * w_s
    cw /= np.sum(w_s)
    idx = np.searchsorted(cw, 0.5)
    b_med = r_s[min(idx, len(r_s) - 1)]
    res["wmedian_p"] = 2 * norm.sf(abs(b_med / np.sqrt(1.0 / np.sum(wm))))
    # Cochran Q
    beta_ref = np.sum(ratio / byse**2) / np.sum(1 / byse**2)
    q = np.sum((ratio - beta_ref)**2 / byse**2)
    res["cochran_q_p"] = chi2.sf(q, len(bx) - 1)
    return res


# ---------- 蛋白层 ----------
print("\n=== 蛋白层 (pQTL × UACR) ===", flush=True)
inst = pd.read_csv(T + r"\D8_pQTL_instruments.csv")
prows = []
for gene, g in inst.groupby("gene"):
    recs = []
    for _, r in g.iterrows():
        o = rows.get(r["snp"])
        if o is None:
            continue
        same = str(r["a1"]).upper() == str(o["a1_out"]).upper()
        recs.append({
            "beta_pqtl": float(r["beta_pqtl"]),
            "se_pqtl": float(r["se_pqtl"]),
            "beta_out": o["beta_out"] if same else -o["beta_out"],
            "se_out": o["se_out"],
            "p_out": o["p_out"],
            "eaf": o["eaf"],
        })
    if len(recs) == 0:
        continue
    d = pd.DataFrame(recs)
    # MAF过滤（eaf 0.01-0.99）
    d = d[(d.eaf > 0.01) & (d.eaf < 0.99)]
    if len(d) == 0:
        continue
    m = mr_battery(d.beta_pqtl.values, d.se_pqtl.values,
                   d.beta_out.values, d.se_out.values)
    prows.append({"gene": gene, "n_matched": len(d),
                  "ivw_beta": m["ivw_beta"], "ivw_p": m["ivw_p"],
                  "egger_p": m["egger_p"], "wmedian_p": m["wmedian_p"],
                  "cochran_q_p": m["cochran_q_p"]})
pdf = pd.DataFrame(prows)
pdf["nominal"] = pdf.ivw_p < 0.05
pdf["bh_q"] = bh_fdr(pdf.ivw_p.values)
pdf["fdr"] = pdf.bh_q < 0.05
pdf.to_csv(OUT_P, index=False)
print(f"  {len(pdf)}个基因测试 | nominal {int(pdf.nominal.sum())} | FDR {int(pdf.fdr.sum())}")
for _, r in pdf[pdf.nominal].sort_values("ivw_p").iterrows():
    print(f"    {r.gene}: p={r.ivw_p:.4g} q={r.bh_q:.3g} (n={r.n_matched})")
for gene in ["NID2", "C2", "TGFB3", "LOXL4"]:
    m = pdf[pdf.gene == gene]
    if len(m):
        r = m.iloc[0]
        print(f"    [{gene}]: p={r.ivw_p:.3f} (n={r.n_matched})")

# ---------- 基因层 ----------

