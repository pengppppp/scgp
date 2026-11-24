# -*- coding: utf-8 -*-
"""Construct-matched test of the NID2 pQTL instruments against the
CKDGen CKD binary outcome (Wuttke 2019; dbgap-format summary stats).
Requires data/raw/ckdgen/CKD_overall_EA.dbgap.txt.gz (see README)."""
import gzip
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, "results_tables")
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

RAW = os.path.join(ROOT, "data", "raw", "ckdgen",
                   "CKD_overall_EA.dbgap.txt.gz")
INST = os.path.join(T, "D8_pQTL_instruments.csv")
OUT = os.path.join(T, "D12_NID2_CKDbinary_MR.csv")

print("读取CKD binary summary stats (dbgap格式: Chr Pos RSID A1 A2 Freq1 Effect StdErr P n)...", flush=True)
rows = {}
with gzip.open(RAW, "rt") as f:
    header = f.readline()
    for line in f:
        p = line.split()
        if len(p) < 10:
            continue
        rsid = p[2]
        try:
            rows[rsid] = {
                "a1": p[3], "a2": p[4],
                "freq1": float(p[5]),
                "beta_out": float(p[6]),
                "se_out": float(p[7]),
                "p_out": float(p[8]),
                "n_out": int(float(p[9])),
            }
        except ValueError:
            continue
print(f"  读取 {len(rows)} 个变异", flush=True)

inst = pd.read_csv(INST)
nid2 = inst[inst.gene == "NID2"].copy()
print(f"NID2仪器: {len(nid2)}", flush=True)

recs = []
for _, r in nid2.iterrows():
    snp = r["snp"]
    if snp not in rows:
        continue
    o = rows[snp]
    same = str(r["a1"]).upper() == str(o["a1"]).upper()
    beta_out = o["beta_out"] if same else -o["beta_out"]
    recs.append({
        "snp": snp,
        "beta_pqtl": float(r["beta_pqtl"]),
        "se_pqtl": float(r["se_pqtl"]),
        "a1": r["a1"], "a1_out": o["a1"],
        "beta_out": beta_out, "se_out": o["se_out"],
        "p_out": o["p_out"], "n_out": o["n_out"],
        "freq1": o["freq1"],
        "allele_aligned": same,
    })
df = pd.DataFrame(recs)
print(f"匹配: {len(df)}/{len(nid2)} 仪器", flush=True)

if len(df) == 0:
    print("无匹配仪器")
    sys.exit(1)

bx = df.beta_pqtl.values
bxse = df.se_pqtl.values
by = df.beta_out.values
byse = df.se_out.values

def ivw(by, byse):
    w = 1.0 / byse**2
    beta = np.sum(w * by) / np.sum(w)
    se = np.sqrt(1.0 / np.sum(w))
    return beta, se, 2 * norm.sf(abs(beta / se))

def egger(bx, by, byse):
    X = np.column_stack([np.ones(len(bx)), bx])
    sw = np.sqrt(1.0 / byse**2)
    WX, Wy = X * sw[:, None], by * sw
    coef, *_ = np.linalg.lstsq(WX, Wy, rcond=None)
    resid = Wy - WX @ coef
    dof = len(bx) - 2
    sigma2 = np.sum(resid**2) / dof
    cov = sigma2 * np.linalg.inv(WX.T @ WX)
    se_int = np.sqrt(cov[0, 0])
    return coef[1], coef[0], 2 * norm.sf(abs(coef[0] / se_int))

def wmedian(bx, by, byse):
    ratio = by / bx
    w = 1.0 / (byse / bx)**2
    order = np.argsort(ratio)
    r_s, w_s = ratio[order], w[order]
    cw = np.cumsum(w_s) - 0.5 * w_s
    cw /= np.sum(w_s)
    idx = np.searchsorted(cw, 0.5)
    beta = r_s[min(idx, len(r_s) - 1)]
    se = np.sqrt(1.0 / np.sum(w))
    return beta, 2 * norm.sf(abs(beta / se))

def cochran_q(bx, by, byse):
    ratio = by / bx
    beta_ref = np.sum(ratio / byse**2) / np.sum(1 / byse**2)
    q = np.sum((ratio - beta_ref)**2 / byse**2)
    return q, chi2.sf(q, len(bx) - 1)

beta, se, p_ivw = ivw(by, byse)
slope, intercept, p_int = egger(bx, by, byse)
beta_wm, p_wm = wmedian(bx, by, byse)
q, q_p = cochran_q(bx, by, byse)
log_or = beta / np.log(10)
log_or_se = se / np.log(10)
OR = 10 ** log_or
lo = 10 ** (log_or - 1.96 * log_or_se)
hi = 10 ** (log_or + 1.96 * log_or_se)

out = pd.DataFrame([{
    "gene": "NID2",
    "outcome": "CKDGen_CKD_binary_EA (Wuttke 2019)",
    "n_cases": 41395, "n_controls": 439303,
    "n_instruments": len(df),
    "ivw_beta_logOR": round(log_or, 4),
    "ivw_OR": round(OR, 3),
    "ivw_CI": f"{lo:.3f}-{hi:.3f}",
    "ivw_p": p_ivw,
    "egger_intercept": round(intercept, 4),
    "egger_intercept_p": p_int,
    "weighted_median_beta": round(beta_wm, 4),
    "weighted_median_p": p_wm,
    "cochran_Q": round(q, 2),
    "cochran_Q_p": q_p,
}])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
out.to_csv(OUT, index=False)
print("\n=== NID2 × CKDGen CKD binary（构念匹配测试）===")
for k, v in out.iloc[0].items():
    print(f"  {k}: {v}")
print("\n=== 单仪器明细 ===")
print(df[["snp", "beta_pqtl", "beta_out", "se_out", "p_out",
          "allele_aligned"]].to_string(index=False))
