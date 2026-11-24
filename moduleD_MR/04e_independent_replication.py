# -*- coding: utf-8 -*-
"""Independent-GWAS replication: 453 lead cis-eQTL instruments x CKDGen eGFRcrea
(Wuttke 2019 Nat Genet; GWAS Catalog GCST008059; European ancestry, 42 studies,
n = 567,460). Completely independent of FinnGen (no shared cohorts).

Direction convention: eGFRcrea is continuous, higher = better kidney function.
A gene protective for kidney failure (FinnGen OR < 1) should therefore show a
POSITIVE eGFR association (wald_beta > 0 = higher eGFR per unit expression).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm, binomtest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
TAB = os.path.join(ROOT, "results", "tables")
HITS = os.path.join(ROOT, "data", "external", "ckdgen",
                    "CKDGen_instrument_hits.tsv")

# 1) 解析命中文件（空格分隔；只保留精确rsID匹配，剔除grep子串污染）
d = pd.read_csv(HITS, sep=r"\s+", header=None, on_bad_lines="skip")
d.columns = ["chr", "pos", "rsid", "a1", "a2", "freq1", "effect", "stderr",
             "pval", "n"][:d.shape[1]]
panel = set(open(os.path.join(ROOT, "data", "external", "ckdgen",
                              "instrument_rsids.txt")).read().split())
d = d[d["rsid"].astype(str).isin(panel)].copy()
d = d.drop_duplicates("rsid")
print(f"CKDGen exact matches: {len(d)}/{len(panel)} instruments")

# 2) 合并暴露（D1有453个lead instrument）
d1 = pd.read_csv(os.path.join(TAB, "D1_instruments.csv"))
d2 = pd.read_csv(os.path.join(TAB, "D2_MR_results_kidney.csv"))
m = d1.merge(d, on=None, left_on="snp", right_on="rsid", how="inner")
print(f"merged: {len(m)} instruments")

# 3) 等位基因协调：eQTLGen z对应a1_x（评估等位）；CKDGen beta对应Allele1(a1_y)
same = m["a1_x"].astype(str).str.lower() == m["a1_y"].astype(str).str.lower()
flip = ~same
m["beta_ckd"] = m["effect"].astype(float)
m.loc[flip[m.index], "beta_ckd"] *= -1
m["se_ckd"] = m["stderr"].astype(float)
m["p_ckd"] = m["pval_y"].astype(float)

# 4) 重建暴露beta/se（eQTLGen官方公式；D1无MAF列，用CKDGen欧洲MAF代理）
p = m["freq1"].astype(float).clip(0.01, 0.99)
N = m["n_samples"].astype(float)
denom = np.sqrt(2 * p * (1 - p) * (N + m["z_exp"] ** 2))
m["beta_exp"] = m["z_exp"] / denom
m["se_exp"] = 1.0 / denom

# 5) Wald ratio vs eGFRcrea
m["wald_beta"] = m["beta_ckd"] / m["beta_exp"]
m["se_wald"] = m["se_ckd"] / m["beta_exp"].abs()
z = m["wald_beta"] / m["se_wald"]
m["wald_p"] = 2 * norm.sf(z.abs())

# 6) FinnGen复制状态 + 方向一致性（疾病尺度）
fin = d2.groupby("gene").agg(
    n_sig=("ivw_p", lambda s: int((s < 0.05).sum())),
    med_or=("OR", "median")).reset_index()
fin = fin[fin["n_sig"] >= 2]
m = m.merge(fin, on="gene", how="left")
m["replicated_in_finngen"] = m["n_sig"] >= 2
m["ckd_nominal"] = m["wald_p"] < 0.05
m["direction_concordant"] = np.where(
    m["replicated_in_finngen"],
    np.sign(m["wald_beta"]) == -np.sign(np.log(m["med_or"].clip(lower=1e-9))),
    np.nan)

keep = ["gene", "snp", "z_exp", "beta_exp", "beta_ckd", "se_ckd", "p_ckd",
        "wald_beta", "se_wald", "wald_p", "n_sig", "med_or",
        "replicated_in_finngen", "ckd_nominal", "direction_concordant"]
out = m[keep].sort_values(["replicated_in_finngen", "wald_p"],
                          ascending=[False, True])
out.to_csv(os.path.join(TAB, "D7_independent_CKDGen_replication.csv"),
           index=False)

rep = out[out["replicated_in_finngen"]]
n_rep = len(rep)
n_nom = int(rep["ckd_nominal"].sum())
n_conc = int((rep["direction_concordant"] == True).sum())  # noqa: E712
print(f"\n=== Independent replication (CKDGen eGFRcrea, GCST008059) ===")
print(f"of the 41 FinnGen-replicated genes, {n_rep} have a CKDGen match")
print(f"nominally significant: {n_nom}/{n_rep} (chance ~{0.05 * n_rep:.1f})")
print(f"direction-concordant: {n_conc}/{n_rep}")
if n_rep:
    print(f"binomial p = "
          f"{binomtest(n_nom, n_rep, 0.05, alternative='greater').pvalue:.4f}")
print()
cols = ["gene", "snp", "wald_beta", "wald_p", "n_sig", "med_or",
        "ckd_nominal", "direction_concordant"]
print(rep[cols].to_string(index=False))
print()
print(f"all {len(out)} instruments: CKDGen nominal "
      f"{int(out['ckd_nominal'].sum())}")
