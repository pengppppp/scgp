# -*- coding: utf-8 -*-
"""
Module D: cis-eQTL Mendelian randomization for the secreted-protein candidate panel.
=========================================================================================
Public data:
  Exposure: eQTLGen phase-1 significant cis-eQTLs (31,684 blood donors)
            https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/
            2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  Outcomes: FinnGen R12 GWAS summary statistics (public, no individual data needed)
            endpoints: N14_CHRONKIDNEYDIS (CKD), N14_RENFAIL (renal failure),
                       N14_RENALTUB (renal tubulo-interstitial diseases),
                       DM_NEPHROPATHY (diabetic nephropathy),
                       N14_IGA_NEPHROPATHY, I9_HYPTENSRD (hypertensive renal disease)

Method notes:
  - For a single-instrument (lead-SNP) MR the Wald-ratio z statistic equals
    sign(beta_exp) * beta_out/se_out, i.e. the p-value is the outcome GWAS p-value
    at that SNP; single-instrument results are therefore treated as prioritisation
    evidence only and are required to replicate across endpoints.
  - Multi-instrument genes additionally get a HEIDI-style heterogeneity test.
  - Replication enrichment is calibrated against an empirical null of random
    MAF-matched SNPs across the same six endpoints.
  - Palindromic SNPs (A/T, C/G) are removed (strand unresolvable without a
    reference panel).
  - MR-Egger regression uses outcome-precision weights (1/se_out^2), per Burgess
    and Thompson (2017).
  - Weighted median uses a weighted bootstrap; weighted mode uses a
    Silverman-bandwidth kernel.
  - MR-PRESSO is an approximation; primary inference relies on IVW plus
    heterogeneity plus replication.

Run:  python 04_MR_causal.py [--multi]
"""
import os
import gzip
import json
import argparse
import urllib.request
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
RAW = os.path.join(ROOT, "data", "raw")
EXT = os.path.join(ROOT, "data", "external")
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
LOGS = os.path.join(ROOT, "results", "logs")
for d in (EXT, PROC, TAB, FIG, LOGS):
    os.makedirs(d, exist_ok=True)

EQTLGEN_URL = ("https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
               "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
EQTLGEN_LOCAL = os.path.join(RAW, "eqtlgen_cis_FDR05.txt.gz")
FINNGEN_BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{}.gz"

KIDNEY_ENDPOINTS = {
    "N14_CHRONKIDNEYDIS": "Chronic kidney disease",
    "N14_RENFAIL": "Renal failure",
    "N14_RENALTUB": "Renal tubulo-interstitial diseases",
    "DM_NEPHROPATHY": "Diabetic nephropathy",
    "N14_IGA_NEPHROPATHY": "IgA nephropathy",
    "I9_HYPTENSRD": "Hypertensive renal disease",
}
P_INSTRUMENT = 5e-8


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def download(url, out):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    log(f"downloading {os.path.basename(out)} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=7200) as r, open(out, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    return out


# ------------------------------------------------------------------ eQTLGen exposure
def load_exposure(genes):
    """Stream eQTLGen file; keep significant cis-eQTLs for candidate genes.

    eQTLGen provides Z-scores (not beta/se). Per the official eQTLGen README
    (SMR_formatted/README_cis_full_SMR), beta and se are reconstructed as:
        beta = z / sqrt(2*p*(1-p)*(N + z^2))
        se   = 1 / sqrt(2*p*(1-p)*(N + z^2))
    with p = allele frequency and N = per-effect sample size (NrSamples column).
    The allele frequency used here is the FinnGen af_alt of the harmonised variant
    (European ancestry; documented approximation in the manuscript).
    """
    download(EQTLGEN_URL, EQTLGEN_LOCAL)
    gene_set = set(genes)
    rows = []
    with gzip.open(EQTLGEN_LOCAL, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {h.strip().lower(): i for i, h in enumerate(header)}

        def c(*names):
            for n in names:
                if n in col:
                    return col[n]
            return None

        i_p = c("pvalue", "p_value")
        i_snp = c("snp")
        i_chr = c("snpchr")
        i_pos = c("snppos")
        i_a1 = c("assessedallele")
        i_a2 = c("otherallele")
        i_z = c("zscore")
        i_gene = c("genesymbol", "gene_symbol")
        i_n = c("nrsamples")
        log(f"eQTLGen header: {header}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                g = parts[i_gene] if i_gene is not None else None
            except Exception:
                continue
            if g not in gene_set:
                continue
            try:
                p = float(parts[i_p]) if i_p is not None else 1.0
            except ValueError:
                continue
            if p > P_INSTRUMENT:
                continue
            try:
                z = float(parts[i_z])
            except (ValueError, IndexError):
                continue
            try:
                n = float(parts[i_n]) if i_n is not None else 31864.0
            except (ValueError, IndexError):
                n = 31864.0
            rows.append({"gene": g, "snp": parts[i_snp], "chr": parts[i_chr],
                         "pos": parts[i_pos], "a1": parts[i_a1], "a2": parts[i_a2],
                         "pval": p, "z_exp": z, "n_samples": n})
    exp = pd.DataFrame(rows)
    log(f"exposure instruments: {len(exp)} eQTLs for "
        f"{exp['gene'].nunique() if len(exp) else 0} genes (p<{P_INSTRUMENT:g})")
    return exp


def select_instruments(exp, multi=False):
    """Lead-SNP per gene (primary) or 1-Mb window pruning (sensitivity)."""
    if exp.empty:
        return exp
    out = []
    for gene, sub in exp.groupby("gene"):
        sub = sub.sort_values("pval")
        if not multi:
            out.append(sub.iloc[[0]])
        else:
            kept, kept_pos = [], []
            for _, r in sub.iterrows():
                pos = int(r["pos"]) if str(r["pos"]).isdigit() else None
                if pos is None or all(abs(pos - p) > 1_000_000 for p in kept_pos):
                    kept.append(r)
                    if pos is not None:
                        kept_pos.append(pos)
            out.append(pd.DataFrame(kept))
    return pd.concat(out, ignore_index=True) if out else exp


# ------------------------------------------------------------------ FinnGen outcome
def load_outcome(endpoint):
    url = FINNGEN_BASE.format(endpoint)
    local = os.path.join(RAW, f"finngen_R12_{endpoint}.gz")
    download(url, local)
    rows = []
    with gzip.open(local, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        col = {h.strip().lower(): i for i, h in enumerate(header)}
        i_rsid = col.get("rsids")
        i_beta = col.get("beta")
        i_se = col.get("sebeta")
        i_p = col.get("pval")
        i_chr = col.get("#chrom", col.get("chrom"))
        i_pos = col.get("pos")
        i_ref = col.get("ref")
        i_alt = col.get("alt")
        i_af = col.get("af_alt")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                rsid = parts[i_rsid] if i_rsid is not None else ";"
                if ";" in rsid:
                    rsid = rsid.split(";")[0]
                af = None
                if i_af is not None:
                    try:
                        af = float(parts[i_af])
                    except (ValueError, IndexError):
                        af = None
                rows.append({"rsid": rsid, "chr": parts[i_chr], "pos": parts[i_pos],
                             "ref": parts[i_ref], "alt": parts[i_alt],
                             "beta_out": float(parts[i_beta]), "se_out": float(parts[i_se]),
                             "pval_out": float(parts[i_p]), "af_alt": af})
            except (ValueError, IndexError):
                continue
    df = pd.DataFrame(rows)
    log(f"  {endpoint}: {len(df)} variants parsed")
    return df


# ------------------------------------------------------------------ MR engine
def harmonize(exp_inst, out_df):
    """Match instruments to outcome by rsID (fallback chr:pos), align alleles to the
    outcome ALT allele, DROP palindromic SNPs, and reconstruct exposure beta/se from
    the eQTLGen z-score using the official eQTLGen formula:
        beta = z / sqrt(2*p*(1-p)*(N + z^2));  se = 1 / sqrt(2*p*(1-p)*(N + z^2))
    p = FinnGen af_alt (European MAF proxy), N = eQTLGen NrSamples.
    """
    if exp_inst.empty or out_df.empty:
        return pd.DataFrame()
    m = exp_inst.merge(out_df, left_on="snp", right_on="rsid", how="inner")
    if m.empty:
        m = exp_inst.merge(out_df, left_on=["chr", "pos"], right_on=["chr", "pos"], how="inner")
    if m.empty:
        return pd.DataFrame()
    m = m.dropna(subset=["z_exp", "beta_out", "se_out"])
    m = m[m["se_out"] > 0]
    # allele alignment: eQTLGen z refers to a1 (assessed allele); FinnGen beta to alt.
    # if the assessed allele equals the outcome REF allele, flip the z sign.
    flip = (m["a1"].astype(str) == m["ref"].astype(str))
    m["z_aligned"] = m["z_exp"].astype(float)
    m.loc[flip, "z_aligned"] = -m.loc[flip, "z_exp"].astype(float)
    # drop palindromic SNPs (strand unresolvable without an LD reference)
    pal = np.array([tuple(sorted([a, b])) in [("A", "T"), ("C", "G")]
                    for a, b in zip(m["a1"].astype(str), m["a2"].astype(str))])
    n_pal = int(pal.sum())
    m = m.loc[~pal].copy()
    # allele frequency: FinnGen af_alt when available, else 0.5 (conservative)
    p = m["af_alt"].astype(float).fillna(0.5)
    p = p.clip(lower=1e-4, upper=0.9999)
    N = m["n_samples"].astype(float).fillna(31864.0)
    z2 = m["z_aligned"] ** 2
    denom = np.sqrt(2 * p * (1 - p) * (N + z2))
    m["beta_exp"] = m["z_aligned"] / denom
    m["se_exp"] = 1.0 / denom
    m = m[(m["se_exp"] > 0) & np.isfinite(m["beta_exp"])]
    m.attrs["n_palindromic_dropped"] = n_pal
    return m


def wald_ratios(m):
    m = m.copy()
    m["beta_xy"] = m["beta_out"] / m["beta_exp"]
    m["se_xy"] = m["se_out"] / m["beta_exp"].abs()
    m["w"] = 1.0 / m["se_xy"] ** 2
    m["F"] = (m["beta_exp"] / m["se_exp"]) ** 2
    return m


def ivw(m):
    w = m["w"].values
    beta = m["beta_xy"].values
    b = np.sum(w * beta) / np.sum(w)
    se_fixed = np.sqrt(1.0 / np.sum(w))
    q = np.sum(w * (beta - b) ** 2)
    df = len(m) - 1
    if df > 0:
        phi = max(q / df, 1e-12)
        w_re = 1.0 / (m["se_xy"].values ** 2 * phi)
        b_re = np.sum(w_re * beta) / np.sum(w_re)
        se_re = np.sqrt(1.0 / np.sum(w_re))
    else:
        b_re, se_re = b, se_fixed
    return {"beta": b, "se": se_fixed, "p": _p_from_beta_se(b, se_fixed),
            "beta_re": b_re, "se_re": se_re, "p_re": _p_from_beta_se(b_re, se_re),
            "Q": q, "Q_df": df, "Q_p": _chi2_sf(q, df) if df > 0 else np.nan}


def mr_egger(m):
    """MR-Egger with OUTCOME-precision weights w = 1/se_out^2 (Burgess & Thompson 2017)."""
    x = m["beta_exp"].values
    y = m["beta_out"].values
    w = 1.0 / m["se_out"].values ** 2
    W = np.diag(w)
    X = np.column_stack([np.ones_like(x), x])
    XtWX = X.T @ W @ X
    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return {"intercept": None, "intercept_se": None, "intercept_p": None,
                "beta": None, "se": None, "p": None}
    coef = XtWX_inv @ X.T @ W @ y
    intercept, slope = coef
    resid = y - X @ coef
    dof = len(m) - 2
    sigma2 = float(resid @ W @ resid / dof) if dof > 0 else np.nan
    cov = sigma2 * XtWX_inv
    se_int = np.sqrt(cov[0, 0])
    se_slope = np.sqrt(cov[1, 1])
    return {"intercept": intercept, "intercept_se": se_int,
            "intercept_p": _p_from_beta_se(intercept, se_int),
            "beta": slope, "se": se_slope, "p": _p_from_beta_se(slope, se_slope)}


def weighted_median(m):
    """Weighted median with a WEIGHTED bootstrap SE (weights = ratio precision)."""
    b = m["beta_xy"].values
    w = m["w"].values
    order = np.argsort(b)
    b, w = b[order], w[order]
    cw = np.cumsum(w) / np.sum(w)
    idx = np.searchsorted(cw, 0.5)
    idx = min(idx, len(b) - 1)
    beta = b[idx]
    rng = np.random.default_rng(42)
    probs = w / w.sum()
    boots = [np.median(rng.choice(b, size=len(b), replace=True, p=probs))
             for _ in range(1000)]
    se = float(np.std(boots))
    return {"beta": float(beta), "se": se, "p": _p_from_beta_se(beta, se)}


def weighted_mode(m):
    """Weighted mode with Silverman-bandwidth kernel density on the Wald-ratio scale."""
    b = m["beta_xy"].values
    w = m["w"].values
    sd = float(np.std(b))
    iqr = float(np.percentile(b, 75) - np.percentile(b, 25))
    bw = 1.06 * min(sd, iqr / 1.34) * len(b) ** (-0.2)
    if not np.isfinite(bw) or bw <= 0:
        bw = max(sd, 1e-6)
    num = np.sum(w * np.exp(-((b[:, None] - b[None, :]) / bw) ** 2), axis=1)
    i = int(np.argmax(num))
    beta = float(b[i])
    se = float(np.sqrt(1.0 / np.sum(w)))
    return {"beta": beta, "se": se, "p": _p_from_beta_se(beta, se)}


def heidi_test(m):
    """HEIDI-style heterogeneity test: Cochran's Q of the Wald ratios (df = n-1).
    Non-significant Q = ratios homogeneous across instruments = consistent with a
    single shared causal variant. Not computable for single-instrument genes."""
    if len(m) < 3:
        return {"Q": None, "Q_p": None, "df": len(m) - 1}
    r = ivw(m)
    return {"Q": r["Q"], "Q_p": r["Q_p"], "df": r["Q_df"]}


def mr_presso(m, n_sim=2000, alpha=0.05):
    """MR-PRESSO outlier detection (approximation; documented in Methods)."""
    r = ivw(m)
    if r["beta"] is None or np.isnan(r["beta"]):
        return {"outliers": [], "global_p": np.nan}
    expected = r["beta"] * m["beta_exp"].values
    obs = m["beta_out"].values
    sse = (obs - expected) ** 2
    rng = np.random.default_rng(42)
    expected_sim = rng.normal(r["beta"], np.abs(r["se"]), size=(n_sim, len(m))) \
        * m["beta_exp"].values[None, :]
    sse_sim = (expected_sim - m["beta_out"].values[None, :]) ** 2
    ps = (sse_sim > sse[None, :]).mean(axis=0)
    outliers = list(m.loc[ps < alpha, "snp"])
    return {"outliers": outliers, "global_p": float(np.min(ps)) if len(ps) else np.nan}


def leave_one_out(m):
    res = []
    for i in range(len(m)):
        sub = m.drop(m.index[i])
        if len(sub) < 1:
            continue
        r = ivw(sub)
        res.append({"omitted_snp": m.iloc[i]["snp"], "beta": r["beta"],
                    "se": r["se"], "p": r["p"]})
    return pd.DataFrame(res)


def _p_from_beta_se(beta, se):
    from scipy.stats import norm
    try:
        if se is None or se == 0 or np.isnan(se):
            return np.nan
    except (TypeError, ValueError):
        return np.nan
    z = beta / se
    return float(2 * norm.sf(abs(z)))


def _chi2_sf(x, df):
    from scipy.stats import chi2
    return float(chi2.sf(x, df)) if df > 0 else np.nan


def run_mr_for_gene(gene, m, endpoint):
    """Full MR battery for one gene-endpoint pair; returns flat result dict."""
    m = wald_ratios(m)
    r_ivw = ivw(m)
    r_egg = mr_egger(m) if len(m) >= 3 else {}
    r_wm = weighted_median(m) if len(m) >= 3 else {}
    r_wmode = weighted_mode(m) if len(m) >= 3 else {}
    heidi = heidi_test(m)
    presso = mr_presso(m) if len(m) >= 4 else {"outliers": [], "global_p": np.nan}
    res = {
        "gene": gene, "endpoint": endpoint, "n_snps": len(m),
        "ivw_beta": r_ivw["beta"], "ivw_se": r_ivw["se"], "ivw_p": r_ivw["p"],
        "ivw_beta_re": r_ivw["beta_re"], "ivw_p_re": r_ivw["p_re"],
        "Q": r_ivw["Q"], "Q_p": r_ivw["Q_p"],
        "heidi_Q": heidi["Q"], "heidi_p": heidi["Q_p"],
        "egger_intercept": r_egg.get("intercept"), "egger_intercept_p": r_egg.get("intercept_p"),
        "egger_beta": r_egg.get("beta"), "egger_p": r_egg.get("p"),
        "wm_beta": r_wm.get("beta"), "wm_p": r_wm.get("p"),
        "wmode_beta": r_wmode.get("beta"), "wmode_p": r_wmode.get("p"),
        "presso_outliers": ",".join(presso["outliers"]),
        "n_presso_outliers": len(presso["outliers"]),
        "snps": ",".join(m["snp"].astype(str)),
        "mean_F": float(np.nanmean(m["F"].values)) if len(m) else np.nan,
    }
    res["OR"] = float(np.exp(res["ivw_beta"])) if res["ivw_beta"] is not None \
        and not np.isnan(res["ivw_beta"]) else np.nan
    lo, hi = res["ivw_beta"] - 1.96 * res["ivw_se"], res["ivw_beta"] + 1.96 * res["ivw_se"]
    res["OR_lci"], res["OR_uci"] = float(np.exp(lo)), float(np.exp(hi))
    return res, m


def plot_four_panel(gene, m, res, endpoint, figpath):
    """Scatter / forest / funnel / leave-one-out for one gene-endpoint pair."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    ax = axes[0, 0]
    ax.errorbar(m["beta_exp"], m["beta_out"], yerr=1.96 * m["se_out"], fmt="o",
                capsize=3, color="#185FA5", ecolor="#85B7EB")
    be = pd.to_numeric(m["beta_exp"], errors="coerce").astype(float)
    xs = np.linspace(float(be.min()), float(be.max()), 50)
    ax.plot(xs, float(res["ivw_beta"]) * xs, color="#A32D2D", lw=1.5, label="IVW")
    egg_b = res.get("egger_beta")
    egg_i = res.get("egger_intercept")
    try:
        egg_b = float(egg_b) if egg_b is not None else np.nan
        egg_i = float(egg_i) if egg_i is not None else 0.0
    except (TypeError, ValueError):
        egg_b, egg_i = np.nan, 0.0
    if np.isfinite(egg_b):
        ax.plot(xs, egg_b * xs + egg_i, color="#854F0B", lw=1, ls="--", label="Egger")
    ax.set_xlabel("eQTL effect on gene expression (beta)")
    ax.set_ylabel(f"GWAS effect on {endpoint} (beta)")
    ax.set_title(f"{gene} - {endpoint}: scatter")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    methods, betas, ses = [], [], []
    for m_name, b_key, s_key in [("IVW", "ivw_beta", "ivw_se"),
                                 ("Egger", "egger_beta", "egger_se"),
                                 ("WM", "wm_beta", "wm_se"),
                                 ("WMode", "wmode_beta", "wmode_se")]:
        b = res.get(b_key, None)
        s = res.get(s_key, None)
        try:
            b = float(b) if b is not None else np.nan
            s = float(s) if s is not None else np.nan
        except (TypeError, ValueError):
            b, s = np.nan, np.nan
        if np.isfinite(b) and np.isfinite(s) and s > 0:
            methods.append(m_name)
            betas.append(b)
            ses.append(s)
    ypos = np.arange(len(methods))
    if methods:
        ax.errorbar(betas, ypos, xerr=1.96 * np.array(ses, dtype=float), fmt="s",
                    capsize=4, color="#0F6E56", ecolor="#5DCAA5", markersize=7)
        ax.set_yticks(ypos)
        ax.set_yticklabels(methods)
        ax.set_ylim(len(methods) - 0.5, -0.5)
    # honest annotation: multi-instrument estimators that are UNDEFINED at this
    # instrument count are named in the empty area instead of leaving it blank
    missing = [nm for nm, bk, sk in [("Egger", "egger_beta", "egger_se"),
                                     ("WM", "wm_beta", "wm_se"),
                                     ("WMode", "wmode_beta", "wmode_se")]
               if nm not in methods]
    if missing:
        ax.text(0.02, 0.05 if len(methods) <= 2 else 0.85,
                "not defined for n=" + str(len(m)) + " instruments://n" +
                ", ".join(missing),
                transform=ax.transAxes, fontsize=7.5, va="bottom", ha="left",
                color="#666666",
                bbox=dict(boxstyle="round,pad=0.3", fc="#F5F5F5", ec="#CCCCCC"))
    ax.axvline(0, color="gray", lw=0.8, ls=":")
    ax.set_title("MR estimates (beta, 95% CI)")
    ax = axes[1, 0]
    ax.scatter(1.0 / m["se_xy"], m["beta_xy"], color="#7F77DD")
    ax.axvline(1.0 / res["ivw_se"], color="#A32D2D", lw=1, ls="--")
    ax.axhline(res["ivw_beta"], color="#A32D2D", lw=1)
    ax.set_xlabel("1/SE (precision)")
    ax.set_ylabel("Wald ratio")
    ax.set_title("funnel plot")
    ax = axes[1, 1]
    loo = leave_one_out(m)
    if len(loo):
        ypos = np.arange(len(loo))
        ax.errorbar(loo["beta"], ypos, xerr=1.96 * loo["se"], fmt="o", capsize=4,
                    color="#993C1D", ecolor="#F0997B", markersize=9,
                    elinewidth=1.6)
        ax.axvline(res["ivw_beta"], color="#A32D2D", lw=1, ls="--")
        ax.set_yticks(ypos)
        ax.set_yticklabels(loo["omitted_snp"], fontsize=6)
        ax.set_title("leave-one-out")
    else:
        ax.text(0.5, 0.5, "leave-one-out not applicable\n(single instrument)",
                ha="center", va="center", fontsize=9)
        ax.set_title("leave-one-out")
    try:
        or_txt = f"OR={float(res['OR']):.2f}"
    except (TypeError, ValueError):
        or_txt = "OR=NA"
    try:
        p_txt = f"p={float(res['ivw_p']):.2g}"
    except (TypeError, ValueError):
        p_txt = "p=NA"
    plt.suptitle(f"{gene} -> {endpoint}  (IVW {or_txt}, {p_txt})", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(figpath, bbox_inches="tight", dpi=150)
    plt.close("all")


# ------------------------------------------------------------------ empirical null
def empirical_null(null_lookup, seed=42):
    """Empirical null for cross-endpoint replication.

    null_lookup: {endpoint: {rsid: pval}} built from ONE candidate rsID set drawn
    once (MAF-stratified to the instrument allele-frequency distribution) and looked
    up in all six FULL FinnGen files. Returns the empirical P(>=k endpoints nominal)
    over the candidate SNPs, replacing the independence assumption (FinnGen endpoints
    are correlated and their p-values are not exactly calibrated to 5%).
    """
    eps = list(null_lookup.keys())
    common = set(null_lookup[eps[0]].keys())
    for ep in eps[1:]:
        common &= set(null_lookup[ep].keys())
    common = sorted(common)
    if not common:
        return None
    counts = np.zeros(len(common), dtype=int)
    for ep in eps:
        pv = null_lookup[ep]
        counts += np.array([1 if pv.get(r, 1.0) < 0.05 else 0 for r in common])
    probs = {k: float((counts >= k).mean()) for k in range(1, len(eps) + 1)}
    log(f"empirical null: {len(common)} candidate SNPs present in all {len(eps)} "
        f"endpoints; P(>=1)={probs[1]:.4f}, P(>=2)={probs[2]:.4f}, P(>=3)={probs[3]:.4f}")
    return {"n_snps": int(len(common)), "probs": probs}



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--multi", action="store_true", help="multi-SNP instruments (window pruning)")
    args = ap.parse_args()

    panel_file = os.path.join(TAB, "B1_secreted_candidate_panel.csv")
    if not os.path.exists(panel_file):
        log("FATAL: run 02_secretome.py first (Module B panel missing)")
        raise SystemExit(1)
    panel = pd.read_csv(panel_file)
    genes = panel["gene"].astype(str).tolist()
    log(f"candidate panel: {len(genes)} genes")

    exp = load_exposure(genes)
    if exp.empty:
        log("FATAL: no eQTL instruments found for panel genes")
        raise SystemExit(1)
    exp_inst = select_instruments(exp, multi=args.multi)
    log(f"instruments selected: {len(exp_inst)} (mode={'multi' if args.multi else 'lead-SNP'})")
    _suf = "_MULTI" if args.multi else ""
    exp_inst.to_csv(os.path.join(TAB, f"D1_instruments{_suf}.csv"), index=False)

    inst_rsids = set(exp_inst["snp"].astype(str))
    all_res, null_lookup = [], {}
    candidate_rsids = None
    for ep, desc in KIDNEY_ENDPOINTS.items():
        try:
            out_df = load_outcome(ep)
        except Exception as e:
            log(f"{ep}: outcome download failed ({e}) -> skipped")
            continue
        # empirical-null lookup: draw ONE candidate rsID set from the first endpoint,
        # MAF-STRATIFIED to match the instrument set's own allele-frequency
        # distribution (nominal-significance probability depends on MAF, so an
        # unmatched null is biased). Then look the candidates up in each full file.
        if candidate_rsids is None:
            # memory-safe: only the two columns needed for sampling are materialised
            # (copying the full 21M-row frame caused MemoryError)
            full = out_df[["rsid", "af_alt"]].copy()
            inst_maf = full.loc[full["rsid"].astype(str).isin(inst_rsids),
                                "af_alt"].dropna()
            inst_maf = inst_maf[(inst_maf > 0.01) & (inst_maf < 0.99)]
            if len(inst_maf) >= 20:
                # instrument MAF distribution -> sampling weights per MAF bin
                bins = np.quantile(inst_maf.values,
                                   np.linspace(0, 1, 11))
                bins = np.unique(bins)
                if len(bins) >= 3:
                    inst_bin = np.digitize(inst_maf.values, bins[1:-1])
                    bin_prop = pd.Series(inst_bin).value_counts(normalize=True)
                    pool = full[full["rsid"].astype(str).str.startswith("rs")]
                    pool = pool[(pool["af_alt"] > 0.01) & (pool["af_alt"] < 0.99)]
                    pool_bin = np.digitize(pool["af_alt"].values, bins[1:-1])
                    n_target = 50000
                    take = []
                    for b, prop in bin_prop.items():
                        n_b = int(round(prop * n_target))
                        sub = pool[pool_bin == b]
                        if len(sub) and n_b > 0:
                            take.append(sub.sample(min(n_b, len(sub)),
                                                   random_state=42)["rsid"])
                    if take:
                        candidate_rsids = set(pd.concat(take).astype(str))
                        log(f"empirical-null candidates: {len(candidate_rsids)} "
                            f"MAF-matched SNPs ({len(bins)-1} bins from instrument "
                            f"MAF distribution)")
            if not candidate_rsids:
                valid = full.loc[full["rsid"].astype(str).str.startswith("rs"),
                                 "rsid"].unique()
                candidate_rsids = set(pd.Series(valid).sample(
                    min(50000, len(valid)), random_state=42))
        sub = out_df[out_df["rsid"].isin(candidate_rsids)]
        null_lookup[ep] = dict(zip(sub["rsid"].astype(str), sub["pval_out"].astype(float)))
        # MR uses only the instrument-matched rows
        out_df = out_df[out_df["rsid"].astype(str).isin(inst_rsids)].copy()
        log(f"{ep}: {len(out_df)}/{len(inst_rsids)} instrument variants present in outcome")
        for gene in genes:
            gi = exp_inst[exp_inst["gene"] == gene]
            if gi.empty:
                continue
            m = harmonize(gi, out_df)
            if m.empty:
                continue
            res, m_full = run_mr_for_gene(gene, m, ep)
            res["endpoint_desc"] = desc
            res["n_palindromic_dropped"] = m_full.attrs.get("n_palindromic_dropped", 0)
            all_res.append(res)
            if res["ivw_p"] < 0.05:
                figp = os.path.join(FIG, f"D1_MR_{gene}_{ep}.png")
                try:
                    plot_four_panel(gene, m_full, res, ep, figp)
                except Exception as e:
                    log(f"  figure failed for {gene}-{ep}: {e}")
        log(f"{ep}: {sum(1 for r in all_res if r['endpoint'] == ep)} gene tests done")
    res_df = pd.DataFrame(all_res)
    if res_df.empty:
        log("FATAL: no MR tests executed")
        raise SystemExit(1)
    from scipy.stats import false_discovery_control, binomtest
    res_df["bonferroni_threshold"] = 0.05 / len(res_df)
    res_df = res_df.sort_values("ivw_p").reset_index(drop=True)
    res_df["bh_fdr"] = false_discovery_control(res_df["ivw_p"].values, method="bh")
    res_df["significant_bonferroni"] = res_df["ivw_p"] < res_df["bonferroni_threshold"]
    res_df.to_csv(os.path.join(TAB, f"D2_MR_results_kidney{_suf}.csv"), index=False)

    # ---- cross-endpoint replication with EMPIRICAL null ----
    n_genes_tested = res_df["gene"].nunique()
    nom = res_df[res_df["ivw_p"] < 0.05]
    rep = nom.groupby("gene")["endpoint"].nunique()
    obs_ge2 = int((rep >= 2).sum())
    null = empirical_null(null_lookup)
    if null:
        p_ge2 = null["probs"].get(2, 0.05)
        exp_ge2 = n_genes_tested * p_ge2
        bt = binomtest(obs_ge2, n_genes_tested, p_ge2, alternative="greater")
        log(f"REPLICATION (empirical null): {obs_ge2} genes nominal in >=2 endpoints; "
            f"expected {exp_ge2:.1f} (P(>=2)={p_ge2:.4f}); binomial p={bt.pvalue:.3g}")
        pd.DataFrame([
            {"metric": "genes_tested", "value": n_genes_tested},
            {"metric": "genes_nominal_ge2_endpoints", "value": obs_ge2},
            {"metric": "empirical_P_ge2", "value": p_ge2},
            {"metric": "expected_ge2_by_chance", "value": exp_ge2},
            {"metric": "binom_p_greater", "value": bt.pvalue},
        ]).to_csv(os.path.join(TAB, f"D3_replication_enrichment{_suf}.csv"), index=False)
        with open(os.path.join(PROC, "moduleD_params.json"), "w", encoding="utf-8") as f:
            json.dump({"p_instrument": P_INSTRUMENT,
                       "mode": "multi" if args.multi else "lead-SNP",
                       "endpoints": KIDNEY_ENDPOINTS, "n_genes_tested": n_genes_tested,
                       "empirical_null": null, "palindromic_removed": True,
                       "empirical_null_maf_matched": True,
                       "egger_weights": "1/se_out^2 (Burgess & Thompson 2017)",
                       "eQTLGen_url": EQTLGEN_URL,
                       "run": datetime.now().isoformat(timespec="seconds")}, f, indent=1)
    log(f"Module D done: {len(res_df)} gene-endpoint tests; "
        f"{(res_df['ivw_p'] < 0.05).sum()} nominal, "
        f"{res_df['significant_bonferroni'].sum()} Bonferroni")


if __name__ == "__main__":
    main()
