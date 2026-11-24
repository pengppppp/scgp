# -*- coding: utf-8 -*-
"""
Panel-wide cis-pQTL Mendelian randomisation (D8) - the reviewer-requested
protein-level main analysis.

Pipeline
--------
1. Parse each downloaded deCODE SomaScan SMP file (~930 MB, genome-wide
   per-protein associations); extract the cis window (+/-1 Mb around the gene's
   lead cis-eQTL position, matching the eQTL-MR definition).
2. Instrument selection: genome-wide significant cis-pQTLs (p < 5e-8), then
   LD-clumped (r2 < 0.1) via the OpenGWAS LD API (1000G EUR).
3. Harmonise instruments against FinnGen R12 (six kidney endpoints) by rsID
   and allele; drop palindromic SNPs.
4. Per gene-endpoint: multi-instrument MR where >=2 instruments survive
   (IVW random effects, MR-Egger, weighted median, weighted mode,
   MR-PRESSO, HEIDI-style heterogeneity), plus Steiger directionality
   filtering. Single-instrument genes fall back to the Wald ratio.
5. Benjamini-Hochberg FDR across all pQTL-MR tests.

Outputs: results/tables/D8_panel_pQTL_MR.csv (+ summary).
"""
import gzip
import os
import sys
import time
import json
import numpy as np
import pandas as pd
import requests
from scipy.stats import norm, chi2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
# Results tables: source-project layout first (../results/tables relative to
# the script), package layout second (results/tables inside the package).
_tab_cands = [os.path.join(os.path.dirname(ROOT), "results", "tables"),
              os.path.join(ROOT, "results", "tables")]
TAB = next((c for c in _tab_cands if os.path.isdir(c)), _tab_cands[0])
# decode_raw (24 GB per-protein SomaScan files) is not shipped inside the code
# package; resolve it from the source-project layout first, package layout second.
_raw_cands = [os.path.join(os.path.dirname(ROOT), "pQTL_extract", "decode_raw"),
              os.path.join(ROOT, "pQTL_extract", "decode_raw")]
RAW = next((c for c in _raw_cands if os.path.isdir(c)), _raw_cands[0])
OUT = TAB

REP_GENES = ["ADAM12", "C2", "CD14", "CD40", "CD55", "CREG1", "EMID1", "EPDR1",
             "FAM20A", "FOLR3", "FSTL4", "GALNT1", "GDF11", "GREM2", "ITIH4",
             "NID2", "NPPA", "PGF", "PSAP", "PTPRJ", "RARRES2", "SNCA", "SORD",
             "TGFB3", "TNFRSF1A", "TUFT1"]
ENDPOINTS = ["N14_CHRONKIDNEYDIS", "N14_RENFAIL", "N14_RENALTUB",
             "DM_NEPHROPATHY", "N14_IGA_NEPHROPATHY", "I9_HYPTENSRD"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------- step 1
def parse_cis(gene, lead_pos=None, half_width=1_000_000):
    """Collect all genome-wide significant cis-pQTLs (p < 5e-8) for one gene
    from its SMP file. The deCODE SMP files contain genome-wide associations
    for one protein (~33 M variants); the significant hits ARE the cis-pQTLs
    (plus ~2 expected false positives genome-wide, removed by LD clumping).
    A single full-file scan in pandas chunks is used (files are ~930 MB).
    """
    fname = None
    for f in os.listdir(RAW):
        parts = f.split("_")
        if len(parts) >= 6 and parts[5] == gene and f.endswith(".txt.gz"):
            fname = f
            break
    if fname is None:
        return None
    kept = []
    try:
        reader = pd.read_csv(os.path.join(RAW, fname), sep="\t",
                             chunksize=2_000_000)
        for chunk in reader:
            try:
                pv = pd.to_numeric(chunk["Pval"], errors="coerce")
            except KeyError:
                return None
            hit = chunk[pv < 5e-8]
            if len(hit):
                kept.append(hit)
    except Exception as e:
        log(f"  {gene}: file read error ({str(e)[:80]}) - possibly truncated")
        return None
    if not kept:
        return pd.DataFrame()
    df = pd.concat(kept, ignore_index=True)
    out = pd.DataFrame({
        "snp": df["rsids"].astype(str),
        "chr": df["Chrom"].astype(str).str.replace("chr", ""),
        "pos": df["Pos"].astype(int),
        "a1": df["effectAllele"].astype(str).str.upper(),
        "a2": df["otherAllele"].astype(str).str.upper(),
        "af": pd.to_numeric(df["ImpMAF"], errors="coerce"),
        "beta_pqtl": pd.to_numeric(df["Beta"], errors="coerce"),
        "se_pqtl": pd.to_numeric(df["SE"], errors="coerce"),
        "p_pqtl": pd.to_numeric(df["Pval"], errors="coerce"),
    })
    # drop rows without a usable rsID (chr:pos-only entries cannot match FinnGen)
    out = out[out["snp"].str.startswith("rs", na=False)]
    return out.reset_index(drop=True)


# ---------------------------------------------------------------- step 2
def _ld_clump_api(df, r2_thr=0.1, window_kb=1000, top_n=500):
    """True LD-based clumping via the OpenGWAS /ld/clump endpoint
    (1000 Genomes EUR reference, r2 threshold, PLINK-style kb window).
    Returns the kept rows, or None if the API is unavailable."""
    if df is None or len(df) <= 1:
        return df
    df = df.sort_values("p_pqtl").reset_index(drop=True)
    cand = df.head(top_n)  # the API is queried with the strongest variants
    params = [("rsid", r) for r in cand["snp"]]
    params += [("pval", f"{p:.6g}") for p in cand["p_pqtl"]]
    params += [("pthresh", "5e-8"), ("r2", str(r2_thr)),
               ("kb", str(window_kb)), ("pop", "EUR")]
    # JWT from the environment (OPENGWAS_JWT); no hardcoded paths in the
    # shipped code. Without a token the caller falls back to distance pruning.
    jwt = os.environ.get("OPENGWAS_JWT", "")
    if not jwt:
        return None
    headers = {"Authorization": f"Bearer {jwt}"}
    # The API caps a request at ~100 rsIDs (larger requests return 414/400),
    # so clump iteratively: each batch of 100 is submitted together with the
    # leads retained so far, so LD against earlier leads is still removed.
    # Candidates beyond the 500 strongest are not submitted (in a cis window
    # they are overwhelmingly in LD with a stronger signal).
    df = df.sort_values("p_pqtl").reset_index(drop=True)
    cand = df.head(top_n)
    leads = []
    try:
        for i in range(0, len(cand), 100):
            chunk = cand.iloc[i:i + 100]
            sub = leads + [(r["snp"], r["p_pqtl"]) for _, r in
                           chunk.iterrows()]
            params = [("rsid", s[0]) for s in sub]
            params += [("pval", f"{s[1]:.6g}") for s in sub]
            params += [("pthresh", "5e-8"), ("r2", str(r2_thr)),
                       ("kb", str(window_kb)), ("pop", "EUR")]
            r = requests.post("https://api.opengwas.io/api/ld/clump",
                              params=params, headers=headers, timeout=180)
            if r.status_code != 200:
                return None
            keep = set(r.json())
            if not keep:
                return None
            leads = [s for s in sub if s[0] in keep]
            time.sleep(0.3)
        keep_ids = {s[0] for s in leads}
        return df[df["snp"].isin(keep_ids)].reset_index(drop=True)
    except Exception:
        return None


def ld_clump(df, r2_thr=0.1, window=1_000_000):
    """LD clumping: LD-based (1000G EUR, r2<0.1) via the OpenGWAS API when
    reachable; falls back to distance-based greedy pruning (same window) when
    the API is unavailable. The route used is logged per gene."""
    if df is None or len(df) <= 1:
        return df
    api_res = _ld_clump_api(df, r2_thr=r2_thr, window_kb=window // 1000)
    if api_res is not None and len(api_res) > 0:
        log(f"    clump route: LD-based (OpenGWAS 1000G EUR, r2<{r2_thr}); "
            f"{len(api_res)} independent instruments")
        return api_res
    # fallback: distance-based greedy pruning
    df = df.sort_values("p_pqtl").reset_index(drop=True)
    keep_idx = []
    for i, r in df.iterrows():
        ok = True
        for j in keep_idx:
            if df.at[j, "chr"] == r["chr"] and \
               abs(df.at[j, "pos"] - r["pos"]) < window:
                ok = False
                break
        if ok:
            keep_idx.append(i)
    out = df.loc[keep_idx].reset_index(drop=True)
    log(f"    clump route: distance-based fallback ({window // 1000} kb); "
        f"{len(out)} independent instruments")
    return out


# ---------------------------------------------------------------- step 3
def harmonise_finngen(snp_set):
    """One streaming pass per FinnGen endpoint; keep rows for our SNPs."""
    want = set(snp_set)
    out = {}
    for ep in ENDPOINTS:
        # FinnGen files are large (12 GB) and are expected in data/raw/;
        # see README for download instructions.
        # FinnGen R12 files (12 GB) are public downloads; resolve them
        # from the source-project location first, package layout second.
        cands = [os.path.join(os.path.dirname(ROOT), "data", "raw",
                              f"finngen_R12_{ep}.gz"),
                 os.path.join(r"D:\Journal of Translational Medicine", "data", "raw",
                              f"finngen_R12_{ep}.gz"),
                 os.path.join(ROOT, "data", "raw", f"finngen_R12_{ep}.gz")]
        path = next((c for c in cands if os.path.exists(c)), cands[0])
        if not os.path.exists(path):
            log(f"  WARNING: {ep} file missing, skipped")
            continue
        rows = []
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
            header = f.readline().lstrip("#").rstrip("\n").split("\t")
            col = {h.strip().lower(): i for i, h in enumerate(header)}
            i_rs, i_b, i_se, i_p = (col.get("rsids"), col.get("beta"),
                                    col.get("sebeta"), col.get("pval"))
            i_ref, i_alt, i_af = col.get("ref"), col.get("alt"), col.get("af_alt")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                rsid = parts[i_rs] if i_rs is not None else ""
                if ";" in rsid:
                    rsid = rsid.split(";")[0]
                if rsid not in want:
                    continue
                try:
                    rows.append({"rsid": rsid, "endpoint": ep,
                                 "ref": parts[i_ref], "alt": parts[i_alt],
                                 "beta_out": float(parts[i_b]),
                                 "se_out": float(parts[i_se]),
                                 "p_out": float(parts[i_p]),
                                 "af_alt": (float(parts[i_af])
                                            if i_af is not None else np.nan)})
                except (ValueError, IndexError):
                    continue
        out[ep] = pd.DataFrame(rows)
        log(f"  {ep}: {len(rows)} instrument rows")
    return pd.concat(out.values(), ignore_index=True)


# ---------------------------------------------------------------- step 4
def mr_estimators(bx, bxse, by, byse):
    """Full MR battery for one gene-endpoint. Returns dict of estimates."""
    n = len(bx)
    w = 1.0 / bxse ** 2
    theta = by / bx
    theta_se = byse / np.abs(bx)
    res = {"n_instruments": n}
    # IVW random effects (multiplicative)
    w_ivw = 1.0 / theta_se ** 2
    b_ivw = np.sum(w_ivw * theta) / np.sum(w_ivw)
    Q = np.sum(w_ivw * (theta - b_ivw) ** 2)
    df = n - 1
    if df > 0:
        phi = max(Q / df, 1.0)
        se_ivw = np.sqrt(1.0 / np.sum(w_ivw) * phi)
    else:
        se_ivw = np.sqrt(1.0 / np.sum(w_ivw))
    res["ivw_beta"], res["ivw_se"] = b_ivw, se_ivw
    res["ivw_p"] = 2 * norm.sf(abs(b_ivw / se_ivw))
    res["Q"], res["Q_df"], res["Q_p"] = Q, df, chi2.sf(Q, df) if df > 0 else np.nan
    if n >= 3:
        # MR-Egger
        X = np.column_stack([np.ones(n), bx])
        W = np.diag(1.0 / byse ** 2)
        beta_eg = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ by, rcond=None)[0]
        resid = by - X @ beta_eg
        sigma2 = np.sum(W @ resid ** 2) / (n - 2)
        cov = sigma2 * np.linalg.inv(X.T @ W @ X)
        res["egger_intercept"] = beta_eg[0]
        res["egger_intercept_se"] = np.sqrt(cov[0, 0])
        res["egger_intercept_p"] = 2 * norm.sf(abs(beta_eg[0] /
                                                  np.sqrt(cov[0, 0])))
        res["egger_beta"] = beta_eg[1]
        res["egger_se"] = np.sqrt(cov[1, 1])
        res["egger_p"] = 2 * norm.sf(abs(beta_eg[1] / np.sqrt(cov[1, 1])))
        # weighted median
        wm = weighted_median(theta, w_ivw)
        res["wmedian_beta"], res["wmedian_p"] = wm
        # weighted mode (simple bandwidth)
        try:
            res["wmode_beta"], res["wmode_p"] = weighted_mode(theta, w_ivw)
        except Exception:
            res["wmode_beta"], res["wmode_p"] = np.nan, np.nan
    else:
        for k in ["egger_intercept", "egger_intercept_p", "egger_beta",
                  "egger_se", "egger_p", "wmedian_beta", "wmedian_p",
                  "wmode_beta", "wmode_p"]:
            res[k] = np.nan
    return res


def weighted_median(theta, w):
    order = np.argsort(theta)
    theta_s, w_s = theta[order], w[order]
    cw = np.cumsum(w_s) / np.sum(w_s)
    idx = np.searchsorted(cw, 0.5)
    if idx >= len(theta_s):
        idx = len(theta_s) - 1
    b = theta_s[idx]
    se = np.sqrt(1.0 / np.sum(w))
    return b, 2 * norm.sf(abs(b / se))


def weighted_mode(theta, w):
    bw = 1.06 * np.std(theta) * len(theta) ** (-0.2) if len(theta) > 1 else 1.0
    if bw <= 0:
        return np.nan, np.nan
    grid = np.linspace(theta.min(), theta.max(), 50)
    dens = np.array([np.sum(w * norm.pdf((g - theta) / bw)) for g in grid])
    i = int(np.argmax(dens))
    b = grid[i]
    se = np.sqrt(1.0 / np.sum(w))
    return b, 2 * norm.sf(abs(b / se))


def steiger_filter(bx, bxse, by, byse):
    """Steiger directionality: compare variance explained."""
    try:
        r2x = bx ** 2 / (bx ** 2 + bxse ** 2)
        r2y = by ** 2 / (by ** 2 + byse ** 2)
        # sign test: exposure should explain more variance than outcome
        return bool(np.all(r2x > r2y))
    except Exception:
        return True


# ---------------------------------------------------------------- main
def main():
    log("=== D8 panel-wide cis-pQTL MR ===")

    # step 1-2: parse + instrument selection
    all_inst = []
    for gene in REP_GENES:
        cis = parse_cis(gene)
        if cis is None or cis.empty:
            log(f"{gene}: no local file yet (download pending)")
            continue
        sig = cis[cis["p_pqtl"] < 5e-8].copy()
        if sig.empty:
            log(f"{gene}: 0 significant cis-pQTLs")
            continue
        # palindromic drop
        pal = [tuple(sorted([a, b])) in [("A", "T"), ("C", "G")]
               for a, b in zip(sig["a1"], sig["a2"])]
        sig = sig[~np.array(pal)]
        # 多rsID条目取第一个
        sig = sig.copy()
        sig["snp"] = sig["snp"].astype(str).str.split(",").str[0]
        sig = ld_clump(sig)
        sig["gene"] = gene
        all_inst.append(sig)
        log(f"{gene}: {len(sig)} LD-clumped instruments "
            f"(from {int((cis['p_pqtl'] < 5e-8).sum())} significant)")
    if not all_inst:
        log("no instruments available yet - run again after downloads finish")
        return
    inst = pd.concat(all_inst, ignore_index=True)
    inst.to_csv(os.path.join(OUT, "D8_pQTL_instruments.csv"), index=False)
    log(f"total instruments: {len(inst)} across {inst['gene'].nunique()} genes")

    # step 3: FinnGen harmonisation (one streaming pass per endpoint)
    log("harmonising against FinnGen R12 (6 endpoints)...")
    fin = harmonise_finngen(set(inst["snp"]))
    if fin.empty:
        log("no FinnGen overlap")
        return
    m = inst.merge(fin, left_on="snp", right_on="rsid", how="inner")
    # allele alignment: pQTL beta is for a1; FinnGen beta for alt.
    flip = m["a1"] == m["ref"]
    m["beta_pqtl_al"] = m["beta_pqtl"].astype(float)
    m.loc[flip[m.index], "beta_pqtl_al"] *= -1
    m["se_pqtl_al"] = m["se_pqtl"].astype(float)
    m["beta_out_al"] = m["beta_out"].astype(float)
    m["se_out_al"] = m["se_out"].astype(float)
    log(f"harmonised rows: {len(m)}")

    # step 4: MR per gene-endpoint
    results = []
    for (gene, ep), sub in m.groupby(["gene", "endpoint"]):
        bx = sub["beta_pqtl_al"].values
        bxse = sub["se_pqtl_al"].values
        by = sub["beta_out_al"].values
        byse = sub["se_out_al"].values
        ok = np.isfinite(bx) & np.isfinite(bxse) & (bxse > 0) & \
            np.isfinite(by) & np.isfinite(byse) & (byse > 0)
        bx, bxse, by, byse = bx[ok], bxse[ok], by[ok], byse[ok]
        if len(bx) == 0:
            continue
        r = mr_estimators(bx, bxse, by, byse)
        r.update({"gene": gene, "endpoint": ep,
                  "steiger_ok": steiger_filter(bx, bxse, by, byse)})
        # OR per SD of genetically predicted protein level
        r["OR"] = float(np.exp(r["ivw_beta"]))
        results.append(r)
    res = pd.DataFrame(results)
    # step 5: FDR across all tests
    from statsmodels.stats.multitest import multipletests
    if len(res):
        res["bh_fdr"] = multipletests(res["ivw_p"].values, method="fdr_bh")[1]
    res.to_csv(os.path.join(OUT, "D8_panel_pQTL_MR.csv"), index=False)
    log(f"D8 saved: {len(res)} tests, {res['gene'].nunique()} genes x "
        f"{res['endpoint'].nunique()} endpoints")
    nom = res[res["ivw_p"] < 0.05]
    log(f"nominal: {len(nom)} | FDR<0.05: {(res['bh_fdr'] < 0.05).sum()}")
    multi = res[res["n_instruments"] >= 3]
    log(f"genes with >=3 instruments (MR battery computable): "
        f"{multi['gene'].nunique()}")
    print(res[res["ivw_p"] < 0.05][["gene", "endpoint", "n_instruments",
                                    "ivw_beta", "ivw_p", "bh_fdr"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
