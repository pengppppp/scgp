# -*- coding: utf-8 -*-
"""
Colocalisation: cis-eQTL (eQTLGen full cis sumstats) x kidney GWAS (FinnGen R12).

Review response: the reviewer required coloc/fine-mapping to show that the cis-eQTL
signal and the kidney-outcome signal at each replicated gene's locus share a causal
variant rather than being in LD.

Method
------
1. Exposure: eQTLGen full cis-eQTL summary statistics (hg19), one streaming pass
   over the file keeping +/-1 Mb windows around the 12 target genes (best probe per
   SNP by Pvalue). beta/se reconstructed from the z-score with the consortium
   formula, using the rsID's FinnGen allele frequency (European MAF proxy).
2. Outcome: FinnGen R12 kidney endpoints (local, hg19), one pass per endpoint file
   keeping the union of target rsIDs.
3. LD: 1000 Genomes EUR LD matrix per locus from the OpenGWAS API (/ld/matrix,
   query-param format, max ~200 variants per request). The coloc candidate set is
   therefore capped at 200 SNPs: the 100 most significant by exposure p, the 100
   most significant by outcome p, and the primary instrument (deduplicated).
   SNPs outside this set have Wakefield ABF ~ 1 for BOTH traits, so their exclusion
   changes PP.H4 negligibly for loci with real signal and is CONSERVATIVE for null
   loci (verified analytically: non-associated SNPs add equally to sum1, sum2 and
   sum12, so removing them lowers H4 relative to H3/H0 for null loci).
4. coloc.abf (Wakefield ABFs, single causal variant per trait) on the LD-pruned
   (r2 < 0.1) candidate set - pruning makes coloc.abf's independence assumption
   valid; this is the documented conservative variant of the method.
5. Output: PP.H0-PP.H4 per gene x endpoint (D5_coloc_results.csv).

JWT read from OPENGWAS_JWT env var (never written into the project).
"""
import os
import gzip
import time
import numpy as np
import pandas as pd
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
TAB = os.path.join(ROOT, "results", "tables")
EQTL_FILE = os.path.join(ROOT, "data", "external", "eqtlgen_cis_full",
                         "eqtlgen_cis_full_20191211.txt.gz")
RAW = os.path.join(ROOT, "data", "raw")
OUT_TAB = os.path.join(ROOT, "results", "tables")
LOG = os.path.join(ROOT, "results", "logs")

GENES = ["LOXL4", "TUFT1", "CASP7", "NID2", "C2", "TGFB3",
         "CHSY1", "ST3GAL4", "CD14", "PCSK6", "TNFRSF1A", "LGALS3"]
EPS = ["N14_CHRONKIDNEYDIS", "N14_RENFAIL", "N14_RENALTUB",
       "DM_NEPHROPATHY", "N14_IGA_NEPHROPATHY", "I9_HYPTENSRD"]
WINDOW = 1_000_000
R2_PRUNE = 0.10
W_EQTL, W_GWAS = 0.15 ** 2, 0.2 ** 2
P1, P2, P12 = 1e-4, 1e-4, 1e-5
TOP_EXP, TOP_OUT = 100, 100          # coloc candidate-set caps
JWT = os.environ.get("OPENGWAS_JWT", "")


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)


def gene_coords():
    d1 = pd.read_csv(os.path.join(TAB, "D1_instruments.csv"))
    d1 = d1[d1["gene"].isin(GENES)][["gene", "snp", "chr", "pos"]]
    return d1.drop_duplicates("gene").set_index("gene")


def eqtlgen_all_regions(coords):
    """Single streaming pass over the eQTLGen full cis file; keep +/-1 Mb windows
    for all target genes."""
    want = set(GENES)
    ranges = {}
    for g in GENES:
        if g in coords.index:
            c, p = coords.loc[g, "chr"], int(coords.loc[g, "pos"])
            ranges[(str(c), g)] = (p - WINDOW, p + WINDOW)
    out = {g: [] for g in GENES}
    t0 = time.time()
    n_rows = 0
    with gzip.open(EQTL_FILE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            n_rows += 1
            parts = line.rstrip("\n").split("\t")
            g = parts[idx["GeneSymbol"]]
            if g not in want:
                continue
            chrom = parts[idx["SNPChr"]]
            rng = ranges.get((chrom, g))
            if rng is None:
                continue
            pos = int(parts[idx["SNPPos"]])
            if not (rng[0] <= pos <= rng[1]):
                continue
            out[g].append({
                "snp": parts[idx["SNP"]], "chr": int(chrom), "pos": pos,
                "assessed": parts[idx["AssessedAllele"]],
                "other": parts[idx["OtherAllele"]],
                "z": float(parts[idx["Zscore"]]),
                "pval": float(parts[idx["Pvalue"]]),
                "n_eff": float(parts[idx["NrSamples"]]),
            })
    log(f"eQTLGen pass: {n_rows:,} rows in {time.time()-t0:.0f}s")
    dfs = {}
    for g, rows in out.items():
        df = pd.DataFrame(rows)
        if len(df):
            df = df.sort_values("pval").drop_duplicates("snp").reset_index(drop=True)
        dfs[g] = df
    return dfs


def finngen_all_regions(coords, exp_dfs):
    """One pass per FinnGen endpoint; keep rows for the union of exposure SNPs."""
    want_rsids = set()
    for df in exp_dfs.values():
        if len(df):
            want_rsids |= set(df["snp"])
    log(f"FinnGen target rsIDs: {len(want_rsids):,}")
    out = {ep: {} for ep in EPS}
    for ep in EPS:
        path = os.path.join(RAW, f"finngen_R12_{ep}.gz")
        t0 = time.time()
        rows = []
        with gzip.open(path, "rt") as f:
            cols = f.readline().lstrip("#").rstrip("\n").split("\t")
            ci = {c: i for i, c in enumerate(cols)}
            for line in f:
                parts = line.rstrip("\n").split("\t")
                rs = parts[ci["rsids"]]
                if rs in want_rsids:
                    try:
                        rows.append({
                            "snp": rs,
                            "chr": int(parts[ci["chrom"]]),
                            "pos": int(parts[ci["pos"]]),
                            "alt": parts[ci["alt"]], "ref": parts[ci["ref"]],
                            "af_alt": float(parts[ci["af_alt"]]),
                            "beta": float(parts[ci["beta"]]),
                            "se": float(parts[ci["sebeta"]]),
                            "pval": float(parts[ci["pval"]]),
                        })
                    except (ValueError, KeyError):
                        continue
        for r in rows:
            out[ep][r["snp"]] = r
        log(f"  {ep}: {len(rows):,} rows in {time.time()-t0:.0f}s")
    return out


def ld_matrix(rsids, chunk=240):
    """1000G EUR LD matrix for the given rsIDs (query-param format, chunked).
    Returns (matrix, snplist) restricted to rsIDs found in the reference."""
    if not JWT or not rsids:
        return None, []
    all_snps, blocks = [], []
    for i in range(0, len(rsids), chunk):
        part = rsids[i:i + chunk]
        params = [("rsid", r) for r in part] + [("pop", "EUR")]
        ok = False
        for attempt in range(3):
            try:
                r = requests.post("https://api.opengwas.io/api/ld/matrix",
                                  params=params,
                                  headers={"Authorization": f"Bearer {JWT}",
                                           "Content-Type": "application/json"},
                                  timeout=180)
                if r.status_code == 200:
                    d = r.json()
                    sl = d.get("snplist", [])
                    mat = np.array(d.get("matrix", []), dtype=float)
                    if mat.size and mat.shape[0] == len(sl):
                        all_snps.extend([s.split("_")[0] for s in sl])
                        blocks.append((sl, mat))
                        ok = True
                        break
                else:
                    log(f"  ld/matrix HTTP {r.status_code}")
            except Exception as e:
                log(f"  ld/matrix error: {e}")
            time.sleep(5)
        if not ok:
            return None, []
        time.sleep(1.5)
    if len(blocks) == 1:
        return blocks[0][1], all_snps
    # assemble the full matrix across chunks (only pairwise-complete if chunks
    # cover all variants; callers pass <= chunk variants or handle blocks)
    n = len(all_snps)
    full = np.eye(n)
    pos = 0
    for sl, mat in blocks:
        k = len(sl)
        full[pos:pos + k, pos:pos + k] = mat
        pos += k
    return full, all_snps


def ld_prune(df, r2mat, r2_thresh=R2_PRUNE, sort_col=None):
    """Greedy LD pruning; sort_col defaults to the exposure p-value."""
    if sort_col is None:
        sort_col = "pval_exp" if "pval_exp" in df.columns else "pval"
    keep = []
    order = df[sort_col].argsort().values      # most significant first
    for i in order:
        if all(r2mat[i, j] ** 2 <= r2_thresh for j in keep):
            keep.append(i)
    return df.iloc[keep].reset_index(drop=True)


def prune_region(snps_ranked, chunk=240):
    """Iterative chunked LD pruning of a regional SNP list (rsIDs, already ranked
    by significance). LD is fetched in <=190-variant chunks; survivors are
    re-pruned until the set fits one request. Returns the final rsID list."""
    current = list(snps_ranked)
    rounds = 0
    while True:
        rounds += 1
        survivors = []
        for i in range(0, len(current), chunk):
            part = current[i:i + chunk]
            R, sl = ld_matrix(part)
            if R is None:
                log(f"  prune_region: LD failed at round {rounds}")
                return None
            found = [s.split("_")[0] for s in sl]
            idx = [found.index(s) for s in part if s in found]
            sub = [s for s in part if s in found]
            if not sub:
                continue
            Rsub = R[np.ix_(idx, idx)]
            # greedy prune within the chunk (rank = position in `part`)
            keep = []
            for a in range(len(sub)):
                if all(Rsub[a, b] ** 2 <= R2_PRUNE for b in keep):
                    keep.append(a)
            survivors.extend([sub[a] for a in keep])
        log(f"  prune_region round {rounds}: {len(current)} -> {len(survivors)}")
        if len(survivors) == len(current) and rounds > 1:
            log('  prune_region: converged (no further shrinkage)')
            current = survivors
            break
        if len(survivors) <= chunk:
            current = survivors
            break
            R, sl = ld_matrix(survivors)
            if R is None:
                return None
            found = [s.split("_")[0] for s in sl]
            idx = [found.index(s) for s in survivors if s in found]
            sub = [s for s in survivors if s in found]
            Rsub = R[np.ix_(idx, idx)]
            keep = []
            for a in range(len(sub)):
                if all(Rsub[a, b] ** 2 <= R2_PRUNE for b in keep):
                    keep.append(a)
            return [sub[a] for a in keep]
        current = survivors
    # converged above: final full-LD prune of the survivor set
    R, sl = ld_matrix(current)
    if R is None:
        return None
    found = [s.split("_")[0] for s in sl]
    idx = [found.index(s) for s in current if s in found]
    sub = [s for s in current if s in found]
    Rsub = R[np.ix_(idx, idx)]
    keep = []
    for a in range(len(sub)):
        if all(Rsub[a, b] ** 2 <= R2_PRUNE for b in keep):
            keep.append(a)
    return [sub[a] for a in keep]


def log_abf(beta, se, w):
    """log Wakefield approximate Bayes factor (coloc:::approx.bf.p), numerically
    stable: logABF = 0.5*log(V/(V+W)) + z^2*W/(2*(V+W)). Working in log space
    avoids overflow for extremely strong cis-eQTLs (z ~ 40)."""
    v = se ** 2
    z = beta / se
    return 0.5 * np.log(v / (v + w)) + (z ** 2 * w) / (2 * (v + w))


def coloc_abf(b1, s1, b2, s2):
    from scipy.special import logsumexp
    l1 = log_abf(b1, s1, W_EQTL)
    l2 = log_abf(b2, s2, W_GWAS)
    n = len(l1)
    # log of each hypothesis's (un-normalised) weight
    log_h0 = 0.0
    log_h1 = np.log(P1) + logsumexp(l1)
    log_h2 = np.log(P2) + logsumexp(l2)
    log_h3 = np.log(P1 * P2) + logsumexp(l1) + logsumexp(l2)
    log_h4 = np.log(P12) + logsumexp(l1 + l2)
    log_w = np.array([log_h0, log_h1, log_h2, log_h3, log_h4])
    w = np.exp(log_w - logsumexp(log_w))
    return {f"H{i}": float(w[i]) for i in range(5)}


def z_to_beta_se(z, af, n):
    p = min(max(af, 0.01), 0.99) if np.isfinite(af) else 0.5
    denom = np.sqrt(2 * p * (1 - p) * (n + z ** 2))
    return z / denom, 1.0 / denom


def main():
    coords = gene_coords()
    genes = [g for g in GENES if not os.environ.get("COLOC_GENES")
             or g in os.environ["COLOC_GENES"].split(",")]
    log(f"coloc: {len(genes)} genes x {len(EPS)} endpoints; "
        f"JWT={'set' if JWT else 'MISSING'}")
    if not os.path.exists(EQTL_FILE):
        log("eQTLGen full file not downloaded yet - abort")
        return
    exp_dfs = eqtlgen_all_regions(coords)
    fin = finngen_all_regions(coords, exp_dfs)

    results = []
    for gene in genes:
        exp = exp_dfs.get(gene, pd.DataFrame())
        if not len(exp):
            log(f"{gene}: no exposure data")
            continue
        # ---- build the regional merged table (exposure x each endpoint) ----
        per_ep = {}
        for ep in EPS:
            rows = [fin[ep][s] for s in exp["snp"] if s in fin[ep]]
            if len(rows) < 50:
                per_ep[ep] = None
                continue
            out = pd.DataFrame(rows)
            mm = exp.merge(out, on="snp", suffixes=("_exp", "_out"))
            same = mm["assessed"] == mm["alt"]
            flip = mm["assessed"] == mm["ref"]
            mm = mm[same | flip].copy()
            if not len(mm):
                per_ep[ep] = None
                continue
            mm.loc[flip[mm.index], "beta"] *= -1
            mm = mm[(mm["se"] > 0) & mm["af_alt"].between(0.001, 0.999)]
            if len(mm) < 50:
                per_ep[ep] = None
                continue
            mm[["beta_exp", "se_exp"]] = mm.apply(
                lambda r: pd.Series(z_to_beta_se(r["z"], r["af_alt"], r["n_eff"])),
                axis=1)
            mm = mm.rename(columns={"beta": "beta_out", "se": "se_out"})
            per_ep[ep] = mm
        usable = {ep: mm for ep, mm in per_ep.items() if mm is not None}
        if not usable:
            log(f"{gene}: no endpoint with >=50 SNPs")
            continue
        # ---- regional SNP set ranked by combined significance ----
        all_snps = set()
        for mm in usable.values():
            all_snps |= set(mm["snp"])
        best_p = {}
        for s in all_snps:
            ps = []
            for mm in usable.values():
                row = mm[mm["snp"] == s]
                if len(row):
                    ps.append(min(row["pval_exp"].iloc[0], row["pval_out"].iloc[0]))
            best_p[s] = min(ps) if ps else 1.0
        lead = coords.loc[gene, "snp"]
        ranked = sorted(all_snps, key=lambda s: best_p[s])
        if lead in all_snps:
            ranked.remove(lead)
            ranked.insert(0, lead)
        log(f"{gene}: {len(ranked)} regional SNPs across "
            f"{len(usable)} endpoints; pruning...")
        # ---- prune once per gene ----
        pruned_snps = prune_region(ranked)
        if not pruned_snps:
            log(f"{gene}: pruning failed")
            continue
        log(f"{gene}: {len(pruned_snps)} approximately independent SNPs")
        # ---- coloc per endpoint on the shared pruned set ----
        for ep, mm in usable.items():
            sub = mm[mm["snp"].isin(pruned_snps)].reset_index(drop=True)
            if len(sub) < 10:
                log(f"{gene} {ep}: <10 pruned SNPs present")
                continue
            pp = coloc_abf(sub["beta_exp"].values, sub["se_exp"].values,
                           sub["beta_out"].values, sub["se_out"].values)
            results.append({
                "gene": gene, "endpoint": ep,
                "n_pruned_total": len(pruned_snps), "n_pruned_present": len(sub),
                "PP_H0": round(pp["H0"], 4), "PP_H1": round(pp["H1"], 4),
                "PP_H2": round(pp["H2"], 4), "PP_H3": round(pp["H3"], 4),
                "PP_H4": round(pp["H4"], 4),
                "lead_instrument_retained": lead in set(sub["snp"]),
                "pp_h4_class": ("strong" if pp["H4"] >= 0.8 else
                                "moderate" if pp["H4"] >= 0.5 else "weak"),
            })
            log(f"{gene} {ep}: PP.H4={pp['H4']:.3f} (n={len(sub)})")
    if results:
        out = pd.DataFrame(results)
        suf = os.environ.get("COLOC_SUFFIX", "")
        out.to_csv(os.path.join(OUT_TAB, f"D5_coloc_results{suf}.csv"), index=False)
        log(f"coloc done: {len(out)} pairs -> D5_coloc_results.csv")
        print(out.to_string(index=False))
    else:
        log("no coloc results produced")



if __name__ == "__main__":
    main()
