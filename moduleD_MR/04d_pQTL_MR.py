# -*- coding: utf-8 -*-
r"""
pQTL-MR replication: plasma-protein instruments (deCODE SomaScan cis-pQTL) x
FinnGen kidney endpoints.

Review response: the reviewer asked for protein-level (pQTL) replication of the
eQTL-MR findings, because blood cis-eQTLs may not track plasma protein abundance
for secreted proteins.

Inputs
------
- pQTL extract files produced by _work/convert_decode_pqtl.py from the deCODE
  SomaScan release (bucket largescaleplasma-2023, prefix final_somascan_smp/,
  SMP-normalised, hg38), one CSV per gene in D:\133\pQTL_extract\.
  Columns: rsid, pos, a1, a2, beta_pqtl, se_pqtl, p_pqtl, mlog10p_pqtl,
  imp_maf, source_file.  Cis window = gene +/- 1 Mb (deCODE definition) on
  GRCh38 coordinates.  No UK Biobank data are used (the Olink-UKB prefixes in
  the bucket are excluded by project policy; FinnGen outcomes only).
- FinnGen R12 summary statistics (read-only source in the JTM data folder).

Analysis
--------
1. Instrument pool per gene: biallelic SNPs (single-nucleotide alleles, a
   single rsID) inside the cis window, ranked by deCODE -log10(p).  The Pval
   column underflows to 0 for the strongest signals, so -log10(p) is the
   ranking statistic.  Indels and multi-rsID rows stay in the extract for
   audit but are excluded from the instrument pool.
2. Lead instrument: the most significant pool variant that is present and
   allele-consonant in ALL six FinnGen endpoints.  Variant matching is by
   rsID; when an rsID is absent from FinnGen (dbSNP-version drift between the
   deCODE 2023 and FinnGen R12 releases), a chromosome+position match is
   used, with the complement-allele convention as a strand fallback.  Match
   route, the lead's rank in the pool, and any per-endpoint fallback are
   logged and carried into the output table.
3. Allele harmonisation vs the biallelic FinnGen GWAS (flip exposure beta
   when the pQTL effect allele is the FinnGen ref allele).
4. Wald ratio per gene x endpoint: beta = b_out / b_pqtl, se = se_out /
   |b_pqtl|, p from z = beta/se; OR = exp(beta) (direction/rank only).
5. Concordance with the eQTL-MR primary analysis (D2 table): same-sign
   effects and significance pattern.
6. Output: D6_pQTL_MR_replication.csv + a console summary.

Caveats logged with the results: single-SNP Wald p == outcome-association p
(mathematical identity); palindromic lead SNPs are flagged; for proteins
whose strongest cis variants are tied at the deCODE p-value resolution limit
(Pval underflowed to 0 - NID2, C2, LGALS3) the instrument is the most
significant tied variant available in FinnGen and the tie-group size is
reported per row; genes absent from the deCODE SomaScan panel (LOXL4, CASP7,
ST3GAL4, PCSK6, CHSY1) are reported as missing, not imputed.
"""
import os
import gzip
import json
import re
import time

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)                                  # this working copy (writable)
TAB = os.path.join(ROOT, "results", "tables")
RAW = r"D:\Journal of Translational Medicine\data\raw"  # FinnGen: read-only source
PQTL_DIR = os.path.join(ROOT, "pQTL_extract")      # deCODE extracts (D:\133)
COORDS = r"D:\133\_work\gene_coords_hg38.json"     # gene -> hg38 chromosome
OUT_TAB = os.path.join(ROOT, "results", "tables")
LOG = os.path.join(ROOT, "logs")

GENES = ["LOXL4", "TUFT1", "CASP7", "NID2", "C2", "TGFB3",
         "CHSY1", "ST3GAL4", "CD14", "PCSK6", "TNFRSF1A", "LGALS3"]
EPS = ["N14_CHRONKIDNEYDIS", "N14_RENFAIL", "N14_RENALTUB",
       "DM_NEPHROPATHY", "N14_IGA_NEPHROPATHY", "I9_HYPTENSRD"]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}    # complement (strand fallback)
RS_RE = re.compile(r"^rs\d+$")


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(LOG, "moduleD_v10_pQTL_MR.log"), "a",
              encoding="utf-8") as f:
        f.write(line + "\n")


# Columns that must NEVER be picked by substring matching (identifier /
# positional columns).  Without this guard the candidate "p" substring-matches
# "pos" and the lead SNP silently becomes the lowest-coordinate SNP.
NON_VALUE_COLS = {"pos", "position", "chrom", "chr", "chromosome",
                  "rsid", "rsids", "snp", "variant", "a1", "a2",
                  "effectallele", "otherallele", "alt", "ref",
                  "source_file", "gene", "name"}


def find_col(df, candidates):
    # 1) exact match, in candidate priority order
    for c in candidates:
        for col in df.columns:
            if col.lower() == c.lower():
                return col
    # 2) substring match, skipping identifier / positional columns
    for c in candidates:
        for col in df.columns:
            if col.lower() in NON_VALUE_COLS:
                continue
            if c.lower() in col.lower():
                return col
    return None


def load_pqtl(gene):
    """Load the deCODE extract for one gene and build the instrument pool.

    Returns (pool DataFrame ranked by -log10(p), n_all_cis_rows).  The pool is
    restricted to biallelic SNPs carrying a single rsID.
    """
    p = os.path.join(PQTL_DIR, f"{gene}_decode.csv")
    if not os.path.exists(p):
        log(f"{gene}: extract file not found: {p}")
        return None, 0
    df = pd.read_csv(p)
    snp_c = find_col(df, ["rsid", "SNP", "rs", "variant"])
    b_c = find_col(df, ["beta_pqtl", "beta", "effect", "b", "beta_out"])
    se_c = find_col(df, ["se_pqtl", "se", "std_err", "stderr"])
    p_c = find_col(df, ["p_pqtl", "pval", "p_value", "pvalue", "p"])
    m_c = find_col(df, ["mlog10p_pqtl", "mlog10p", "minus_log10_pval"])
    maf_c = find_col(df, ["imp_maf", "impmaf", "maf"])
    a1_c = find_col(df, ["a1", "effect_allele", "effectallele"])
    a2_c = find_col(df, ["a2", "other_allele", "otherallele"])
    pos_c = find_col(df, ["pos", "position"])
    log(f"{gene}: column map snp={snp_c} beta={b_c} se={se_c} p={p_c} "
        f"mlog10p={m_c} maf={maf_c} a1={a1_c} a2={a2_c} pos={pos_c}")
    if snp_c is None or b_c is None or p_c is None or pos_c is None \
            or a1_c is None or a2_c is None:
        log(f"{gene}: columns not recognised {list(df.columns)}")
        return None, 0
    if m_c is None:
        # derive -log10(p) from p (guard against p == 0 underflow)
        df["_mlog10p"] = -np.log10(df[p_c].clip(lower=1e-300))
        m_c = "_mlog10p"
    d = pd.DataFrame({
        "snp": df[snp_c].astype(str).str.strip(),
        "pos": pd.to_numeric(df[pos_c], errors="coerce"),
        "beta_pqtl": pd.to_numeric(df[b_c], errors="coerce"),
        "se_pqtl": pd.to_numeric(df[se_c], errors="coerce") if se_c else np.nan,
        "p_pqtl": pd.to_numeric(df[p_c], errors="coerce"),
        "mlog10p_pqtl": pd.to_numeric(df[m_c], errors="coerce"),
        "imp_maf": pd.to_numeric(df[maf_c], errors="coerce") if maf_c else np.nan,
        "a1": df[a1_c].astype(str).str.strip().str.upper(),
        "a2": df[a2_c].astype(str).str.strip().str.upper(),
    }).dropna(subset=["pos", "beta_pqtl", "p_pqtl", "mlog10p_pqtl"])
    n_all = len(d)
    pool = d[d["snp"].str.match(RS_RE) &
             d["a1"].isin(list("ACGT")) & d["a2"].isin(list("ACGT"))]
    pool = (pool.sort_values("mlog10p_pqtl", ascending=False)
                 .reset_index(drop=True))
    log(f"{gene}: {n_all} cis rows, {len(pool)} biallelic single-rsid SNPs "
        f"in the instrument pool")
    if pool.empty:
        return None, n_all
    lead = pool.iloc[0]
    log(f"{gene}: pool lead (pre-FinnGen) {lead['snp']} "
        f"(-log10p={lead['mlog10p_pqtl']:.4g})")
    return pool, n_all


def finngen_scan(pools):
    """One pass per FinnGen endpoint; collect every pool variant present.

    Matching routes: rsID (exact, incl. comma-split multi-rsID fields), then
    chromosome+position (dbSNP-version drift fallback).  Returns
    {ep: {(gene, key): (alt, ref, af_alt, beta, se, pval)}} with key either
    ("rsid", rsid) or ("pos", pos_int).
    """
    want_rs, want_pos = {}, {}
    for gene, (pool, _) in pools.items():
        chrom = CHROM[gene]
        for pos in pool["pos"].astype(int):
            want_pos[(chrom, str(pos))] = gene
        for snp in pool["snp"]:
            want_rs[snp] = gene
    out = {ep: {} for ep in EPS}
    for ep in EPS:
        path = os.path.join(RAW, f"finngen_R12_{ep}.gz")
        if not os.path.exists(path):
            log(f"FinnGen file missing: {path}")
            continue
        t0, n = time.time(), 0
        with gzip.open(path, "rt") as f:
            cols = f.readline().lstrip("#").rstrip("\n").split("\t")
            ci = {c: i for i, c in enumerate(cols)}
            ichrom, ipos = ci["chrom"], ci["pos"]
            irs, ialt, iref = ci["rsids"], ci["alt"], ci["ref"]
            ib, ise, ipv, iaf = ci["beta"], ci["sebeta"], ci["pval"], ci["af_alt"]
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= iaf:
                    continue
                rs = parts[irs]
                gene, key = want_rs.get(rs), None
                if gene is not None:
                    key = ("rsid", rs)
                elif "," in rs:
                    for tok in rs.split(","):
                        gene = want_rs.get(tok)
                        if gene is not None:
                            key = ("rsid", tok)
                            break
                if gene is None:
                    g = want_pos.get((parts[ichrom], parts[ipos]))
                    if g is not None:
                        gene, key = g, ("pos", int(parts[ipos]))
                if gene is None:
                    continue
                try:
                    info = (parts[ialt], parts[iref], float(parts[iaf]),
                            float(parts[ib]), float(parts[ise]),
                            float(parts[ipv]))
                except ValueError:
                    continue
                out[ep].setdefault((gene, key), info)
                n += 1
        log(f"FinnGen {ep}: {n} pool variants matched "
            f"({time.time() - t0:.0f}s scan)")
    return out


def allele_action(cand, info):
    """How the pQTL exposure allele maps onto the FinnGen alt/ref coding.

    Returns "same" (no flip), "flip" (negate exposure beta) or None (alleles
    not consonant - direct or complement)."""
    a1, a2 = cand["a1"], cand["a2"]
    alt, ref = info[0], info[1]
    if {a1, a2} == {alt, ref}:
        return "flip" if a1 == ref else "same"
    if {COMP[a1], COMP[a2]} == {alt, ref}:
        return "flip" if COMP[a1] == ref else "same"
    return None


def select_lead(gene, pool, fin):
    """Most significant pool variant usable in all six endpoints.

    Returns (lead_idx, {ep: (cand, info, route, action)}) or (None, {}).
    """
    cands = pool.to_dict("records")
    per_ep = {}
    for ep in EPS:
        m = fin.get(ep, {})
        for i, c in enumerate(cands):
            info = m.get((gene, ("rsid", c["snp"])))
            route = "rsid"
            if info is None:
                info = m.get((gene, ("pos", int(c["pos"]))))
                route = "pos"
            if info is None:
                continue
            act = allele_action(c, info)
            if act is None:
                continue
            per_ep[ep] = (i, c, info, route, act)
            break
    if not per_ep:
        return None, {}
    idxs = {v[0] for v in per_ep.values()}
    if len(per_ep) == len(EPS) and len(idxs) == 1:
        i = idxs.pop()
        c = per_ep[EPS[0]][1]
        routes = {v[3] for v in per_ep.values()}
        acts = {v[4] for v in per_ep.values()}
        tie = int((pool["mlog10p_pqtl"] == c["mlog10p_pqtl"]).sum())
        tie_note = (f" (tied with {tie - 1} other pool variant(s) at the same "
                    f"-log10(p); deCODE Pval underflowed to 0)"
                    if tie > 1 else "")
        log(f"{gene}: lead {c['snp']} (pool rank {i + 1}/{len(cands)}, "
            f"-log10p={c['mlog10p_pqtl']:.4g}, MAF={c['imp_maf']:.4g})"
            f"{tie_note} usable in all {len(EPS)} endpoints "
            f"(match route(s): {'/'.join(sorted(routes))}, "
            f"harmonisation: {'/'.join(sorted(acts))})")
        sel = {ep: (per_ep[ep][1], per_ep[ep][2], per_ep[ep][3], per_ep[ep][4],
                    i + 1, tie)
               for ep in EPS}
        return i, sel
    # per-endpoint fallback: endpoints disagree on the first usable variant
    log(f"{gene}: no single variant usable in all endpoints "
        f"(first-usable ranks: "
        f"{ {ep: v[0] + 1 for ep, v in per_ep.items()} }) - "
        f"using the most significant usable variant per endpoint")
    sel = {}
    for ep, v in per_ep.items():
        tie = int((pool["mlog10p_pqtl"] == v[1]["mlog10p_pqtl"]).sum())
        sel[ep] = (v[1], v[2], v[3], v[4], v[0] + 1, tie)
    return None, sel


def main():
    global CHROM
    if not os.path.isdir(PQTL_DIR):
        log(f"pQTL extract directory not found: {PQTL_DIR}")
        return
    CHROM = {g: str(c["chrom"]).strip().lower().replace("chr", "")
             for g, c in json.load(open(COORDS, encoding="utf-8")).items()}

    try:
        d2 = pd.read_csv(os.path.join(TAB, "D2_MR_results_kidney.csv"))
        eqtl_sig = {}
        for _, r in d2[d2["ivw_p"] < 0.05].iterrows():
            eqtl_sig[(r["gene"], r["endpoint"])] = {
                "ivw_p": r["ivw_p"], "OR": r["OR"]}
    except Exception as e:                                   # noqa: BLE001
        log(f"WARNING: could not load D2 eQTL-MR table ({e}); "
            f"concordance columns will be empty")
        eqtl_sig = {}

    pools, n_all = {}, {}
    missing = [g for g in GENES
               if not os.path.exists(os.path.join(PQTL_DIR, f"{g}_decode.csv"))]
    if missing:
        log(f"genes without deCODE pQTL data (not on the SomaScan panel): "
            f"{missing}")
    for gene in GENES:
        pool, n = load_pqtl(gene)
        if pool is not None:
            pools[gene] = (pool, n)
            n_all[gene] = n
    if not pools:
        log("no pQTL instrument pools - aborting")
        return

    log("NOTE: with a single-SNP instrument the Wald-ratio p-value is "
        "mathematically identical to the SNP-outcome association p-value; "
        "treat pQTL-MR p as the outcome-association p at that locus")

    fin = finngen_scan(pools)

    rows = []
    for gene in GENES:
        if gene not in pools:
            continue
        pool = pools[gene][0]
        lead_idx, sel = select_lead(gene, pool, fin)
        if not sel:
            log(f"{gene}: no usable lead variant in any FinnGen endpoint "
                f"- gene skipped")
            continue
        for ep in EPS:
            if ep not in sel:
                log(f"{gene}/{ep}: no usable variant - pair skipped")
                continue
            cand, info, route, act, rank, tie = sel[ep]
            a1, a2 = cand["a1"], cand["a2"]
            alt, ref, af_alt, b_out, se_out, p_out = info
            if act == "flip":
                b_exp = -cand["beta_pqtl"]
            else:
                b_exp = cand["beta_pqtl"]
            if cand["se_pqtl"] <= 0 or se_out <= 0:
                log(f"{gene}/{ep}: non-positive SE - pair skipped")
                continue
            # Wald ratio (exposure: pQTL; outcome: FinnGen)
            b = b_out / b_exp
            se = se_out / abs(cand["beta_pqtl"])
            z = b / se
            p = 2 * norm.sf(abs(z))
            key = (gene, ep)
            e = eqtl_sig.get(key)
            rows.append({
                "gene": gene, "resource": "decode", "endpoint": ep,
                "pqtl_lead_snp": cand["snp"], "pqtl_pos": int(cand["pos"]),
                "pqtl_a1": a1, "pqtl_a2": a2,
                "pqtl_beta": cand["beta_pqtl"], "pqtl_se": cand["se_pqtl"],
                "pqtl_p": cand["p_pqtl"],
                "pqtl_mlog10p": cand["mlog10p_pqtl"],
                "pqtl_imp_maf": cand["imp_maf"],
                "lead_rank_in_pool": rank,
                "lead_tie_group_size": tie,
                "lead_in_all_endpoints": lead_idx is not None,
                "match_route": route,
                "pqtl_lead_palindromic": {a1, a2} in ({"A", "T"}, {"C", "G"}),
                "outcome_alt": alt, "outcome_ref": ref,
                "allele_harmonised": "flipped" if act == "flip" else "same",
                "outcome_beta": b_out, "outcome_p": p_out,
                "wald_ratio": b, "OR": float(np.exp(b)), "wald_p": p,
                "eqtl_mr_nominal_sig": e is not None,
                "eqtl_mr_p": e["ivw_p"] if e else None,
                "direction_concordant": (e is not None and
                                         np.sign(b) == np.sign(np.log(e["OR"]))),
            })
    if not rows:
        log("no pQTL-MR results - check the extract files")
        return
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT_TAB, "D6_pQTL_MR_replication.csv"), index=False)
    log(f"pQTL-MR done: {len(out)} gene-endpoint pairs -> "
        f"D6_pQTL_MR_replication.csv")
    sig = out[out["wald_p"] < 0.05]
    log(f"pQTL-MR nominal significant: {len(sig)} pairs, "
        f"{sig['gene'].nunique()} genes")
    conc = out[out["eqtl_mr_nominal_sig"]]
    if len(conc):
        log(f"among pairs nominal in eQTL-MR ({len(conc)}): "
            f"{int(conc['direction_concordant'].sum())} direction-concordant, "
            f"{int((conc['wald_p'] < 0.05).sum())} also nominal in pQTL-MR")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
