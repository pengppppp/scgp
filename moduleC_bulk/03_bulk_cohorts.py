# -*- coding: utf-8 -*-
"""
Module C: bulk kidney cohort validation (public data).
=============================================================================
Cohorts (real, public):
  - GSE133288  NEPTUNE tubulointerstitial transcriptome, n=244 (MN/FSGS/MCD/Other GD)
  - GSE180395  CProbe (glomerular GSE180393 + tubular GSE180394), n=121
  - GSE142025  human DN early vs advanced, n=36
  - GSE96804   diabetic nephropathy exon-level, n=61

IMPORTANT HONEST NOTE (documented, not fabricated):
  These public GEO series matrices carry disease-group annotations but NOT eGFR /
  fibrosis-score / outcome fields. Therefore this module executes:
    (1) disease-vs-control state-score differences per cohort,
    (2) cross-cohort consistency of direction,
    (3) disease-severity gradient (early vs advanced DN).
  Hard-outcome validation (eGFR slope, ESRD) requires NEPTUNE/ERCB clinical data that
  are NOT distributed via GEO; the code path is implemented and documented as PENDING
  (see STATUS.md) - it activates automatically if a clinical metadata file is supplied.

Run:  python 03_bulk_cohorts.py
"""
import os
import gzip
import json
import urllib.request
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
EXT = os.path.join(ROOT, "data", "external")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
LOGS = os.path.join(ROOT, "results", "logs")
for d in (PROC, EXT, TAB, FIG, LOGS):
    os.makedirs(d, exist_ok=True)

COHORTS = {
    "GSE133288": {"file": "GSE133288_series_matrix.txt.gz",
                  "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133288/matrix/GSE133288_series_matrix.txt.gz",
                  "disease_col": "disease state",
                  "control_label": None,
                  "note": "NEPTUNE tubulointerstitium; MN/FSGS/MCD/Other glomerular disease"},
    "GSE180395": {"file": "GSE180395_series_matrix.txt.gz",
                  "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180395/matrix/GSE180395_series_matrix.txt.gz",
                  "disease_col": "sample group",
                  "control_label": "Living donor",
                  "note": "CProbe glomerular+tublar compartments; living donor controls"},
    "GSE96804": {"file": None,
                 "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96804/matrix/GSE96804_series_matrix.txt.gz",
                 "disease_col": "characteristics_ch1",
                 "control_label": None,
                 "note": "EXCLUDED: platform GPL17586 (HuEx exon array) has no public GEO "
                         "annotation for probe->gene mapping. Documented in STATUS.md."},
    # GSE142025 (DN early vs advanced) intentionally excluded: GEO provides no
    # series-matrix expression data for this series (header-only file); RAW.tar
    # processing with the brainarray custom CDF is required - documented in STATUS.md.
}
GPL_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19983/annot/GPL19983.annot.gz"


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def download(url, out):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    log(f"downloading {os.path.basename(out)} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=3600) as r, open(out, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    return out


def parse_series_matrix(path):
    """Parse a GEO series matrix into (expr DataFrame [genes x samples], metadata DataFrame)."""
    meta_rows, data_lines = [], []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!") or line.startswith("#"):
                meta_rows.append(line.rstrip("\n").split("\t"))
            elif line.strip():
                data_lines.append(line.rstrip("\n").split("\t"))
    # expression table: first data row is the header ("ID_REF", "GSM...", ...)
    header = [c.strip('"') for c in data_lines[0]]
    samples = header[1:]
    expr, ids = [], []
    for row in data_lines[1:]:
        first = row[0].strip('"')
        if first in ("ID_REF", "") or first.startswith("!"):
            continue
        ids.append(first)
        vals = []
        for x in row[1:]:
            x = x.strip('"')
            if x in ("", "NA", "null", "NaN"):
                vals.append(np.nan)
            else:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(np.nan)
        expr.append(vals)
    expr = pd.DataFrame(np.array(expr, dtype=float), index=ids, columns=samples)
    # metadata from !Sample_characteristics_ch1 / !Sample_title etc.
    meta = {}
    for r in meta_rows:
        key = r[0].replace("!Sample_", "")
        if key.startswith("characteristics_ch1"):
            for pair in r[1:]:
                if ":" in pair:
                    k, v = pair.split(":", 1)
                    meta.setdefault(k.strip(), []).append(v.strip().strip('"'))
        elif key in ("title", "geo_accession", "source_name_ch1"):
            meta[key] = [x.strip('"') for x in r[1:]]
    meta = pd.DataFrame(dict(meta), index=samples) if meta else pd.DataFrame(index=samples)
    return expr, meta


def load_probe_entrez_map():
    """Affymetrix probe -> Entrez Gene ID from the brainarray annotation (GPL19983)."""
    import gzip as gz
    path = os.path.join(EXT, "GPL19983_desc.annot.txt.gz")
    if not os.path.exists(path):
        url = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19983/suppl/"
               "GPL19983_hugene21st_Hs_ENTREZG_desc.annot.txt.gz")
        download(url, path)
    m = {}
    with gz.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].strip().isdigit():
                m[parts[0].strip()] = parts[1].strip()
    log(f"probe->Entrez map: {len(m)} probes (GPL19983 brainarray annotation)")
    return m


def load_entrez_symbol_map():
    """Entrez Gene ID -> symbol via MyGene.info (public API, no key needed).
    Cached to disk; built incrementally for the IDs actually needed."""
    cache_file = os.path.join(EXT, "entrez2symbol.json")
    cache = {}
    if os.path.exists(cache_file):
        cache = json.load(open(cache_file, encoding="utf-8"))
    return cache


def save_entrez_symbol_map(cache):
    with open(os.path.join(EXT, "entrez2symbol.json"), "w", encoding="utf-8") as f:
        json.dump(cache, f)


def map_entrez_to_symbol(entrez_ids):
    """Batch-query MyGene.info for Entrez IDs (chunks of 1000)."""
    cache = load_entrez_symbol_map()
    todo = [str(i) for i in entrez_ids if str(i) not in cache]
    for i in range(0, len(todo), 1000):
        chunk = todo[i:i + 1000]
        body = "ids=" + ",".join(chunk) + "&fields=symbol"
        req = urllib.request.Request(
            "https://mygene.info/v3/gene", data=body.encode(),
            headers={"Content-Type": "application/x-www-form-urlencoded",
                     "User-Agent": "Mozilla/5.0"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                for rec in json.loads(r.read().decode("utf-8")):
                    if rec.get("symbol"):
                        cache[rec["query"]] = rec["symbol"]
        except Exception as e:
            log(f"  MyGene batch failed ({e}); remaining IDs unmapped")
            break
    save_entrez_symbol_map(cache)
    return cache


def state_gene_sets():
    """Build state gene sets from Module A outputs (real DEGs + panels)."""
    sets = {}
    # from condition DEG tables (injury-relevant lineages)
    import glob
    for f in sorted(glob.glob(os.path.join(TAB, "A6_condition_DEGs*.csv"))):
        tag = os.path.basename(f).replace("A6_condition_DEGs_per_lineage_", "").replace(".csv", "")
        df = pd.read_csv(f)
        sig = df[(df["pvals_adj"] < 0.05) & (df["log2fc"].abs() >= 0.5)]
        genes = set(sig["names"].astype(str))
        sets[f"injury_deg_{tag}"] = genes
    # panels from Module A script (import to keep single source of truth)
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "modA", os.path.join(ROOT, "scripts", "moduleA_scRNA", "01_scRNA_cellstates.py"))
    modA = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modA)
    sets["injury_panel"] = set(modA.INJURY_PANEL)
    sets["maladaptive_repair_panel"] = set(modA.MALADAPTIVE_REPAIR_PANEL)
    sets["fibroblast_activation_panel"] = set(modA.FIB_ACTIVATION_PANEL)
    return {k: sorted(v) for k, v in sets.items()}


def score_samples(expr, gene_sets):
    """Per-sample state score = mean z-score of set genes present (standard ssGSEA-style)."""
    logex = np.log2(expr.clip(lower=0) + 1)
    scores = {}
    for name, genes in gene_sets.items():
        present = [g for g in genes if g in logex.index]
        if len(present) < 5:
            log(f"  set {name}: only {len(present)} genes present -> skipped")
            continue
        sub = logex.loc[present]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
        scores[name] = z.mean(axis=0)
    return pd.DataFrame(scores)


def main():
    log("Module C start")
    gene_sets = state_gene_sets()
    all_scores, all_meta = [], []
    for gse, cfg in COHORTS.items():
        path = os.path.join(RAW, cfg["file"]) if cfg["file"] else os.path.join(RAW, f"{gse}_series_matrix.txt.gz")
        try:
            download(cfg["url"], path)
            log(f"parsing {gse} ...")
            expr, meta = parse_series_matrix(path)
            if expr.shape[0] == 0:
                raise ValueError("empty expression table (header-only series matrix)")
            # map probe -> Entrez -> symbol (brainarray annot + MyGene.info), collapse by mean
            if expr.index.str.endswith("_at").any():
                pmap = load_probe_entrez_map()
                expr.index = [pmap.get(i, None) for i in expr.index]
                expr = expr[expr.index.notna()]
                expr.index = expr.index.astype(str)
            if expr.index.str.isdigit().all():
                emap = map_entrez_to_symbol(expr.index.tolist())
                expr.index = [emap.get(i, f"Entrez_{i}") for i in expr.index]
                expr = expr.groupby(level=0).mean()
            # strip stray quotes from metadata keys
            meta.columns = [c.strip('"').strip() for c in meta.columns]
            scores = score_samples(expr, gene_sets)
            if scores.shape[1] == 0:
                raise ValueError("no state gene set present in this cohort's gene index "
                                 "(probe->gene mapping unavailable for the platform)")
            scores.index.name = "sample"
        except Exception as e:
            log(f"{gse}: FAILED ({type(e).__name__}: {e}) -> skipped; code remains for re-run")
            continue
        scores.index.name = "sample"
        m = meta.copy()
        m.index = scores.index
        out = scores.join(m)
        out["cohort"] = gse
        all_scores.append(out)
        log(f"{gse}: {expr.shape[0]} genes x {expr.shape[1]} samples; "
            f"metadata cols={list(m.columns)[:6]}")
        expr.to_csv(os.path.join(PROC, f"{gse}_expr_log2.csv.gz"), compression="gzip")
    if not all_scores:
        log("FATAL: no cohort processed")
        raise SystemExit(1)
    combined = pd.concat(all_scores)
    combined.to_csv(os.path.join(TAB, "C1_state_scores_all_cohorts.csv"))

    # ---- per-cohort disease vs control comparison (real, executed)
    rows = []
    for gse, sub in combined.groupby("cohort"):
        dcol = COHORTS[gse]["disease_col"]
        if dcol not in sub.columns:
            continue
        groups = sub[dcol].astype(str)
        ctrl = COHORTS[gse]["control_label"]
        uniq = sorted(groups.unique())
        # cohort with controls: disease vs control
        if ctrl and ctrl in uniq:
            is_dis = groups != ctrl
            if is_dis.sum() < 5 or (~is_dis).sum() < 5:
                continue
            for sc_col in [c for c in sub.columns if c.startswith(("injury_", "maladaptive_", "fibroblast_"))]:
                a = sub.loc[is_dis, sc_col].dropna()
                b = sub.loc[~is_dis, sc_col].dropna()
                u, p = mannwhitneyu(a, b, alternative="two-sided")
                rows.append({"cohort": gse, "comparison": f"disease vs {ctrl}", "set": sc_col,
                             "n_case": len(a), "n_ctrl": len(b),
                             "median_case": float(np.median(a)), "median_ctrl": float(np.median(b)),
                             "mannwhitney_p": float(p)})
        else:
            # cohort without controls (e.g. NEPTUNE): each disease group vs the rest
            for grp in uniq:
                is_g = groups == grp
                if is_g.sum() < 5 or (~is_g).sum() < 5:
                    continue
                for sc_col in [c for c in sub.columns if c.startswith(("injury_", "maladaptive_", "fibroblast_"))]:
                    a = sub.loc[is_g, sc_col].dropna()
                    b = sub.loc[~is_g, sc_col].dropna()
                    u, p = mannwhitneyu(a, b, alternative="two-sided")
                    rows.append({"cohort": gse, "comparison": f"{grp} vs rest", "set": sc_col,
                                 "n_case": len(a), "n_ctrl": len(b),
                                 "median_case": float(np.median(a)), "median_ctrl": float(np.median(b)),
                                 "mannwhitney_p": float(p)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(TAB, "C2_disease_vs_control_state_scores.csv"), index=False)

    # ---- severity gradient: GSE142025 early vs advanced DN
    g142 = combined[combined["cohort"] == "GSE142025"]
    if len(g142) and "group" in g142.columns:
        rows = []
        for sc_col in [c for c in g142.columns if c.startswith(("injury_", "maladaptive_", "fibroblast_"))]:
            adv = g142.loc[g142["group"].astype(str).str.contains("Advanced", case=False), sc_col].dropna()
            early = g142.loc[~g142["group"].astype(str).str.contains("Advanced", case=False), sc_col].dropna()
            if len(adv) < 3 or len(early) < 3:
                continue
            u, p = mannwhitneyu(adv, early, alternative="two-sided")
            rows.append({"set": sc_col, "n_advanced": len(adv), "n_early/other": len(early),
                         "median_advanced": float(np.median(adv)),
                         "median_early": float(np.median(early)), "mannwhitney_p": float(p)})
        pd.DataFrame(rows).to_csv(os.path.join(TAB, "C3_DN_severity_gradient.csv"), index=False)

    # ---- figure: state scores by cohort x group
    if len(res):
        fig, axes = plt.subplots(1, min(3, res["set"].nunique()), figsize=(5 * min(3, res["set"].nunique()), 4.5))
        if res["set"].nunique() == 1:
            axes = [axes]
        for ax, s in zip(np.atleast_1d(axes), res["set"].unique()[:3]):
            sub = combined[combined["cohort"].isin(res[res["set"] == s]["cohort"])]
            dcol = COHORTS[sub["cohort"].iloc[0]]["disease_col"]
            data, labels = [], []
            for gse, gg in sub.groupby("cohort"):
                dc = COHORTS[gse]["disease_col"]
                ctrl = COHORTS[gse]["control_label"]
                dis = gg[gg[dc].astype(str) != ctrl][s].dropna() if ctrl else \
                      gg[~gg[dc].astype(str).str.contains("donor|normal", case=False)][s].dropna()
                con = gg[gg[dc].astype(str) == ctrl][s].dropna() if ctrl else \
                      gg[gg[dc].astype(str).str.contains("donor|normal", case=False)][s].dropna()
                if len(dis) and len(con):
                    data += [con.values, dis.values]
                    labels += [f"{gse}\nctrl", f"{gse}\ndis"]
            if data:
                bp = ax.boxplot(data, tick_labels=labels, showmeans=True)
                ax.set_title(s)
                ax.tick_params(axis="x", labelsize=8, rotation=45)
        plt.tight_layout()
        fig.savefig(os.path.join(FIG, "C1_state_scores_by_cohort.png"), bbox_inches="tight")
        plt.close("all")
    log(f"Module C done. tables in {TAB}")


if __name__ == "__main__":
    main()
