# -*- coding: utf-8 -*-
"""
Module A: Disease/injury cell-state definition from REAL human kidney single-cell RNA-seq.
==========================================================================================
Datasets (real, public):
  A) GSE242425 - human kidney organoids, hypoxia-induced AKI-like injury (d20) and
     maladaptive repair (d25); 3 hypoxia lines (h1-h3) + 3 normoxia lines (n1-n3).
     Bio: "Hypoxic injury triggers maladaptive repair in human kidney organoids" (2025).
  B) GSE209781 - human kidney cortex snRNA-seq, 3 diabetic kidney disease (DKD) vs 3
     normal (NM); 51,315 cells.
     Bio: "A single-cell transcriptomic atlas characterizes the kidney in DKD" (2025).

Usage:
  python 01_scRNA_cellstates.py --dataset GSE242425
  python 01_scRNA_cellstates.py --dataset GSE209781
Outputs: QC report, UMAP, lineage annotation, cluster markers,
condition DEGs per lineage, injury/maladaptive-repair state scores and gene sets
(feed Module B). All parameters logged to params_moduleA.json.
"""
import os
import sys
import json
import glob
import tarfile
import argparse
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# ------------------------------------------------------------------ paths
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
C_CTRL, C_CASE = "#8C8C8C", "#D65F5F"
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "results", "figures")
TAB = os.path.join(ROOT, "results", "tables")
for d in (PROC, FIG, TAB):
    os.makedirs(d, exist_ok=True)

# ------------------------------------------------------------------ parameters
PARAMS = {
    "GSE242425": {
        "replicate_key": "line",   # biological replicate = iPSC line (3 per condition)
        "download_base": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242425/suppl/",
        "min_genes_per_cell": 200,
        "min_cells_per_gene": 3,
        "max_mt_pct": 25,
        "max_counts_per_cell": 40000,
        "n_top_hvg": 3000,
        "n_pcs": 30,
        "neighbors_k": 15,
        "leiden_resolution": 1.0,
        "lineage_margin": 0.02,
        "contrast": ("hypoxia", "normoxia"),
        "contrast_label": "hypoxia-injury/repair vs normoxia",
    },
    "GSE209781": {
        "replicate_key": "sample", # biological replicate = donor (3 per condition)
        "download_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209781/suppl/GSE209781_RAW.tar",
        "min_genes_per_cell": 200,
        "min_cells_per_gene": 3,
        "max_mt_pct": 20,          # snRNA: cytoplasmic mRNA depleted -> relaxed threshold
        "max_counts_per_cell": 20000,
        "n_top_hvg": 3000,
        "n_pcs": 30,
        "neighbors_k": 15,
        "leiden_resolution": 1.0,
        "lineage_margin": 0.02,
        "contrast": ("DKD", "normal"),
        "contrast_label": "DKD vs normal",
    },
}

LINEAGE_MARKERS = {
    "Podocyte":            ["NPHS1", "NPHS2", "PODXL"],
    "Proximal tubule":     ["SLC34A1", "CUBN", "SLC5A2", "SLC13A3"],
    "Loop of Henle":       ["SLC12A1", "UMOD", "CLDN16"],
    "Distal tubule":       ["SLC12A3", "CALD1"],
    "Collecting duct":     ["AQP2", "SCNN1G", "ATP6V0D2"],
    "Endothelial":         ["PECAM1", "VWF", "CLDN5", "EMCN"],
    "Fibroblast":          ["COL1A1", "COL3A1", "PDGFRB", "DCN"],
    "Myofibroblast":       ["ACTA2", "TAGLN", "CCN2"],
    "Macrophage":          ["CD68", "CSF1R", "MRC1", "CD14"],
    "Monocyte":            ["FCN1", "S100A8", "S100A9", "VCAN"],
    "T cell":              ["CD3D", "CD3E", "TRAC"],
    "NK cell":             ["NKG7", "GNLY", "KLRD1"],
    "B cell":              ["MS4A1", "CD79A", "BANK1"],
    "Proliferating":       ["MKI67", "TOP2A", "CENPF"],
    "Nephron progenitor":  ["SIX2", "CITED1", "PAX2"],
}

# injury / maladaptive-repair marker panels (literature-consensus;
# see STATUS report for citations - kept as scoring panels, NOT as outcome claims)
INJURY_PANEL = ["HAVCR1", "VCAM1", "ICAM1", "SOX9", "GDF15", "MMP7", "TNFRSF12A", "NFKBIZ"]
MALADAPTIVE_REPAIR_PANEL = ["CCN2", "TGFB1", "COL1A1", "COL3A1", "S100A8", "S100A9", "C3", "CCL2"]
FIB_ACTIVATION_PANEL = ["CCN2", "COL1A1", "COL3A1", "ACTA2", "TAGLN", "POSTN", "FNDC1", "SPARC"]


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def load_GSE242425():
    """Load the 12 h5 matrices (h1-3/n1-3 x d20/d25) of GSE242425."""
    files = sorted(glob.glob(os.path.join(RAW, "GSE242425", "*.h5"))) or \
            sorted(glob.glob(os.path.join(RAW, "GSE242425_*_sample_filtered_feature_bc_matrix.h5")))
    if not files:
        log("FATAL: no GSE242425 h5 files found. Run scripts/utils/download_parallel.sh first.")
        sys.exit(1)
    adatas = []
    for fp in files:
        base = os.path.basename(fp)
        # e.g. GSE242425_h1_d20_sample_filtered_feature_bc_matrix.h5
        parts = base.replace("GSE242425_", "").split("_")
        line, day = parts[0], parts[1]           # h1 / d20
        cond = "hypoxia" if line.startswith("h") else "normoxia"
        phase = "injury_d20" if day == "d20" else "repair_d25"
        a = sc.read_10x_h5(fp)
        a.var_names_make_unique()
        a.obs["sample"] = f"{line}_{day}"
        a.obs["line"] = line
        a.obs["day"] = day
        a.obs["condition"] = cond
        a.obs["phase"] = phase
        a.obs["dataset"] = "GSE242425"
        adatas.append(a)
        log(f"  {line}_{day} ({cond}/{phase}): {a.n_obs} cells x {a.n_vars} genes")
    adata = sc.concat(adatas, join="outer", fill_value=0, index_unique="-")
    adata.var_names_make_unique()
    return adata


def load_GSE209781(P):
    """Unpack GSE209781_RAW.tar and load the six 10x samples."""
    tar_path = os.path.join(RAW, "GSE209781_RAW.tar")
    if not os.path.exists(tar_path):
        log(f"FATAL: {tar_path} not found. Run scripts/utils/download_geo.sh first.")
        sys.exit(1)
    exdir = os.path.join(RAW, "GSE209781")
    os.makedirs(exdir, exist_ok=True)
    need_extract = not any(f.startswith("GSM") for f in os.listdir(exdir))
    if need_extract:
        log("extracting GSE209781_RAW.tar ...")
        with tarfile.open(tar_path) as tf:
            tf.extractall(exdir)
        for fn in os.listdir(exdir):
            if fn.endswith(".tar.gz"):
                with tarfile.open(os.path.join(exdir, fn)) as tf:
                    tf.extractall(exdir)
    SAMPLES = {
        "GSM6395358": ("NM01", "normal"), "GSM6395359": ("NM02", "normal"),
        "GSM6395360": ("NM03", "normal"), "GSM6395361": ("DKD01", "DKD"),
        "GSM6395362": ("DKD02", "DKD"), "GSM6395363": ("DKD03", "DKD"),
    }
    adatas = []
    for gsm, (sid, cond) in SAMPLES.items():
        # directories are named by sample id (NM01/DKD01...), not GSM accession
        cand = [d for d in os.listdir(exdir)
                if os.path.isdir(os.path.join(exdir, d)) and d == sid]
        if not cand:
            log(f"  WARNING: no dir for {sid}")
            continue
        sdir = os.path.join(exdir, cand[0])
        mtx = None
        for root, _, fs in os.walk(sdir):
            if "matrix.mtx.gz" in fs:
                mtx = root
                break
        if mtx is None:
            log(f"  WARNING: no matrix for {gsm}")
            continue
        a = sc.read_10x_mtx(mtx, var_names="gene_symbols", cache=False)
        # GEO provides the UNFILTERED 10x matrix (all barcodes incl. empty droplets):
        # apply per-sample basic filtering before concatenation (memory control).
        a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False,
                                   inplace=True)
        n0 = a.n_obs
        sc.pp.filter_cells(a, min_genes=P["min_genes_per_cell"])
        sc.pp.filter_cells(a, max_counts=P["max_counts_per_cell"])
        a = a[a.obs["pct_counts_mt"] < P["max_mt_pct"]].copy()
        log(f"  {sid}: {n0} barcodes -> {a.n_obs} cells after per-sample QC")
        a.obs_names = [f"{sid}_{bc}" for bc in a.obs_names]
        a.obs["sample"] = sid
        a.obs["condition"] = cond
        a.obs["dataset"] = "GSE209781"
        a.var_names_make_unique()
        adatas.append(a)
        log(f"  {sid} ({cond}): {a.n_obs} cells x {a.n_vars} genes")
    adata = sc.concat(adatas, join="outer", fill_value=0, index_unique="-")
    adata.var_names_make_unique()
    return adata


def qc_and_filter(adata, P):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    qc = adata.obs.groupby("sample").agg(
        cells=("n_genes_by_counts", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_mt_pct=("pct_counts_mt", "median"),
    )
    tag = adata.obs["dataset"].iloc[0]
    qc.to_csv(os.path.join(TAB, f"A1_qc_summary_{tag}.csv"))
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=P["min_genes_per_cell"])
    sc.pp.filter_cells(adata, max_counts=P["max_counts_per_cell"])
    adata = adata[adata.obs["pct_counts_mt"] < P["max_mt_pct"]].copy()
    sc.pp.filter_genes(adata, min_cells=P["min_cells_per_gene"])
    log(f"QC filter: {n0} -> {adata.n_obs} cells")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    sc.pl.violin(adata, "n_genes_by_counts", groupby="sample", show=False, ax=axes[0])
    axes[0].set_title("genes/cell")
    sc.pl.violin(adata, "total_counts", groupby="sample", show=False, ax=axes[1])
    axes[1].set_title("counts/cell")
    sc.pl.violin(adata, "pct_counts_mt", groupby="sample", show=False, ax=axes[2])
    axes[2].set_title("mt%")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG, f"A1_qc_violin_{tag}.png"), bbox_inches="tight")
    plt.close(fig)
    return adata


def cluster_and_integrate(adata, P, tag):
    """Normalise -> HVG -> batch correction -> PCA -> neighbours -> UMAP -> leiden.

    Batch correction is applied to the scaled expression matrix BEFORE PCA, so the
    embedding used for neighbours/UMAP comes from the corrected data. A
    condition-mixing diagnostic (average silhouette width of the condition label on
    the UMAP graph) is logged: |ASW| close to 0 = well mixed.
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=P["n_top_hvg"], batch_key="sample")
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    # With n=3 donors per condition, each donor is its own sequencing batch, so
    # batch is perfectly collinear with condition. Batch-correction methods
    # (ComBat/Harmony) would then remove the biological signal itself, so no
    # cross-sample correction is applied: clustering is descriptive, and ALL
    # condition inference is done at the donor
    # (pseudobulk) level. Logged explicitly so no silent failure can recur.
    rep_key = P.get("replicate_key", "sample")
    cond_per_rep = adata.obs.groupby(rep_key)["condition"].nunique()
    n_rep_per_cond = adata.obs.groupby("condition")[rep_key].nunique()
    confounded = bool((cond_per_rep == 1).all() and (n_rep_per_cond < 5).all())
    log(f"replicate structure: {n_rep_per_cond.to_dict()} (key={rep_key}); "
        f"batch collinear with condition: {confounded}")
    if confounded:
        sc.tl.pca(adata, n_comps=P["n_pcs"], svd_solver="arpack")
        rep = "X_pca"
        used = ("none (batch collinear with condition at n<5/group; inference at donor "
                "level via pseudobulk)")
    else:
        try:
            import harmonypy  # noqa
            sc.tl.pca(adata, n_comps=P["n_pcs"], svd_solver="arpack")
            sc.external.pp.harmony_integrate(adata, key="sample", max_iter_harmony=30)
            rep = "X_pca_harmony"
            used = "Harmony"
        except Exception:
            sc.pp.combat(adata, key="sample")      # BEFORE PCA (correct order)
            sc.tl.pca(adata, n_comps=P["n_pcs"], svd_solver="arpack")
            rep = "X_pca"
            used = "ComBat(before PCA)"
    log(f"batch correction: {used}")
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=P["neighbors_k"])
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=P["leiden_resolution"], key_added="leiden",
                 flavor="igraph", n_iterations=2)
    # Batch-mixing diagnostic (scIB standard): condition average silhouette width
    # on the PCA embedding with Euclidean distances (Luecken et al., scIB,
    # Nat Methods 2022).
    try:
        from sklearn.metrics import silhouette_score
        asw = silhouette_score(adata.obsm["X_pca"],
                               adata.obs["condition"].astype(str).values,
                               metric="euclidean", sample_size=min(5000, adata.n_obs),
                               random_state=42)
        log(f"condition ASW on PCA embedding: {asw:.3f} "
            f"(0 = perfectly mixed by condition; >0.1 = conditions separated in the "
            f"embedding. With n=3 donors per condition and NO batch correction, a "
            f"positive ASW is EXPECTED and is the reason all condition inference is "
            f"done at the donor/pseudobulk level rather than on embeddings)")
    except Exception as e:
        log(f"ASW diagnostic failed: {e}")
    adata.write(os.path.join(PROC, f"{tag}_integrated.h5ad"))
    for color in ["leiden", "sample", "condition"]:
        if color in adata.obs.columns:
            sc.pl.umap(adata, color=color, show=False, title=f"{tag} - {color}")
            plt.savefig(os.path.join(FIG, f"A2_umap_{tag}_{color}.png"), bbox_inches="tight")
            plt.close("all")
    return adata



def score_lineages(adata, P, tag):
    panels = {k: [g for g in v if g in adata.raw.var_names] for k, v in LINEAGE_MARKERS.items()}
    panels = {k: v for k, v in panels.items() if v}
    score_mat = {}
    for name, genes in panels.items():
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{name}", use_raw=True)
        score_mat[name] = adata.obs[f"score_{name}"].values
    S = pd.DataFrame(score_mat, index=adata.obs_names)
    # thresholded assignment: require a positive score AND a clear margin over the
    # runner-up lineage; otherwise 'Unassigned' (guards against biologically
    # implausible low-confidence calls, e.g. 'Nephron progenitor' in adult kidney).
    sorted_S = np.sort(S.values, axis=1)
    labels = S.idxmax(axis=1).values
    top1 = S.values.max(axis=1)
    margin = top1 - (sorted_S[:, -2] if sorted_S.shape[1] > 1 else 0.0)
    assigned = (top1 > 0.0) & (margin > P.get("lineage_margin", 0.02))
    adata.obs["lineage"] = np.where(assigned, S.idxmax(axis=1).values, "Unassigned")
    log(f"lineage assignment: {int(assigned.sum())}/{len(assigned)} cells assigned "
        f"(threshold: score>0, margin>{P.get('lineage_margin', 0.02)})")
    pd.DataFrame({k: ",".join(v) for k, v in panels.items()}, index=["genes_used"]).T \
        .to_csv(os.path.join(TAB, f"A3_lineage_panels_{tag}.csv"))
    comp = pd.crosstab(adata.obs["lineage"], adata.obs["condition"], normalize="columns") * 100
    comp.to_csv(os.path.join(TAB, f"A4_lineage_composition_{tag}.csv"))
    fig, ax = plt.subplots(figsize=(7, 5))
    comp.plot(kind="bar", ax=ax)
    ax.set_ylabel("cells (%)")
    ax.set_title(f"Lineage composition by condition ({tag})")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG, f"A3_lineage_composition_{tag}.png"), bbox_inches="tight")
    plt.close("all")
    sc.pl.umap(adata, color="lineage", show=False, title=f"{tag} cell lineages")
    plt.savefig(os.path.join(FIG, f"A4_umap_lineage_{tag}.png"), bbox_inches="tight")
    plt.close("all")
    return adata


def _norm_deg_cols(df):
    """scanpy>=1.10 renamed 'log2fc'->'logfoldchange(s)', 'pvals_adj'->'padj'."""
    df = df.copy()
    for old in ("logfoldchanges", "logfoldchange"):
        if old in df.columns and "log2fc" not in df.columns:
            df = df.rename(columns={old: "log2fc"})
    if "padj" in df.columns and "pvals_adj" not in df.columns:
        df = df.rename(columns={"padj": "pvals_adj"})
    return df


def pseudobulk_degs(adata, P, tag):
    """Pseudobulk differential expression (fixes pseudoreplication).

    Cell-level Wilcoxon tests treat cells as independent replicates, which inflates
    significance with n=3 vs 3 donors (the biological replicate). Here we aggregate
    RAW counts per (sample x lineage) into pseudobulk profiles, library-size
    normalise (CPM), log2-transform, and test each gene with a Welch t-test on the
    pseudobulk values (n = number of donors per condition), BH-adjusted. This is the
    accepted pseudobulk paradigm for few-replicate scRNA-seq designs.
    """
    from scipy.stats import ttest_ind
    raw = adata.raw.to_adata()
    case, ctrl = P["contrast"]
    rows = []
    rep_key = P.get("replicate_key", "sample")
    for lin in [l for l in adata.obs["lineage"].unique() if l != "Unassigned"]:
        sub = adata.obs["lineage"] == lin
        reps = adata.obs.loc[sub, rep_key].astype(str)
        conds = adata.obs.loc[sub, "condition"].astype(str)
        # pseudobulk unit = biological replicate (donor / iPSC line), NOT sample:
        # samples sharing a replicate (e.g. d20/d25 of one line) are aggregated first.
        case_s = sorted(set(reps[conds == case]))
        ctrl_s = sorted(set(reps[conds == ctrl]))
        if len(case_s) < 3 or len(ctrl_s) < 3:
            continue
        cnt = pd.DataFrame(raw.X[sub.values].toarray(),
                           index=adata.obs_names[sub.values],
                           columns=raw.var_names).astype(float)
        pb = cnt.groupby(reps.values).sum().T                       # genes x replicates
        pb = pb.reindex(columns=case_s + ctrl_s).fillna(0)
        # expression filter: CPM > 1 in >= half of the samples of either condition
        cpm = pb.div(pb.sum(axis=0), axis=1) * 1e6
        keep = ((cpm[case_s] > 1).sum(axis=1) >= max(1, len(case_s) // 2)) | \
               ((cpm[ctrl_s] > 1).sum(axis=1) >= max(1, len(ctrl_s) // 2))
        pb, cpm = pb.loc[keep], cpm.loc[keep]
        logcpm = np.log2(cpm + 1)
        a = logcpm[case_s].values
        b = logcpm[ctrl_s].values
        with np.errstate(divide="ignore", invalid="ignore"):
            t, p = ttest_ind(a, b, axis=1, equal_var=False)
        l2fc = a.mean(axis=1) - b.mean(axis=1)
        df = pd.DataFrame({"names": logcpm.index, "log2fc": l2fc,
                           "pvals": np.nan_to_num(p, nan=1.0), "lineage": lin,
                           "n_samples_case": len(case_s), "n_samples_ctrl": len(ctrl_s)})
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    degs = pd.concat(rows, ignore_index=True)
    # BH adjustment within lineage
    from scipy.stats import false_discovery_control
    degs["pvals_adj"] = np.nan
    for lin, idx in degs.groupby("lineage").groups.items():
        degs.loc[idx, "pvals_adj"] = false_discovery_control(
            degs.loc[idx, "pvals"].values, method="bh")
    degs.to_csv(os.path.join(TAB, f"A6_condition_DEGs_per_lineage_{tag}.csv"), index=False)
    sig = degs[(degs["pvals"] < 0.05) & (degs["log2fc"].abs() >= P.get("min_log2fc", 0.5))]
    sig_fdr = degs[(degs["pvals_adj"] < 0.05) & (degs["log2fc"].abs() >= P.get("min_log2fc", 0.5))]
    log(f"pseudobulk DEGs ({P['contrast_label']}): {len(sig)} nominal (p<0.05, "
        f"|log2FC|>={P.get('min_log2fc', 0.5)}) of {len(degs)} tests; "
        f"{len(sig_fdr)} survive BH-FDR correction (expected: with n=3 donors per "
        f"condition, donor-level FDR control is extremely stringent - reported "
        f"transparently; downstream claims do not rely on single-DEG significance)")
    return degs


def markers_states_degs(adata, P, tag):
    # cluster markers (descriptive only, not used for inference)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True,
                            key_added="markers_leiden")
    mk = _norm_deg_cols(sc.get.rank_genes_groups_df(adata, group=None, key="markers_leiden"))
    mk = mk[(mk["log2fc"].abs() >= P.get("min_log2fc", 0.5)) & (mk["pvals_adj"] < 0.05)]
    mk.to_csv(os.path.join(TAB, f"A5_cluster_markers_{tag}.csv"), index=False)
    # condition DEGs per lineage via PSEUDOBULK (biological replicates = donors)
    degs = pseudobulk_degs(adata, P, tag)
    # state scores: per-SAMPLE scores (pseudobulk-based), Welch t-test on n donors
    case, ctrl = P["contrast"]
    raw = adata.raw.to_adata()
    state_rows = []
    for panel_name, panel in [("injury", INJURY_PANEL),
                              ("maladaptive_repair", MALADAPTIVE_REPAIR_PANEL),
                              ("fibroblast_activation", FIB_ACTIVATION_PANEL)]:
        genes = [g for g in panel if g in raw.var_names]
        # cell-level score for visualisation only
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"state_{panel_name}", use_raw=True)
        log(f"  state panel {panel_name}: {len(genes)}/{len(panel)} genes present")
        rep_key = P.get("replicate_key", "sample")
        for lin in [l for l in adata.obs["lineage"].unique() if l != "Unassigned"]:
            sub = adata.obs["lineage"] == lin
            reps = adata.obs.loc[sub, rep_key].astype(str)
            conds = adata.obs.loc[sub, "condition"].astype(str)
            cnt = pd.DataFrame(raw.X[sub.values].toarray(),
                               index=adata.obs_names[sub.values],
                               columns=raw.var_names).astype(float)
            pb = cnt.groupby(reps.values).sum().T
            cpm = pb.div(pb.sum(axis=0), axis=1) * 1e6
            logcpm = np.log2(cpm.loc[genes] + 1)
            z = logcpm.sub(logcpm.mean(axis=1), axis=0).div(
                logcpm.std(axis=1).replace(0, 1), axis=0)
            per_sample = z.mean(axis=0)
            case_s = [s for s in per_sample.index if conds[reps == s].iloc[0] == case]
            ctrl_s = [s for s in per_sample.index if conds[reps == s].iloc[0] == ctrl]
            if len(case_s) < 3 or len(ctrl_s) < 3:
                continue
            a, b = per_sample[case_s].values, per_sample[ctrl_s].values
            from scipy.stats import ttest_ind
            t, p = ttest_ind(a, b, equal_var=False)
            state_rows.append({"panel": panel_name, "lineage": lin,
                               "mean_case": float(np.mean(a)), "mean_ctrl": float(np.mean(b)),
                               "welch_t": float(t), "welch_p": float(p),
                               "n_samples_case": len(case_s), "n_samples_ctrl": len(ctrl_s),
                               "per_sample_scores": ";".join(
                                   f"{s}:{v:.3f}" for s, v in per_sample.items()),
                               # condition map so downstream figures can place each
                               # individual replicate under its TRUE condition
                               "condition_map": ";".join(
                                   f"{s}:{conds[reps == s].iloc[0]}" for s in per_sample.index)})
    pd.DataFrame(state_rows).to_csv(os.path.join(TAB, f"A7_state_scores_{tag}.csv"), index=False)
    # figures: cell-level UMAP (visualisation) + per-sample violin
    for col in ["state_injury", "state_maladaptive_repair", "state_fibroblast_activation"]:
        sc.pl.umap(adata, color=col, show=False, title=f"{tag} {col}", cmap="magma")
        plt.savefig(os.path.join(FIG, f"A5_umap_{tag}_{col}.png"), bbox_inches="tight")
        plt.close("all")
        srows = pd.read_csv(os.path.join(TAB, f"A7_state_scores_" + tag + ".csv"))
        # map this figure's column back to its state panel name
        panel_of_col = {"state_injury": "injury",
                        "state_maladaptive_repair": "maladaptive_repair",
                        "state_fibroblast_activation": "fibroblast_activation"}
        this_panel = panel_of_col.get(col)
        # replicate -> condition map rebuilt from the object (no leaked variables)
        rep_key = P.get("replicate_key", "sample")
        cond_map = adata.obs.groupby(rep_key)["condition"].first().to_dict()
        case, ctrl = P["contrast"]
        sub = srows[(srows["panel"] == this_panel)].copy()
        sub = sub.sort_values("welch_p").head(3)
        if len(sub):
            fig, axes = plt.subplots(1, len(sub), figsize=(3.4 * len(sub), 3.4))
            axes = np.atleast_1d(axes)
            for ax, (_, r) in zip(axes, sub.iterrows()):
                per_rep = {}
                for pair in str(r["per_sample_scores"]).split(";"):
                    s, v = pair.split(":")
                    per_rep[s] = float(v)
                for i, cond in enumerate([ctrl, case]):
                    vals = [v for s, v in per_rep.items() if cond_map.get(s) == cond]
                    if vals:
                        ax.bar(i, np.mean(vals), yerr=np.std(vals), capsize=4,
                               color=C_CTRL if i == 0 else C_CASE, alpha=0.75,
                               width=0.5)
                        ax.scatter(np.full(len(vals), i) +
                                   np.random.uniform(-0.09, 0.09, len(vals)),
                                   vals, color="black", s=16, zorder=5)
                        for v in vals:
                            ax.text(i, v + 0.04, format(v, ".2f"), ha="center",
                                    fontsize=6)
                ax.set_xticks([0, 1])
                ax.set_xticklabels([ctrl, case], fontsize=8)
                ax.set_ylabel("state score (z)", fontsize=8)
                ax.set_title(r["lineage"] + " (Welch p=" +
                             format(r["welch_p"], ".3g") + ")", fontsize=8)
            plt.suptitle(tag + ": " + col + " (per-replicate pseudobulk, "
                         "mean +/- s.d., n = " +
                         str(sum(1 for s in per_rep if cond_map.get(s) == case)) +
                         " vs " +
                         str(sum(1 for s in per_rep if cond_map.get(s) == ctrl)) +
                         " replicates)", fontsize=8)
            plt.tight_layout(rect=[0, 0, 1, 0.93])
            fig.savefig(os.path.join(FIG, "A6_violin_" + tag + "_" + col + ".png"),
                        bbox_inches="tight")
            plt.close("all")
    return adata, degs



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=["GSE242425", "GSE209781"])
    args = ap.parse_args()
    tag = args.dataset
    P = dict(PARAMS[tag])
    P["run_timestamp"] = datetime.now().isoformat(timespec="seconds")
    P["scanpy_version"] = "1.12.4"
    log(f"Module A start | dataset={tag} | contrast={P['contrast_label']}")
    integrated = os.path.join(PROC, f"{tag}_integrated.h5ad")
    if os.path.exists(integrated):
        log(f"resuming from existing {integrated}")
        adata = sc.read(integrated)
        if "leiden" not in adata.obs.columns:
            log("integrated file lacks clustering -> reprocessing from scratch")
            adata = None
    else:
        adata = None
    if adata is None:
        adata = load_GSE242425() if tag == "GSE242425" else load_GSE209781(P)
        log(f"loaded {adata.n_obs} cells x {adata.n_vars} genes")
        adata = qc_and_filter(adata, P)
        adata = cluster_and_integrate(adata, P, tag)
    adata = score_lineages(adata, P, tag)
    adata, degs = markers_states_degs(adata, P, tag)
    adata.write(os.path.join(PROC, f"{tag}_annotated.h5ad"))
    with open(os.path.join(PROC, f"params_moduleA_{tag}.json"), "w", encoding="utf-8") as f:
        json.dump(P, f, ensure_ascii=False, indent=1)
    log(f"Module A done for {tag}: cells={adata.n_obs}, lineages={adata.obs['lineage'].nunique()}, "
        f"clusters={adata.obs['leiden'].nunique()}")
    log(f"outputs -> {PROC} | {TAB} | {FIG}")


if __name__ == "__main__":
    main()
