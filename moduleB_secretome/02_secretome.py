# -*- coding: utf-8 -*-
"""
Module B: state-restricted secretome extraction.
==================================================================================
Input : Module A outputs (condition DEGs per lineage + state panels), public
        annotation resources:
        - Human Protein Atlas (HPA) predicted secretome (API, CC BY)
        - UniProt reviewed human "Secreted" subcellular-location annotation (REST API)
        - eQTLGen measured gene list (from the public cis-eQTL significant file)
Output: candidate secreted-protein panel with cell-state provenance + ranking,
        full audit trail of every filtering step (numbers traceable).

Run:  python 02_secretome.py
"""
import os
import io
import numpy as np
import gzip
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
PROC = os.path.join(ROOT, "data", "processed")
EXT = os.path.join(ROOT, "data", "external")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
for d in (EXT, TAB):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------- helpers
def http_get(url, timeout=120, retries=3):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------- 1. state genes from Module A
def collect_state_genes():
    """Union of (a) condition-DEGs in injury/fibroblast lineages, (b) state panel genes,
    across all processed Module A datasets."""
    deg_files = sorted(f for f in os.listdir(TAB) if f.startswith("A6_condition_DEGs"))
    state_files = sorted(f for f in os.listdir(TAB) if f.startswith("A7_state_scores"))
    if not deg_files:
        log("FATAL: no Module A DEG output. Run 01_scRNA_cellstates.py first.")
        raise SystemExit(1)
    frames, provenance = [], []
    for f in deg_files:
        df = pd.read_csv(os.path.join(TAB, f))
        tag = f.replace("A6_condition_DEGs_per_lineage_", "").replace(".csv", "")
        df["dataset"] = tag
        # donor-level (pseudobulk) DEGs. note: with n=3 donors per condition no gene
        # survives BH-FDR correction (reported transparently in STATUS/manuscript); the
        # state-gene set therefore uses NOMINAL p<0.05 with |log2FC|>=0.5, and all
        # downstream claims rest on the independent validation layers (bulk cohorts,
        # MR replication), never on single-DEG significance.
        sig = df[(df["pvals"] < 0.05) & (df["log2fc"].abs() >= 0.5)]
        # keep DEGs from injury-relevant lineages (epithelial/fibroblast/myeloid)
        keep_lin = [l for l in sig["lineage"].unique()
                    if any(k in l.lower() for k in
                           ["tubule", "fibro", "macro", "mono", "endothel", "podocyte",
                            "collecting", "loop", "distal"])]
        sig = sig[sig["lineage"].isin(keep_lin)]
        frames.append(sig)
        provenance.append({"dataset": tag, "file": f,
                           "n_deg_lineages_kept": sig["lineage"].nunique(),
                           "n_sig_degs": len(sig)})
    degs = pd.concat(frames, ignore_index=True)
    # genes significant in >=1 dataset
    gene_counts = degs.groupby("names")["dataset"].nunique().rename("n_datasets")
    gene_dir = degs.groupby("names")["log2fc"].mean().rename("mean_log2fc")
    # direction tracking (review point): the state set contains both up- and
    # down-regulated genes; record direction and cross-dataset consistency so the
    # panel composition is transparent instead of silently averaged.
    up = degs[degs["log2fc"] > 0].groupby("names")["dataset"].nunique().rename("n_up")
    down = degs[degs["log2fc"] < 0].groupby("names")["dataset"].nunique().rename("n_down")
    state_genes = pd.concat([gene_counts, gene_dir, up, down], axis=1).fillna(0)
    state_genes["direction"] = np.where(
        (state_genes["n_up"] > 0) & (state_genes["n_down"] > 0), "mixed",
        np.where(state_genes["n_up"] > 0, "up", "down"))
    state_genes["direction_consistent"] = (
        (state_genes["n_up"] == 0) | (state_genes["n_down"] == 0))
    state_genes = state_genes.reset_index()
    state_genes.columns = ["gene", "n_datasets", "mean_log2fc", "n_up", "n_down",
                           "direction", "direction_consistent"]
    n_up_tot = int((degs["log2fc"] > 0).sum())
    n_dn_tot = int((degs["log2fc"] < 0).sum())
    log(f"state-associated DEGs (injury-relevant lineages): {len(state_genes)} genes "
        f"across {degs['dataset'].nunique()} datasets; direction: {n_up_tot} up / "
        f"{n_dn_tot} down DEG entries")
    return state_genes, provenance


# ---------------------------------------------------------------- 2. HPA secretome 
def fetch_hpa_secretome():
    """HPA 'secreted' search -> gene list (predicted secreted proteins)."""
    url = ("https://www.proteinatlas.org/api/search_download.php?search=secreted"
           "&format=json&columns=g&compress=no")
    raw = http_get(url, timeout=180)
    data = json.loads(raw.decode("utf-8"))
    genes = sorted({d["Gene"] for d in data if d.get("Gene")})
    with open(os.path.join(EXT, "hpa_secretome.json"), "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    log(f"HPA predicted secretome: {len(genes)} genes (source: proteinatlas.org API)")
    return genes


# ---------------------------------------------------------------- 3. UniProt secreted 
def fetch_uniprot_secreted():
    """UniProt reviewed human proteins annotated 'Secreted', streamed TSV."""
    query = ('organism_id:9606 AND reviewed:true AND '
             'cc_subcellular_location:"Secreted"')
    url = ("https://rest.uniprot.org/uniprotkb/stream?fields=gene_primary&format=tsv&query="
           + urllib.parse.quote(query))
    raw = http_get(url, timeout=300)
    lines = raw.decode("utf-8").strip().split("\n")
    genes = sorted({l.split("\t")[0].strip() for l in lines[1:] if l.strip()})
    with open(os.path.join(EXT, "uniprot_secreted.tsv"), "w", encoding="utf-8") as f:
        f.write("\n".join(genes))
    log(f"UniProt reviewed human 'Secreted' proteins: {len(genes)} genes (source: UniProt REST)")
    return genes


# ---------------------------------------------------------------- 4. eQTLGen measured genes
EQTLGEN_URL = ("https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
               "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
# note: same local cache path as Module D (04_MR_causal.py) - downloaded once, reused.
EQTLGEN_LOCAL = os.path.join(ROOT, "data", "raw", "eqtlgen_cis_FDR05.txt.gz")


def fetch_eqtlgen_genes():
    """Gene list from the public eQTLGen significant cis-eQTL file (real download)."""
    gene_file = os.path.join(EXT, "eqtlgen_genes.txt")
    if os.path.exists(gene_file):
        genes = [g.strip() for g in open(gene_file) if g.strip()]
        log(f"eQTLGen measured genes (cached): {len(genes)}")
        return genes
    if not os.path.exists(EQTLGEN_LOCAL):
        log("downloading eQTLGen significant cis-eQTL file (323 MB) ...")
        http_get(EQTLGEN_URL, timeout=7200)  # verify reachable first
        import urllib.request
        req = urllib.request.Request(EQTLGEN_URL, headers={"User-Agent": "Mozilla/5.0 (research)"})
        with urllib.request.urlopen(req, timeout=7200) as r, open(EQTLGEN_LOCAL, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
    log(f"  cached at {EQTLGEN_LOCAL} ({os.path.getsize(EQTLGEN_LOCAL)/1e6:.0f} MB)")
    genes = set()
    with gzip.open(EQTLGEN_LOCAL, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
        # note: header contains BOTH 'Gene' (Ensembl ID) and 'GeneSymbol' - we want the symbol
        lower = [h.strip().lower() for h in header]
        if "genesymbol" in lower:
            sym_col = lower.index("genesymbol")
        elif "gene_symbol" in lower:
            sym_col = lower.index("gene_symbol")
        elif "gene" in lower:
            sym_col = lower.index("gene")
        else:
            log(f"FATAL: cannot find gene column in header: {header}")
            raise SystemExit(1)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > sym_col:
                genes.add(parts[sym_col].strip())
    genes = sorted(g for g in genes if g)
    with open(gene_file, "w", encoding="utf-8") as f:
        f.write("\n".join(genes))
    log(f"eQTLGen measured genes: {len(genes)} (source: eqtlgen.org)")
    return genes


# ---------------------------------------------------------------- main
def main():
    state_genes, prov = collect_state_genes()
    pd.DataFrame(prov).to_csv(os.path.join(TAB, "B0_state_gene_provenance.csv"), index=False)

    hpa = fetch_hpa_secretome()
    uni = fetch_uniprot_secreted()
    eqt = fetch_eqtlgen_genes()

    df = state_genes.copy()
    df["in_HPA_secretome"] = df["gene"].isin(hpa)
    df["in_UniProt_secreted"] = df["gene"].isin(uni)
    df["in_eQTLGen"] = df["gene"].isin(eqt)
    # consensus secretome: HPA predicted OR UniProt annotated
    df["secreted_consensus"] = df["in_HPA_secretome"] | df["in_UniProt_secreted"]
    panel = df[df["secreted_consensus"] & df["in_eQTLGen"]].copy()
    panel = panel.sort_values(["n_datasets", "mean_log2fc"],
                              ascending=[False, False]).reset_index(drop=True)
    panel.insert(0, "rank", panel.index + 1)

    # audit trail (every number traceable)
    audit = {
        "state_genes_total": int(len(df)),
        "in_HPA_secretome": int(df["in_HPA_secretome"].sum()),
        "in_UniProt_secreted": int(df["in_UniProt_secreted"].sum()),
        "secreted_consensus": int(df["secreted_consensus"].sum()),
        "in_eQTLGen": int(df["in_eQTLGen"].sum()),
        "final_panel": int(len(panel)),
        "hpa_n": len(hpa), "uniprot_n": len(uni), "eqtlgen_n": len(eqt),
        "run_timestamp": datetime.now().isoformat(timespec="seconds"),
        "sources": {
            "HPA": "https://www.proteinatlas.org/api/search_download.php?search=secreted",
            "UniProt": 'organism_id:9606 AND reviewed:true AND cc_subcellular_location:"Secreted"',
            "eQTLGen": "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
        },
    }
    with open(os.path.join(PROC, "moduleB_audit.json"), "w", encoding="utf-8") as f:
        json.dump(audit, f, ensure_ascii=False, indent=1)
    panel.to_csv(os.path.join(TAB, "B1_secreted_candidate_panel.csv"), index=False)
    log(f"FILTER FUNNEL: state genes {audit['state_genes_total']} -> secretome "
        f"{audit['secreted_consensus']} -> +eQTLGen instruments {audit['final_panel']}")
    log(f"panel -> {os.path.join(TAB, 'B1_secreted_candidate_panel.csv')}")
    if len(panel):
        log("top 15 candidates:\n" + panel.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
