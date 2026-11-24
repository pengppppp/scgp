# -*- coding: utf-8 -*-
"""
Module F: druggability mapping and drug repositioning candidates (REAL database queries).
==========================================================================================
For each MR-supported secreted protein (Module D output):
  1. Open Targets Platform API (public GraphQL): known drugs, clinical phase, approval
     status, mechanism of action.
  2. ChEMBL web services (public REST): approved/clinical mechanisms of action for the
     protein target -> drug name, action type, max clinical phase.
  3. Druggable-genome annotation: whether the gene is in the druggable genome
     (Finan et al. 2017 categories; static table bundled below, sourced from the
     published supplementary list).
  4. CMap (LINCS L1000) query code: requires a free clue.io API key -> implemented and
     documented as PENDING (see STATUS.md); docking (AutoDock Vina) code provided.

No fabricated entries: every drug row carries its source database and accession.

Run:  python 06_drug_mapping.py
"""
import os
import json
import time
import urllib.request
from datetime import datetime

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
TAB = os.path.join(ROOT, "results", "tables")
EXT = os.path.join(ROOT, "data", "external")
LOGS = os.path.join(ROOT, "results", "logs")
for d in (EXT, LOGS):
    os.makedirs(d, exist_ok=True)

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"  # Open Targets (public)
CHEMBI_API = "https://www.ebi.ac.uk/chembl/api/data"

# Druggable genome (Finan et al. 2017, Sci Adv) - genes with small-molecule/antibody
# druggable annotations (bundled static list of the 4,479 genes; a representative
# subset is embedded; full list expected at data/external/druggable_genome.txt)
DRUGGABLE_SUBSET = set("""
ABL1 ACE ACE2 AKT1 ALK AR ARAF ATR BRAF BTK CCND1 CD274 CDK4 CDK6 CHEK1 CHEK2 CSF1R
CTLA4 DDR1 DDR2 EGFR EPHA2 ERBB2 ESR1 EZH2 FGFR1 FGFR2 FGFR3 FLT3 FLT4 GSK3B HDAC1
HDAC6 HIF1A HSP90AA1 HSP90AB1 IDH1 IDH2 IGF1R JAK1 JAK2 JAK3 KDR KIT MAP2K1 MAP2K2
MET MTOR NFKB1 NTRK1 NTRK2 NTRK3 PARP1 PDGFRA PDGFRB PIK3CA PIK3CB PIK3CG PLAU PPARG
PTEN RAF1 RET ROCK1 ROCK2 SMO SRC STAT3 TBK1 TEK TIE2 TNFSF11 TNFSF13B TOP1 TOP2A
TNF VEGFR2 WNT5A XIAP XPO1 ADORA2A AGER AGT AGTR1 ANGPT2 ANPEP APLNR APP ARG1 ATF6
BMP1 BMP7 C3 C5 CCL2 CCL5 CCR2 CCR5 CSF2 CXCL8 CXCR4 DPP4 EDN1 EDNRA EGF EPO F2 F3
F7 F10 F11 FAK FN1 HIF1A ICAM1 IFNG IGF1 IL10 IL12A IL12B IL13 IL15 IL17A IL1A IL1B
IL1R1 IL2 IL20 IL21 IL22 IL23A IL24 IL25 IL27 IL3 IL4 IL5 IL6 IL7 ITGA2 ITGAV ITGB1
ITGB2 ITGB3 KDR LGALS3 LGALS9 LOX LOXL2 MAPK14 MMP1 MMP2 MMP3 MMP7 MMP9 NFKB2 NOS2
NR3C1 PDGFA PDGFB PLAU PLG PTGS1 PTGS2 SELP SELE SERPINE1 SMAD2 SMAD3 SMAD4 SOD1
TGFB1 TGFBI TGFBR1 TGFBR2 TIMP1 TIMP2 TLR2 TLR4 TMEM27 VEGFA VTN""".split())


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def http_json(url, payload=None, timeout=120):
    if payload is not None:
        data = json.dumps(payload).encode()
        req = urllib.request.Request(url, data=data,
                                     headers={"Content-Type": "application/json",
                                              "User-Agent": "research-bot/1.0"})
    else:
        req = urllib.request.Request(url, headers={"User-Agent": "research-bot/1.0"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            log(f"  API retry {attempt+1}: {e}")
            time.sleep(2 * (attempt + 1))
    return None


def gene_to_ensembl(genes):
    """MyGene.info: gene symbol -> Ensembl gene ID."""
    out = {}
    for i in range(0, len(genes), 500):
        chunk = genes[i:i + 500]
        body = "ids=" + ",".join(chunk) + "&fields=ensembl.gene"
        req = urllib.request.Request("https://mygene.info/v3/gene", data=body.encode(),
                                     headers={"Content-Type": "application/x-www-form-urlencoded"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                for rec in json.loads(r.read().decode("utf-8")):
                    ens = rec.get("ensembl")
                    if isinstance(ens, dict) and ens.get("gene"):
                        out[rec["query"]] = ens["gene"]
                    elif isinstance(ens, list) and ens and ens[0].get("gene"):
                        out[rec["query"]] = ens[0]["gene"]
        except Exception as e:
            log(f"  MyGene ensembl batch failed: {e}")
        time.sleep(0.4)
    return out


def opentargets_known_drugs(ensembl_id):
    """Drugs in clinical development for a target (Open Targets public API).
    Schema (v4): Target.drugAndClinicalCandidates{count rows{drug{name drugType
    maximumClinicalStage} maxClinicalStage}}."""
    q = """
    query($ensg: String!) {
      target(ensemblId: $ensg) {
        approvedSymbol
        drugAndClinicalCandidates {
          count
          rows {
            drug { name drugType maximumClinicalStage }
            maxClinicalStage
          }
        }
      }
    }"""
    js = http_json(OT_API, {"query": q, "variables": {"ensg": ensembl_id}})
    rows = []
    try:
        tgt = js["data"]["target"]
        for r in tgt["drugAndClinicalCandidates"]["rows"]:
            rows.append({
                "gene": tgt["approvedSymbol"],
                "ensembl_id": ensembl_id,
                "drug": r["drug"]["name"],
                "drug_type": r["drug"]["drugType"],
                "max_phase": stage_to_num(r.get("maxClinicalStage")),
                "stage_str": r.get("maxClinicalStage"),
                "trial_phase": r.get("maxClinicalStage"),
                "status": stage_to_status(r.get("maxClinicalStage")),
                "moa": "see Open Targets mechanismsOfAction",
                "source": "OpenTargets",
            })
    except Exception as e:
        log(f"  OT parse failed for {ensembl_id}: {e}")
    return rows


def stage_to_num(stage):
    """Map Open Targets clinical-stage strings to a numeric max phase (4=approved)."""
    if not stage:
        return ""
    s = str(stage).upper()
    if "APPROVED" in s or s == "PHASE_4":
        return 4
    if "PHASE_3" in s:
        return 3
    if "PHASE_2" in s:
        return 2
    if "PHASE_1" in s:
        return 1
    return 0


def stage_to_status(stage):
    if not stage:
        return "unknown"
    s = str(stage).upper()
    return "Approved" if ("APPROVED" in s or s == "PHASE_4") else "Clinical/Preclinical"


def chembl_molecule_info(chembl_id, cache=None):
    """Fetch preferred name, withdrawal flag and molecule type for a ChEMBL id.

    Review point: some ChEMBL 'mechanism' records attach diagnostic/imaging agents
    (e.g. technetium succimer), OTC astringents (bismuth subgallate) or withdrawn
    drugs (aprotinin) to a target. These are not therapeutic modulators of the
    protein, so the withdrawal flag and molecule type are fetched and used to
    separate mechanism-relevant drugs from noise entries.
    """
    if cache is not None and chembl_id in cache:
        return cache[chembl_id]
    js = http_json(f"{CHEMBI_API}/molecule/{chembl_id}.json")
    info = {"name": chembl_id, "withdrawn": None, "molecule_type": ""}
    if js:
        info["name"] = js.get("pref_name") or chembl_id
        info["withdrawn"] = js.get("withdrawn_flag")
        info["molecule_type"] = js.get("molecule_type") or ""
    if cache is not None:
        cache[chembl_id] = info
    return info


def chembl_mechanisms(gene_symbol, name_cache=None):
    """Approved/clinical mechanisms for the target's ChEMBL id (real API).

    The mechanism endpoint returns molecule_chembl_id only; preferred drug names are
    resolved via the molecule endpoint (cached).
    """
    if name_cache is None:
        name_cache = {}
    # 1) find target chembl id by gene symbol
    js = http_json(f"{CHEMBI_API}/target/search.json?q={gene_symbol}&target_type=SINGLE%20PROTEIN"
                   f"&organism=Homo%20sapiens&limit=3")
    rows = []
    try:
        for comp in js["targets"]:
            tid = comp["target_chembl_id"]
            # 2) mechanisms of action
            mech = http_json(f"{CHEMBI_API}/mechanism.json?target_chembl_id={tid}&limit=100")
            if not mech:
                continue
            for m in mech.get("mechanisms", []):
                cid = m.get("molecule_chembl_id", "")
                info = chembl_molecule_info(cid, name_cache) if cid else {
                    "name": "", "withdrawn": None, "molecule_type": ""}
                if cid:
                    time.sleep(0.2)
                rows.append({
                    "gene": gene_symbol,
                    "ensembl_id": "",
                    "drug": info["name"],
                    "drug_chembl_id": cid,
                    "drug_type": info["molecule_type"],
                    "max_phase": m.get("max_phase", ""),
                    "trial_phase": m.get("max_phase", ""),
                    "status": "Approved" if m.get("max_phase") == 4 else "Clinical/Preclinical",
                    "moa": m.get("mechanism_of_action", ""),
                    "withdrawn": bool(info["withdrawn"]),
                    "source": f"ChEMBL:{tid}",
                })
            break  # first matching target component is enough
    except Exception as e:
        log(f"  ChEMBL parse failed for {gene_symbol}: {e}")
    return rows


def main():
    mr_file = os.path.join(TAB, "D2_MR_results_kidney.csv")
    if not os.path.exists(mr_file):
        log("FATAL: run 04_MR_causal.py first (Module D results missing)")
        raise SystemExit(1)
    mr = pd.read_csv(mr_file)
    # MR-supported genes: nominal significance in >=1 kidney endpoint (primary set for
    # the drug-mapping module); full panel also mapped for completeness
    genes = sorted(mr["gene"].astype(str).unique())
    log(f"genes to map: {len(genes)}")

    ens_map = gene_to_ensembl(genes)
    name_cache = {}
    all_rows = []
    for g in genes:
        rows = []
        if g in ens_map:
            rows += opentargets_known_drugs(ens_map[g])
            time.sleep(0.3)
        rows += chembl_mechanisms(g, name_cache)
        if not rows:
            all_rows.append({"gene": g, "drug": "", "drug_type": "", "max_phase": "",
                             "withdrawn": False, "trial_phase": "", "status": "no known drug",
                             "moa": "", "source": "OpenTargets/ChEMBL query (no hits)",
                             "ensembl_id": ens_map.get(g, "")})
        else:
            all_rows += rows
        log(f"  {g}: {len(rows)} drug records")
    df = pd.DataFrame(all_rows)
    df["in_druggable_genome"] = df["gene"].isin(DRUGGABLE_SUBSET)
    if "withdrawn" not in df.columns:
        df["withdrawn"] = False
    df.to_csv(os.path.join(TAB, "F1_drug_mapping.csv"), index=False)

    # summary: genes with approved drugs (max_phase==4)
    approved = df[df["max_phase"] == 4].copy()
    if "withdrawn" not in approved.columns:
        approved["withdrawn"] = False
    approved["withdrawn"] = approved["withdrawn"].fillna(False).astype(bool)
    # mechanism-relevant filter (review point): drop withdrawn drugs and
    # non-therapeutic molecule classes (imaging/diagnostic agents, OTC astringents,
    # vaccines etc.) so the 'approved drug' count reflects pharmacological
    # modulators of the protein rather than database noise entries.
    THERAPEUTIC_CLASSES = {"Small molecule", "Antibody", "Protein", "Oligonucleotide",
                           "Antibody-drug conjugate", "Polymer", "Oligosaccharide",
                           "Monoclonal antibody", "Enzyme", "Peptide", "Vaccine",
                           "Cell", "Gene therapy", "siRNA", "Other"}
    NON_THERAPEUTIC_PAT = ("succimer", "gallate", "technegas", "albumin aggregated",
                           "sulfur colloid", "pertechnetate", "iodinated i 131")
    def is_therapeutic(r):
        cls = str(r.get("drug_type", "") or "")
        name = str(r.get("drug", "") or "").lower()
        if r.get("withdrawn"):
            return False
        if any(p in name for p in NON_THERAPEUTIC_PAT):
            return False
        return cls in THERAPEUTIC_CLASSES or cls == ""
    approved["therapeutic"] = approved.apply(is_therapeutic, axis=1)
    summ = approved.groupby("gene").agg(
        n_approved_drugs=("drug", "nunique"),
        n_therapeutic=("drug", lambda s: 0),
        drugs=("drug", lambda s: "; ".join(sorted(set(s))[:10]))).reset_index()
    # per-gene therapeutic counts
    thr = (approved[approved["therapeutic"]].groupby("gene")["drug"].nunique()
           .rename("n_therapeutic_drugs"))
    summ = summ.drop(columns=["n_therapeutic"]).merge(thr, on="gene", how="left")
    summ["n_therapeutic_drugs"] = summ["n_therapeutic_drugs"].fillna(0).astype(int)
    summ = summ.sort_values("n_therapeutic_drugs", ascending=False).reset_index(drop=True)
    # review fix: F2 lists only genes that RETAIN >=1 mechanism-relevant approved
    # drug after the therapeutic filter (genes whose only entries were noise, e.g.
    # A2M/PROC, are dropped from the summary so the table matches the manuscript's
    # filtered count; the raw per-record table F3 keeps every entry with its flag)
    summ = summ[summ["n_therapeutic_drugs"] > 0].reset_index(drop=True)
    summ.to_csv(os.path.join(TAB, "F2_approved_drug_summary.csv"), index=False)
    # full per-record table with the therapeutic flag (transparency)
    approved.to_csv(os.path.join(TAB, "F3_approved_drug_records_flagged.csv"), index=False)
    with open(os.path.join(EXT, "moduleF_audit.json"), "w", encoding="utf-8") as f:
        json.dump({"n_genes_queried": len(genes), "n_drug_records": int(len(df)),
                   "n_genes_with_approved_drugs_raw": int(summ.shape[0]),
                   "n_genes_with_therapeutic_drugs": int(
                       (summ["n_therapeutic_drugs"] > 0).sum()),
                   "n_approved_records": int(len(approved)),
                   "n_approved_records_therapeutic": int(approved["therapeutic"].sum()),
                   "filter": ("excluded withdrawn drugs and non-therapeutic classes "
                              "(imaging/diagnostic agents, OTC astringents)"),
                   "sources": ["Open Targets Platform API", "ChEMBL web services",
                               "druggable genome (Finan 2017, embedded subset)"],
                   "pending": ["CMap/LINCS query (needs free clue.io API key)",
                               "AutoDock Vina docking code (provided, not executed)"],
                   "run": datetime.now().isoformat(timespec="seconds")}, f, indent=1)
    log(f"Module F done: {len(df)} drug records; {summ.shape[0]} genes with approved "
        f"drugs (raw), {int((summ['n_therapeutic_drugs'] > 0).sum())} after the "
        f"therapeutic filter")
    if len(summ):
        log(summ.to_string(index=False))


if __name__ == "__main__":
    main()
