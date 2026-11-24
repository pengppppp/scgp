# SCGP: Single-Cell anchored Genetic Prioritisation

## Quick start

```bash
bash scgp.sh        # 14 steps, fixed seed 42
```

Steps (modules):

| Module | Content |
|---|---|
| `moduleA_scRNA/` | Single-cell QC, lineage annotation, donor-level pseudobulk DEGs, state scores (GSE242425, GSE209781) |
| `moduleB_secretome/` | State-restricted secreted candidate panel (486 proteins), plasma-detectability annotation |
| `moduleC_bulk/` | Bulk validation of state panels (GSE133288/NEPTUNE, GSE180395/CProbe) |
| `moduleD_MR/` | cis-eQTL instrument selection and MR, endpoint correlation, colocalisation, focused pQTL-MR, CKDGen/BBJ/UACR independent replication, panel-wide multi-instrument cis-pQTL MR, cross-layer summary table |
| `moduleE_population/` | Population replication layer |
| `moduleF_drug/` | Drug mapping (Open Targets + ChEMBL) and approved-drug filtering |
| `utils/make_fig9.py` | Figure 9 (panel-wide pQTL-MR forest/replication/coverage) |
| `utils/consolidate_tables.py` | Bundles `results_tables/` CSVs into the S1-S7 Excel workbooks deposited in GigaDB |

Supplementary analyses that use files too large for the code archive
(`make_kuppe_table.py`, `nid2_ckd_binary.py`, `uacr_replication.py`) read
from `data/raw/` after a one-time download; see "Large input files" below.

## Requirements

- bash; Python 3.13 (pandas, numpy, scipy, statsmodels, matplotlib,
  openpyxl); R 4.5.3 (VariantAnnotation). Versions pinned in
  `requirements.txt`.
- All random procedures use fixed seed 42.

## Repository layout

```
.
├── scgp.sh                  one-command reproduction (13 steps)
├── LICENSE                  MIT
├── requirements.txt         pinned Python dependencies
├── moduleA_scRNA/ ... moduleF_drug/
├── utils/                   figure and helper scripts
├── results_tables/          50 result tables (CSV, one per analysis step;
│                             consolidated into S1-S7 Excel workbooks for GigaDB
│                             by utils/consolidate_tables.py)
├── results/logs/            per-step run logs
├── data/raw/                large input files (download separately, see below)
```

## Public data inputs

| Resource | Accession / location |
|---|---|
| Organoid injury model scRNA-seq | GEO GSE242425 |
| DKD snRNA-seq | GEO GSE209781 |
| NEPTUNE | GEO GSE133288 |
| CProbe | GEO GSE180395 |
| eQTLGen phase-1 cis-eQTL | eqtlgen.org/cis-eqtls.html |
| FinnGen R12 | finn.fi |
| deCODE SomaScan cis-pQTL | decode.is (largescaleplasma bucket) |
| CKDGen eGFRcrea | GCST008059 |
| CKDGen eGFRcys / eGFRcrea trans-ethnic / CKD / UACR | CKDGen server (ckdgen.imbi.uni-freiburg.de) |
| Biobank Japan eGFR | GWAS Catalog bbj-a-60 |
| Kuppe 2021 atlas supplement | Supplementary Data 2 (Nature 2021;589:281-286) |
| HPA / UniProt / MyGene.info / Open Targets / ChEMBL | as cited in the manuscript |

## Large input files

Three supplementary analyses need files that are too large to ship in the
code archive. Download once into `data/raw/` and re-run the scripts:

| Script | File to place |
|---|---|
| `moduleD_MR/make_kuppe_table.py` | `data/raw/kuppe/MOESM4_gene_expression.xlsx` (Kuppe 2021 Supplementary Data 2) |
| `moduleD_MR/nid2_ckd_binary.py` | `data/raw/ckdgen/CKD_overall_EA.dbgap.txt.gz` (CKDGen CKD binary, Wuttke 2019) |
| `moduleD_MR/uacr_replication.py` | `data/raw/ckdgen/UACR_overall_ALL.tbl.rsid.gz` (CKDGen UACR, Teumer 2019) |

The CKDGen files are also mirrored in the project's data deposit; the
Kuppe supplement is available from the publisher's article page.

## Reproducibility notes

- A SHA-256 manifest of all scripts is included at the package root.
- Run logs record environment, versions and per-step row counts.
- Manuscript figures are produced by the module scripts plus
  `utils/make_fig9.py`; the manuscript text is authored separately and
  reports the numbers in `results_tables/`.

## License

MIT (see `LICENSE`).
