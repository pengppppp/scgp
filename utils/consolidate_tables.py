# -*- coding: utf-8 -*-
"""Consolidate the 50 result CSVs into 7 layer-based Excel workbooks for
GigaDB deposition (one sheet per original table). CSVs stay in the code
package for machine-readable provenance."""
import os
import pandas as pd

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(CODE, "results_tables")
DST = os.path.join(CODE, "..", "03_tables_for_GigaDB")

# layer -> (xlsx name, [csv stems])
LAYOUT = [
    ("S1_single_cell.xlsx",
     ["A1_qc_summary_GSE242425", "A1_qc_summary_GSE209781",
      "A3_lineage_panels_GSE242425", "A3_lineage_panels_GSE209781",
      "A4_lineage_composition_GSE242425", "A4_lineage_composition_GSE209781",
      "A5_cluster_markers_GSE242425", "A5_cluster_markers_GSE209781",
      "A6_condition_DEGs_per_lineage_GSE242425",
      "A6_condition_DEGs_per_lineage_GSE209781",
      "A7_state_scores_GSE242425", "A7_state_scores_GSE209781"]),
    ("S2_secretome.xlsx",
     ["B0_state_gene_provenance", "B1_secreted_candidate_panel",
      "B2_panel_plasma_detectability", "B2_plasma_detectability_summary"]),
    ("S3_bulk_validation.xlsx",
     ["C1_state_scores_all_cohorts", "C2_disease_vs_control_state_scores"]),
    ("S4_genetic_MR.xlsx",
     ["D1_instruments", "D2_MR_results_kidney", "D3_replication_enrichment",
      "D4_endpoint_genetic_correlation", "D5_coloc_results",
      "D6_pQTL_MR_replication", "D8_pQTL_instruments",
      "D8_panel_pQTL_MR"]),
    ("S5_independent_replication.xlsx",
     ["D7_independent_CKDGen_replication",
      "D9_eGFRcys_eQTL_wald", "D9_eGFRcys_pQTL_MR",
      "D9_eGFRcrea_trans_eQTL_wald", "D9_eGFRcrea_trans_pQTL_MR",
      "D9_D10_independent_outcome_summary", "D9_D10_background_rate_analysis",
      "D10_BBJ_eGFR_eQTL_wald", "D10_BBJ_eGFR_pQTL_MR",
      "D11_kuppe_celltype_expression", "D12_NID2_CKDbinary_MR",
      "D13_UACR_pQTL_MR", "D13_UACR_eQTL_wald", "D13_UACR_summary",
      "E1_replication_MR_results", "E2_replication_consistency"]),
    ("S6_druggability.xlsx",
     ["F1_drug_mapping", "F2_approved_drug_summary",
      "F3_approved_drug_records_flagged", "F4_druggability_enrichment",
      "F5_full_universe_drugs"]),
    ("S7_cross_layer.xlsx",
     ["T0_cross_layer_evidence", "T0_evidence_tiers",
      "T1_known_protein_controls"]),
]


# Excel sheet names are capped at 31 chars; these need shortening.
SHORTNAME = {
    "A4_lineage_composition_GSE209781": "A4_lineage_comp_GSE209781",
    "A4_lineage_composition_GSE242425": "A4_lineage_comp_GSE242425",
    "A6_condition_DEGs_per_lineage_GSE242425": "A6_DEGs_GSE242425",
    "A6_condition_DEGs_per_lineage_GSE209781": "A6_DEGs_GSE209781",
    "C2_disease_vs_control_state_scores": "C2_disease_vs_control_scores",
    "D7_independent_CKDGen_replication": "D7_CKDGen_replication",
    "D9_D10_independent_outcome_summary": "D9_D10_outcome_summary",
    "F3_approved_drug_records_flagged": "F3_drug_records_flagged",
}


# Excel sheet names are capped at 31 chars; these need shortening.
SHORTNAME = {
    "A4_lineage_composition_GSE209781": "A4_lineage_comp_GSE209781",
    "A4_lineage_composition_GSE242425": "A4_lineage_comp_GSE242425",
    "A6_condition_DEGs_per_lineage_GSE242425": "A6_DEGs_GSE242425",
    "A6_condition_DEGs_per_lineage_GSE209781": "A6_DEGs_GSE209781",
    "C2_disease_vs_control_state_scores": "C2_disease_vs_control_scores",
    "D7_independent_CKDGen_replication": "D7_CKDGen_replication",
    "D9_D10_independent_outcome_summary": "D9_D10_outcome_summary",
    "F3_approved_drug_records_flagged": "F3_drug_records_flagged",
}


# Excel sheet names are capped at 31 chars; these need shortening.
SHORTNAME = {
    "A4_lineage_composition_GSE209781": "A4_lineage_comp_GSE209781",
    "A4_lineage_composition_GSE242425": "A4_lineage_comp_GSE242425",
    "A6_condition_DEGs_per_lineage_GSE242425": "A6_DEGs_GSE242425",
    "A6_condition_DEGs_per_lineage_GSE209781": "A6_DEGs_GSE209781",
    "C2_disease_vs_control_state_scores": "C2_disease_vs_control_scores",
    "D7_independent_CKDGen_replication": "D7_CKDGen_replication",
    "D9_D10_independent_outcome_summary": "D9_D10_outcome_summary",
    "F3_approved_drug_records_flagged": "F3_drug_records_flagged",
}


def main():
    os.makedirs(DST, exist_ok=True)
    total = 0
    for xlsx, stems in LAYOUT:
        sheets = {}
        for s in stems:
            for ext in (".csv", ".csv.gz"):
                p = os.path.join(SRC, s + ext)
                if os.path.exists(p):
                    sheets[s] = pd.read_csv(p)
                    break
            else:
                print(f"  WARN missing {s}")
        out = os.path.join(DST, xlsx)
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            for name, df in sheets.items():
                sn = SHORTNAME.get(name, name)
                if len(sn) > 31:
                    sn = sn[:31]
                df.to_excel(w, sheet_name=sn, index=False)
        total += len(sheets)
        print(f"  {xlsx}: {len(sheets)} sheets")
    print(f"total {total} sheets across 7 workbooks")


if __name__ == "__main__":
    main()
