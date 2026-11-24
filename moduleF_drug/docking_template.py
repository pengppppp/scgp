# -*- coding: utf-8 -*-
"""
Optional module: molecular docking of candidate small molecules against the top
MR-supported secreted protein targets (AutoDock Vina).

STATUS: code template provided, NOT executed (requires protein PDB structures and
ligand SDF files prepared locally; see README). Included so the full therapeutic
translation workflow is reproducible end-to-end.

Requirements: vina (pip install vina), openbabel (for ligand preparation), and PDB
structures for the target proteins (AlphaFold EBI API used here to fetch predicted
structures automatically - public, no key needed).
"""
import os
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # package root, derived from this file's location (portable)
TAB = os.path.join(ROOT, "results", "tables")
EXT = os.path.join(ROOT, "data", "external")


def fetch_alphafold_structure(uniprot_id, out_path):
    """Fetch AlphaFold predicted structure (public EBI API)."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(out_path, "wb") as f:
        f.write(r.read())
    return out_path


def dock(target_pdb, ligand_sdf, center, box, out_prefix, exhaustiveness=8, n_poses=9):
    """Run AutoDock Vina docking (target_pdb/ligand_sdf pre-prepared)."""
    from vina import Vina
    v = Vina(sf_name="vina")
    v.set_receptor(target_pdb)
    v.set_ligand_from_file(ligand_sdf)
    v.compute_vina_maps(center=center, box_size=box)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = v.energies(n_poses=n_poses)
    v.write_poses(f"{out_prefix}_poses.pdbqt", n_poses=n_poses, overwrite=True)
    return energies


def main():
    """
    Expected workflow (fill in once Module D/F identify the top targets):
      1. For each top target gene: UniProt -> accession (rest.uniprot.org).
      2. AlphaFold structure -> PDB (fetch_alphafold_structure).
      3. Ligand SDF from PubChem/ChEMBL (approved drugs from F1/F2 tables).
      4. define binding box from the ligand-binding site annotation / literature.
      5. dock(...) -> table of binding affinities; optional 100 ns MD (GROMACS)
         and MM/PBSA for the top complexes.
    All steps are real-data executable; nothing is pre-computed or fabricated.
    """
    raise SystemExit(
        "This is a template. Execute after Module D/F define the top targets; "
        "see README section 4 and STATUS.md.")


if __name__ == "__main__":
    main()
