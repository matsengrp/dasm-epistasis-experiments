#!/usr/bin/env python3
"""
Contact Analysis Script: distance-based alternative to the ΔRSA burial analysis

## Overview

For every heavy-chain residue, measures proximity to the three genetically uncoupled
partners used in the solvent accessibility analysis (`run_sasa_analysis.py`):

| Partner  | Atoms                                                  |
|----------|--------------------------------------------------------|
| antigen  | all antigen chains                                     |
| light    | the light chain                                        |
| cdr3     | heavy-chain Chothia sites 95-102 (same range removed   |
|          | in the `heavy_no_vdj_junction` ΔRSA scenario)          |

Structures, chain assignments, and filtering are identical to `run_sasa_analysis.py`
(its `StructureFilter` and `ProteinMetadataExtractor` are reused).

## Metrics (per heavy-chain residue and partner)

- `min_dist_<partner>`: minimum heavy-atom distance (Å) between the residue and the partner.
- `n_atoms_5A_<partner>`: number of partner heavy atoms within 5 Å of any residue heavy atom.
- `n_residues_5A_<partner>`: number of partner residues with a heavy atom within 5 Å.

Only standard amino acid residues are considered (waters, ligands, and non-protein antigens
are ignored), and hydrogens are excluded.

For the CDR3 partner, sites 95-102 themselves get NaN (as in the ΔRSA analysis, where they are
removed), and partner residues within `MIN_SEQ_SEPARATION - 1` positions in chain order are
excluded so that covalently bonded neighbours (e.g. site 94 next to 95) do not register as
trivial contacts.

## Output

Per-residue CSV keyed on `pdb_id, chain_id, residue_number, insertion_code`, the same keys as the
SASA output, so the two can be merged.

## Usage

    python run_contact_analysis.py --scheme opig-chothia --organism "homo sapiens" --output _output/contacts_human_chothia.csv
"""

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
from scipy.spatial import cKDTree
from tqdm import tqdm

from run_sasa_analysis import (
    HEAVY_REMOVAL_SCENARIOS,
    PDB_BASE_DIR,
    ProteinMetadataExtractor,
    StructureFilter,
)

CONTACT_CUTOFF = 5.0  # Å
CDR3_RANGE = HEAVY_REMOVAL_SCENARIOS["heavy_no_vdj_junction"]  # (95, 102)
MIN_SEQ_SEPARATION = 3  # CDR3 partner residues must be >= 3 positions away in chain order


def _standard_residues(chain):
    return [res for res in chain if res.id[0] == " " and is_aa(res, standard=True)]


def _heavy_atom_coords(residues):
    """Return (coords, residue_index) for all non-hydrogen atoms of the residues."""
    coords, owner = [], []
    for i, res in enumerate(residues):
        for atom in res:
            if atom.element in ("H", "D"):
                continue
            coords.append(atom.coord)
            owner.append(i)
    return np.array(coords, dtype=float).reshape(-1, 3), np.array(owner, dtype=int)


def _proximity(res_coords, partner_tree, partner_owner):
    """Minimum distance, atom count and residue count within the cutoff."""
    min_dist = partner_tree.query(res_coords, k=1)[0].min()
    hits = partner_tree.query_ball_point(res_coords, r=CONTACT_CUTOFF)
    atom_idx = np.unique(np.concatenate([np.array(h, dtype=int) for h in hits]))
    n_residues = len(np.unique(partner_owner[atom_idx])) if len(atom_idx) else 0
    return float(min_dist), len(atom_idx), n_residues


def compute_contacts(pdb_path, heavy_id, light_id, antigen_ids):
    """Per-residue proximity of the heavy chain to antigen, light chain and CDR3."""
    pdb_id = os.path.basename(pdb_path).removesuffix(".pdb")
    model = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(pdb_id, pdb_path)[0]

    missing = [c for c in [heavy_id, light_id, *antigen_ids] if c not in model]
    if missing:
        raise ValueError(f"{pdb_id}: chains {missing} not found in structure")

    heavy = _standard_residues(model[heavy_id])
    light = _standard_residues(model[light_id])
    antigen = [res for c in antigen_ids for res in _standard_residues(model[c])]
    if not heavy or not light or not antigen:
        raise ValueError(f"{pdb_id}: heavy, light or antigen chain has no amino acid residues")

    trees = {}
    for name, residues in [("antigen", antigen), ("light", light)]:
        coords, owner = _heavy_atom_coords(residues)
        trees[name] = (cKDTree(coords), owner)

    start, end = CDR3_RANGE
    cdr3_idx = [i for i, res in enumerate(heavy) if start <= res.id[1] <= end]
    if not cdr3_idx:
        raise ValueError(f"{pdb_id}: no heavy-chain residues in CDR3 range {CDR3_RANGE}")
    cdr3_set = set(cdr3_idx)

    rows = []
    for i, res in enumerate(heavy):
        res_coords, _ = _heavy_atom_coords([res])
        row = {
            "pdb_id": pdb_id,
            "chain_id": heavy_id,
            "residue_number": res.id[1],
            "insertion_code": res.id[2].strip(),
            "amino_acid": protein_letters_3to1[res.get_resname()],
        }
        for name in ["antigen", "light"]:
            tree, owner = trees[name]
            (row[f"min_dist_{name}"], row[f"n_atoms_5A_{name}"],
             row[f"n_residues_5A_{name}"]) = _proximity(res_coords, tree, owner)

        partner = [j for j in cdr3_idx if abs(j - i) >= MIN_SEQ_SEPARATION]
        if i in cdr3_set or not partner:
            row.update({"min_dist_cdr3": np.nan, "n_atoms_5A_cdr3": np.nan,
                        "n_residues_5A_cdr3": np.nan})
        else:
            coords, owner = _heavy_atom_coords([heavy[j] for j in partner])
            (row["min_dist_cdr3"], row["n_atoms_5A_cdr3"],
             row["n_residues_5A_cdr3"]) = _proximity(res_coords, cKDTree(coords), owner)
        rows.append(row)

    return pd.DataFrame(rows)


def _worker(job):
    try:
        return compute_contacts(*job), None
    except Exception as e:
        return None, f"{os.path.basename(job[0])}: {e}"


def main():
    parser = argparse.ArgumentParser(
        description="Per-residue heavy-atom distances and 5 Å contact counts between the heavy "
        "chain and antigen, light chain, and CDR3 (sites 95-102).",
    )
    parser.add_argument("--scheme", default="opig-chothia",
                        choices=["opig-imgt", "opig-chothia", "rcsb"])
    parser.add_argument("--pdb-dir", type=str, help="Custom PDB directory (overrides --scheme)")
    parser.add_argument("--output", "-o", type=str, required=True)
    parser.add_argument("--organism", type=str)
    parser.add_argument("--max-files", type=int)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    pdb_directory = args.pdb_dir or f"{PDB_BASE_DIR}/pdb/{args.scheme}"
    if not os.path.exists(pdb_directory):
        sys.exit(f"Error: PDB directory not found: {pdb_directory}")

    metadata = ProteinMetadataExtractor(verbose=args.verbose)
    pdb_files = StructureFilter(metadata, verbose=args.verbose).apply_filters(
        pdb_directory, args.organism, args.max_files
    )
    if not pdb_files:
        sys.exit("Error: no PDB files passed the filters")

    jobs = []
    for f in pdb_files:
        info = metadata.get_chain_info(f.stem)
        jobs.append((str(f), info["heavy_chain"], info["light_chain"], info["antigen_chains"]))

    results, failures = [], []
    with ProcessPoolExecutor(max_workers=args.n_jobs) as pool:
        for df, err in tqdm(pool.map(_worker, jobs, chunksize=8), total=len(jobs)):
            if err:
                failures.append(err)
            else:
                results.append(df)

    print(f"Processed {len(results)} structures, {len(failures)} failed")
    for err in failures[:20]:
        print(f"  {err}")
    if not results:
        sys.exit("Error: no structures processed")

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(args.output, index=False)
    print(f"Saved {len(combined)} residues to {args.output}")


if __name__ == "__main__":
    main()
