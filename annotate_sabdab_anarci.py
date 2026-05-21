#!/usr/bin/env python3
"""
Annotate SAbDab structures with V/J gene assignments using ANARCI.

Reads the SAbDab bulk download TSV and PDB files, extracts chain sequences,
runs ANARCI with germline assignment, and produces an annotation file
compatible with run_sasa_analysis.py.

Usage:
    python annotate_sabdab_anarci.py --output data/sabdab/sabdab_anarci_annotations.tsv
    python annotate_sabdab_anarci.py --max-pdbs 10 --output test_annotations.tsv
"""

import argparse
import os
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from anarci import anarci
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SABDAB_TSV = SCRIPT_DIR / "data/sabdab/sabdab_summary_all_2024-01-26.tsv"
DEFAULT_PDB_DIR = Path("/fh/fast/matsen_e/shared/bcr-mut-sel/sabdab/pdb-db/pdb/opig-imgt")

# Map SAbDab heavy_species to ANARCI allowed_species
SPECIES_MAP = {
    "homo sapiens": "human",
    "mus musculus": "mouse",
    "rattus norvegicus": "rat",
    "oryctolagus cuniculus": "rabbit",
    "macaca mulatta": "rhesus",
    "sus scrofa": "pig",
    "vicugna pacos": "alpaca",
    "lama glama": "alpaca",
    "bos taurus": "cow",
}


def extract_chain_sequence(pdb_path, chain_id):
    """Extract amino acid sequence for a specific chain from a PDB file.

    SAbDab uses lowercase chain IDs to distinguish duplicate copies of the
    same chain. PDB files only have uppercase IDs, so we fall back to
    uppercase if the exact chain ID isn't found.
    """
    try:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure("pdb", str(pdb_path))
        model = structure[0]
    except Exception:
        return None

    # Try exact chain ID first, then uppercase fallback
    if chain_id in model:
        chain = model[chain_id]
    elif chain_id.upper() in model:
        chain = model[chain_id.upper()]
    else:
        return None

    ppb = PPBuilder()
    peptides = ppb.build_peptides(chain)
    if not peptides:
        return None

    return "".join(str(pp.get_sequence()) for pp in peptides)


def run_anarci_on_sequence(sequence, allowed_species=None):
    """Run ANARCI on a single sequence and return V/J gene assignments.

    Returns (v_gene, j_gene) or (None, None) on failure.
    """
    kwargs = dict(scheme="imgt", output=False, assign_germline=True)
    if allowed_species:
        kwargs["allowed_species"] = allowed_species

    try:
        numbering, alignment_details, hit_tables = anarci(
            [("query", sequence)], **kwargs
        )
    except Exception:
        return None, None

    if not alignment_details or not alignment_details[0]:
        return None, None

    detail = alignment_details[0][0]
    germlines = detail.get("germlines", {})

    v_gene = germlines.get("v_gene", [None])[0]
    j_gene = germlines.get("j_gene", [None])[0]

    # v_gene/j_gene are tuples like ('human', 'IGHV3-23*03') or None
    if isinstance(v_gene, tuple):
        v_gene = v_gene[1]
    if isinstance(j_gene, tuple):
        j_gene = j_gene[1]

    return v_gene, j_gene


def annotate_pdb_chains(args):
    """Annotate all entries for a single PDB file.

    Groups entries by unique (Hchain, Lchain) pairs to avoid redundant ANARCI
    calls when the same antibody appears with different antigens. All rows are
    kept in the output.

    Returns a list of result dicts, or empty list on failure.
    """
    pdb_id, entries, pdb_dir = args

    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        return []

    # Cache: (chain_id) -> sequence
    seq_cache = {}
    # Cache: (h_chain, l_chain, species) -> (va, ja, vb, jb)
    gene_cache = {}

    results = []
    for row in entries:
        h_chain = row["Hchain"]
        l_chain = row["Lchain"]
        species = row["heavy_species"]

        # Extract sequences (cached per chain)
        if h_chain not in seq_cache:
            seq_cache[h_chain] = extract_chain_sequence(pdb_path, h_chain)
        if l_chain not in seq_cache:
            seq_cache[l_chain] = extract_chain_sequence(pdb_path, l_chain)

        h_seq = seq_cache[h_chain]
        l_seq = seq_cache[l_chain]
        if not h_seq or not l_seq:
            continue

        # Run ANARCI (cached per unique chain pair + species)
        cache_key = (h_chain, l_chain, species)
        if cache_key not in gene_cache:
            anarci_species = SPECIES_MAP.get(species)
            allowed = [anarci_species] if anarci_species else None
            vb, jb = run_anarci_on_sequence(h_seq, allowed_species=allowed)
            va, ja = run_anarci_on_sequence(l_seq, allowed_species=allowed)
            gene_cache[cache_key] = (va, ja, vb, jb)

        va, ja, vb, jb = gene_cache[cache_key]
        abid = f"{pdb_id}{h_chain}{l_chain}"

        results.append({
            "organism": species if pd.notna(species) else "",
            "pdbid": pdb_id,
            "abid": abid,
            "va": va if va else "unknown",
            "ja": ja if ja else "unknown",
            "vb": vb if vb else "unknown",
            "jb": jb if jb else "unknown",
            "chainseq_a": l_seq,
            "chainseq_b": h_seq,
        })

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Annotate SAbDab structures with ANARCI V/J gene assignments."
    )
    parser.add_argument(
        "--sabdab-tsv",
        type=str,
        default=str(DEFAULT_SABDAB_TSV),
        help="Path to SAbDab bulk download TSV",
    )
    parser.add_argument(
        "--pdb-dir",
        type=str,
        default=str(DEFAULT_PDB_DIR),
        help="Directory containing PDB files",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output TSV file path",
    )
    parser.add_argument(
        "--max-pdbs",
        type=int,
        default=None,
        help="Maximum number of entries to process (for testing)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, cpu_count()),
        help="Number of parallel workers (default: min(8, cpu_count))",
    )
    parser.add_argument(
        "--save-every",
        type=int,
        default=500,
        help="Save intermediate results every N entries",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from partial results file if it exists",
    )
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    if not pdb_dir.exists():
        print(f"Error: PDB directory not found: {pdb_dir}")
        sys.exit(1)

    # Load SAbDab summary
    sabdab = pd.read_table(args.sabdab_tsv)
    print(f"Loaded {len(sabdab)} entries from SAbDab summary")

    # Filter to entries with both H and L chains
    sabdab = sabdab[sabdab.Hchain.notna() & sabdab.Lchain.notna()].copy()
    print(f"Entries with H+L chains: {len(sabdab)}")

    # Group entries by PDB ID for efficient processing
    grouped = {}
    for _, row in sabdab.iterrows():
        pdb_id = row["pdb"]
        if pdb_id not in grouped:
            grouped[pdb_id] = []
        grouped[pdb_id].append(row.to_dict())

    pdb_ids = list(grouped.keys())
    if args.max_pdbs:
        pdb_ids = pdb_ids[:args.max_pdbs]
        print(f"Limited to {len(pdb_ids)} PDBs")

    # Resume from partial results if requested
    results = []
    done_pdbs = set()
    partial_path = Path(f"{args.output}.partial")
    if args.resume and partial_path.exists():
        partial_df = pd.read_csv(partial_path, sep="\t")
        results = partial_df.to_dict("records")
        done_pdbs = set(partial_df["pdbid"].unique())
        print(f"Resuming: loaded {len(results)} entries from {len(done_pdbs)} PDBs")
        pdb_ids = [p for p in pdb_ids if p not in done_pdbs]
        print(f"Remaining: {len(pdb_ids)} PDBs")

    work_items = [(pdb_id, grouped[pdb_id], pdb_dir) for pdb_id in pdb_ids]
    print(f"Processing {len(work_items)} PDBs ({sum(len(g) for _, g, _ in work_items)} entries) with {args.workers} workers...")

    if args.workers == 1:
        for item in tqdm(work_items, desc="Annotating"):
            batch = annotate_pdb_chains(item)
            results.extend(batch)
            if args.save_every and len(results) >= args.save_every and len(results) % args.save_every < len(batch):
                pd.DataFrame(results).to_csv(f"{args.output}.partial", sep="\t", index=False)
    else:
        with Pool(args.workers) as pool:
            for batch in tqdm(
                pool.imap(annotate_pdb_chains, work_items),
                total=len(work_items),
                desc="Annotating",
            ):
                results.extend(batch)
                if args.save_every and len(results) >= args.save_every and len(results) % args.save_every < len(batch):
                    pd.DataFrame(results).to_csv(f"{args.output}.partial", sep="\t", index=False)

    if not results:
        print("No entries were successfully annotated")
        sys.exit(1)

    df = pd.DataFrame(results)
    df.to_csv(args.output, sep="\t", index=False)

    # Clean up partial file
    partial_path = Path(f"{args.output}.partial")
    if partial_path.exists():
        partial_path.unlink()

    # Summary
    print(f"\n=== Annotation Summary ===")
    print(f"Total annotated: {len(df)}")
    print(f"  with V heavy gene: {(df.vb != 'unknown').sum()} ({100*(df.vb != 'unknown').mean():.1f}%)")
    print(f"  with J heavy gene: {(df.jb != 'unknown').sum()} ({100*(df.jb != 'unknown').mean():.1f}%)")
    print(f"  with V light gene: {(df.va != 'unknown').sum()} ({100*(df.va != 'unknown').mean():.1f}%)")
    print(f"  with J light gene: {(df.ja != 'unknown').sum()} ({100*(df.ja != 'unknown').mean():.1f}%)")
    print(f"Unique organisms: {df.organism.nunique()}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
