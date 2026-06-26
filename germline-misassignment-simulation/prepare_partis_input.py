"""Prepare FASTA files for partis annotation of original and simulated MRCAs.

Writes two pooled FASTAs (all subjects combined):
  - original_mrcas.fasta: original MRCAs (child_heavy from naive-to-MRCA PCPs)
  - simulated_mrcas.fasta: simulated MRCAs (all replicates)

Each sequence gets a unique ID encoding sample_id, family, and replicate,
so annotations can be mapped back after partis runs.

Usage:
    python prepare_partis_input.py
    python prepare_partis_input.py --simulated-mrcas _output/simulated_mrcas.csv.gz --output-dir _output/partis_input
"""

import argparse
import os
import sys

# Add repo root to path so dnsmex can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pandas as pd

from dnsmex.dxsm_data import filter_valid_pcps, dataset_dict


RODRIGUEZ_WITH_NAIVE = dataset_dict["v1rodriguezWithNaive"]
RODRIGUEZ_NO_NAIVE = dataset_dict["v1rodriguez"]

# Default paths (relative to repo root)
DEFAULT_SIMULATED_MRCAS = "_output/simulated_mrcas.csv.gz"
DEFAULT_OUTPUT_DIR = "_output/partis_input"


def load_naive_to_mrca_pcps(with_naive_path, no_naive_path):
    """Load naive-to-MRCA PCPs, filtered to families present in the no-naive dataset.

    Same logic as simulate_mrca.py to ensure we use the same set of PCPs.
    """
    with_naive = pd.read_csv(with_naive_path)
    no_naive = pd.read_csv(no_naive_path)

    naive_pcps = with_naive[with_naive["parent_is_naive"] == True].copy()

    nn_families = set(
        zip(no_naive["sample_id"], no_naive["family"].astype(str))
    )
    naive_pcps["family_str"] = naive_pcps["family"].astype(str)
    mask = [
        (sid, fam) in nn_families
        for sid, fam in zip(naive_pcps["sample_id"], naive_pcps["family_str"])
    ]
    naive_pcps = naive_pcps[mask].copy()

    # Apply same filtering as simulate_mrca.py
    for col in ["parent_light", "child_light"]:
        if col not in naive_pcps.columns:
            naive_pcps[col] = ""
    naive_pcps = filter_valid_pcps(naive_pcps).reset_index(drop=True)
    naive_pcps["family_str"] = naive_pcps["family"].astype(str)

    return naive_pcps


def write_fasta(sequences, filepath):
    """Write a dict of {seq_id: sequence} to a FASTA file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


def make_seq_id(sample_id, family, suffix=""):
    """Create a unique sequence ID for partis."""
    base = f"{sample_id}_fam{family}"
    if suffix:
        base = f"{base}_{suffix}"
    return base


def main():
    parser = argparse.ArgumentParser(description="Prepare FASTA files for partis annotation")
    parser.add_argument("--simulated-mrcas", type=str, default=DEFAULT_SIMULATED_MRCAS,
                        help="Path to simulated MRCAs CSV")
    parser.add_argument("--output-dir", type=str, default=DEFAULT_OUTPUT_DIR,
                        help="Output directory for FASTA files")
    args = parser.parse_args()

    print("Loading naive-to-MRCA PCPs (same filtering as simulation)...")
    naive_pcps = load_naive_to_mrca_pcps(RODRIGUEZ_WITH_NAIVE, RODRIGUEZ_NO_NAIVE)
    print(f"  {len(naive_pcps)} valid PCPs")

    print(f"Loading simulated MRCAs from {args.simulated_mrcas}...")
    sim_df = pd.read_csv(args.simulated_mrcas)
    print(f"  {len(sim_df)} rows ({sim_df['replicate'].nunique()} replicates)")

    # Verify the PCPs match
    sim_pcps = set(zip(sim_df["sample_id"], sim_df["family"].astype(str)))
    orig_pcps = set(zip(naive_pcps["sample_id"], naive_pcps["family_str"]))
    assert sim_pcps == orig_pcps, (
        f"PCP mismatch: {len(sim_pcps)} simulated vs {len(orig_pcps)} original"
    )

    # Build pooled FASTAs (all subjects combined)
    original_seqs = {}
    simulated_seqs = {}

    for _, row in naive_pcps.iterrows():
        sid = row["sample_id"]
        fam = row["family_str"]
        original_seqs[make_seq_id(sid, fam, "orig")] = row["child_heavy"]

    for _, row in sim_df.iterrows():
        sid = row["sample_id"]
        fam = str(row["family"])
        simulated_seqs[make_seq_id(sid, fam, f"rep{row['replicate']}")] = row["simulated_mrca_heavy"]

    print(f"\nWriting pooled FASTAs to {args.output_dir}/")
    write_fasta(original_seqs, os.path.join(args.output_dir, "original_mrcas.fasta"))
    write_fasta(simulated_seqs, os.path.join(args.output_dir, "simulated_mrcas.fasta"))

    print(f"  original_mrcas.fasta: {len(original_seqs)} sequences")
    print(f"  simulated_mrcas.fasta: {len(simulated_seqs)} sequences")
    print(f"\nDone.")


if __name__ == "__main__":
    main()
