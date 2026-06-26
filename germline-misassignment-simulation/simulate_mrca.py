"""Simulate MRCA sequences from naive-to-MRCA PCPs using Thrifty+DASM.

For each naive-to-MRCA PCP, refits the branch length using the DASM crepe
(so it's calibrated to Thrifty+DASM, not just neutral Thrifty), then simulates
a new child sequence. Repeats for multiple replicates to estimate per-allele
partis misassignment rates.

Usage:
    python simulate_mrca.py --replicates 100 --output _output/simulated_mrcas.csv.gz
    python simulate_mrca.py --replicates 2 --max-pcps 50 --output _output/simulated_mrcas_test.csv.gz
"""

import argparse
import os
import sys

# Add repo root to path so dnsmex can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pandas as pd
import torch
from tqdm.auto import tqdm

from netam.common import chunked
from netam.framework import (
    load_crepe,
    codon_probs_of_parent_seq,
    sample_sequence_from_codon_probs,
    add_shm_model_outputs_to_pcp_df,
)
from netam import pretrained
from netam.dasm import DASMBurrito, DASMDataset
from netam.dnsm import DNSMBurrito, DNSMDataset
from dnsmex.dxsm_data import filter_valid_pcps, dataset_dict


# --- Paths ---
RODRIGUEZ_WITH_NAIVE = dataset_dict["v1rodriguezWithNaive"]
RODRIGUEZ_NO_NAIVE = dataset_dict["v1rodriguez"]
DASM_CREPE_PATH = "/home/nharel/re/dnsm-experiments-1/dasm-train/trained_models/dasm_4m-v1jaffeCC+v1tangCC-joint"

SEED = 42

# --- Branch length fitting (extracted from dnsmex/simulation.py) ---

_dxsm_classes_of_name = {
    "dasm": (DASMDataset, DASMBurrito),
    "dnsm": (DNSMDataset, DNSMBurrito),
}


def apply_func_to_chunked_df(df, func, batch_size, concat_func=torch.concat):
    if batch_size is None:
        chunk_size = len(df)
    else:
        n_chunks = int(len(df) / batch_size + 1)
        chunk_size = (len(df) // n_chunks) + 1
    position_chunks = list(chunked(range(len(df)), chunk_size))
    n_chunks = len(position_chunks)

    df_chunks = (df.iloc[position_chunk] for position_chunk in position_chunks)

    results = []
    for idx, df_chunk in enumerate(df_chunks):
        if n_chunks > 1:
            print(f"Processing batch {idx + 1} of {n_chunks}")
        results.append(func(df_chunk))
        del df_chunk
    if concat_func is not None:
        return concat_func(results)
    else:
        return results


def fit_branch_lengths(selection_crepe, pcp_df, batch_size=500_000):
    """Fit branch lengths to a pcp_df using the provided selection crepe.

    Extracted from dnsmex/simulation.py in dnsm-experiments-1.
    """
    model_type = selection_crepe.model.hyperparameters["model_type"]
    print("model type: ", model_type)
    dataset_cls, burrito_cls = _dxsm_classes_of_name[model_type]
    known_token_count = selection_crepe.model.hyperparameters["known_token_count"]
    neutral_crepe = pretrained.load(selection_crepe.model.neutral_model_name)
    multihit_model = pretrained.load_multihit(selection_crepe.model.multihit_model_name)
    pcp_df = add_shm_model_outputs_to_pcp_df(pcp_df, neutral_crepe)
    pcp_df["in_train"] = False

    def _branch_length_fit_helper(df_chunk):
        _, val_dataset = dataset_cls.train_val_datasets_of_pcp_df(
            df_chunk, known_token_count, multihit_model=multihit_model
        )
        burrito = burrito_cls(
            None,
            val_dataset,
            selection_crepe.model,
        )
        burrito.standardize_and_optimize_branch_lengths()
        result = burrito.val_dataset.branch_lengths
        del burrito, val_dataset
        return result

    return apply_func_to_chunked_df(pcp_df, _branch_length_fit_helper, batch_size)


# --- Data loading ---

def load_naive_to_mrca_pcps(with_naive_path, no_naive_path):
    """Load naive-to-MRCA PCPs, filtered to families present in the no-naive dataset."""
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

    return naive_pcps


# --- Simulation ---

def simulate_mrcas(naive_pcps, selection_crepe, replicates, seed):
    """Simulate new MRCA sequences from naive-to-MRCA PCPs.

    For each PCP, simulates `replicates` child sequences using Thrifty+DASM,
    each with a different random seed.
    """
    neutral_crepe = pretrained.load(selection_crepe.model.neutral_model_name)
    multihit_model = pretrained.load_multihit(
        selection_crepe.model.multihit_model_name
    )

    results = []

    for rep in range(replicates):
        rep_seed = seed + rep
        torch.manual_seed(rep_seed)

        for _, row in tqdm(
            naive_pcps.iterrows(),
            total=len(naive_pcps),
            desc=f"Replicate {rep + 1}/{replicates} (seed={rep_seed})",
        ):
            naive_seq_heavy = row["parent_heavy"]
            branch_length = row["branch_length"]

            nt_parent_pair = (naive_seq_heavy, "")
            codon_probs = codon_probs_of_parent_seq(
                selection_crepe,
                nt_parent_pair,
                branch_length,
                neutral_crepe,
                multihit_model=multihit_model,
            )
            simulated_heavy = sample_sequence_from_codon_probs(codon_probs[0])

            results.append(
                {
                    "replicate": rep,
                    "replicate_seed": rep_seed,
                    "sample_id": row["sample_id"],
                    "family": row["family_str"],
                    "v_gene_heavy": row["v_gene_heavy"],
                    "j_gene_heavy": row["j_gene_heavy"],
                    "branch_length": branch_length,
                    "naive_heavy": naive_seq_heavy,
                    "simulated_mrca_heavy": simulated_heavy,
                }
            )

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="Simulate MRCA sequences from naive-to-MRCA PCPs")
    parser.add_argument("--replicates", type=int, default=100, help="Number of simulation replicates per PCP")
    parser.add_argument("--seed", type=int, default=SEED, help="Base random seed")
    parser.add_argument("--output", type=str, default="_output/simulated_mrcas.csv.gz", help="Output file path")
    parser.add_argument("--max-pcps", type=int, default=None, help="Limit number of PCPs (for testing)")
    parser.add_argument("--crepe-path", type=str, default=DASM_CREPE_PATH, help="Path to DASM crepe model")
    args = parser.parse_args()

    print("Loading naive-to-MRCA PCPs...")
    naive_pcps = load_naive_to_mrca_pcps(RODRIGUEZ_WITH_NAIVE, RODRIGUEZ_NO_NAIVE)
    print(f"  {len(naive_pcps)} naive-to-MRCA PCPs")

    if args.max_pcps is not None:
        naive_pcps = naive_pcps.head(args.max_pcps)
        print(f"  Limited to {len(naive_pcps)} PCPs for testing")

    print(f"\nLoading DASM crepe from {args.crepe_path}...")
    selection_crepe = load_crepe(args.crepe_path)

    print("\nRefitting branch lengths with DASM crepe...")
    # Add empty light chain columns for heavy-only data
    for col in ["parent_light", "child_light"]:
        if col not in naive_pcps.columns:
            naive_pcps[col] = ""
    naive_pcps = filter_valid_pcps(naive_pcps).reset_index(drop=True)
    print(f"  {len(naive_pcps)} valid PCPs after filtering")
    new_bls = fit_branch_lengths(selection_crepe, naive_pcps)
    naive_pcps["branch_length"] = new_bls.detach().numpy()
    print(f"  Branch lengths refitted")

    print(f"\nSimulating {args.replicates} replicates per PCP (seed={args.seed})...")
    sim_df = simulate_mrcas(naive_pcps, selection_crepe, args.replicates, args.seed)

    print(f"\nSaving {len(sim_df)} simulated MRCAs to {args.output}...")
    sim_df.to_csv(args.output, index=False)
    print("Done.")


if __name__ == "__main__":
    main()
