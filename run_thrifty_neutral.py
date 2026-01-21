#!/usr/bin/env python
"""
Simplified script to generate neutral mutation probabilities using CachedNeutralMutabilityDataset.

This script replaces the old NeutralMutationProbability implementation with the new
CachedNeutralMutabilityDataset from dnsm-experiments-1/dnsmex/neutral_mutability.py

Usage:
    python run_thrifty_neutral.py <dataset> [options]

Examples:
    # Generate with default settings (IMGT, mutation_frequency)
    python run_thrifty_neutral.py v1rodriguez

    # Chothia numbering with synonymous frequency
    python run_thrifty_neutral.py v1rodriguez --numbering-scheme chothia --branch-length-mode synonymous_frequency

    # Parallel execution for multiple datasets
    parallel python run_thrifty_neutral.py {} --numbering-scheme chothia ::: v1rodriguez v1tang v1jaffeCC
"""

import argparse
from dnsmex.neutral_mutability import CachedNeutralMutabilityDataset


def main():
    parser = argparse.ArgumentParser(
        description='Generate neutral mutation probabilities for a dataset using CachedNeutralMutabilityDataset',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s v1rodriguez
  %(prog)s v1rodriguez --numbering-scheme chothia
  %(prog)s v1rodriguez --branch-length-mode synonymous_frequency
  %(prog)s v1tang --branch-length-scale-factor 1.5

Parallel execution:
  parallel python %(prog)s {} --numbering-scheme chothia ::: v1rodriguez v1tang v1jaffeCC
        """
    )

    parser.add_argument('dataset_nickname', type=str,
                        help='Dataset name (e.g., v1rodriguez, v1tang, v1jaffeCC)')

    parser.add_argument('--branch-length-mode', type=str,
                        default='mutation_frequency',
                        choices=['constant', 'mutation_frequency', 'synonymous_frequency'],
                        help='Branch length calculation mode (default: mutation_frequency)')

    parser.add_argument('--branch-length', type=float, default=0.1,
                        help='Branch length value for constant mode (default: 0.1)')

    parser.add_argument('--branch-length-scale-factor', type=float, default=1.0,
                        help='Multiplicative scale factor for branch lengths (default: 1.0)')

    parser.add_argument('--numbering-scheme', type=str, default='imgt',
                        help='Antibody numbering scheme (default: imgt, options: imgt, chothia)')

    parser.add_argument('--subset-size', type=int, default=None,
                        help='Use only a subset of sequences for testing (default: None = use all)')

    args = parser.parse_args()

    print(f"\n{'='*70}")
    print(f"Generating neutral mutation probabilities")
    print(f"{'='*70}")
    print(f"Dataset: {args.dataset_nickname}")
    print(f"Numbering scheme: {args.numbering_scheme}")
    print(f"Branch length mode: {args.branch_length_mode}")
    if args.branch_length_mode == 'constant':
        print(f"Branch length: {args.branch_length}")
    print(f"Branch length scale factor: {args.branch_length_scale_factor}")
    if args.subset_size:
        print(f"Subset size: {args.subset_size}")
    print(f"{'='*70}\n")

    # Create the cached dataset - this will either load from cache or create new data
    neutral_dataset = CachedNeutralMutabilityDataset(
        dataset_nickname=args.dataset_nickname,
        branch_length_mode=args.branch_length_mode,
        branch_length=args.branch_length,
        branch_length_scale_factor=args.branch_length_scale_factor,
        numbering_scheme=args.numbering_scheme,
        subset_size=args.subset_size,
    )

    print(f"\n{'='*70}")
    print(f"✓ Complete! Data available in cache.")
    print(f"{'='*70}")
    print(f"Available DataFrames:")
    print(f"  - nucleotide_df: {len(neutral_dataset.nucleotide_df):,} rows")
    print(f"  - aa_neutral_df: {len(neutral_dataset.aa_neutral_df):,} rows")
    print(f"  - aa_to_any_neutral_df: {len(neutral_dataset.aa_to_any_neutral_df):,} rows")
    print(f"  - codon_neutral_df: {len(neutral_dataset.codon_neutral_df):,} rows")
    print(f"  - codon_to_any_neutral_df: {len(neutral_dataset.codon_to_any_neutral_df):,} rows")
    print(f"  - pcp_df: {len(neutral_dataset.pcp_df):,} rows")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
