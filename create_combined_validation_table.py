"""
Create a combined validation table with results from both validation approaches
plus DASM selection factors, unfiltered.

Approach 1: Productive (v1rodriguez) vs Non-productive (tangshm/out-of-frame)
Approach 2: Productive (v1jaffe+v1tang) vs Thrifty expected

This script reads pre-computed merged data from the notebook outputs in figures/
rather than recomputing from raw data.
"""

import os
import pandas as pd

from utils import load_entrenched_sites

# Configuration
NUMBERING_SCHEME = 'chothia'
PSEUDOCOUNT = 0.5
MIN_EXPECTED_COUNTS = 5

# Pre-computed notebook output files
THRIFTY_CSV = 'figures/rates_productive_w_thrifty_v1jaffe+v1tang_compare_dasm_rates_unfiltered.csv'
OOF_CSV = 'figures/rates_productive_non_productive_compare_dasm_rates_unfiltered.csv'


def load_approach1_data():
    """
    Load Approach 1 data (productive vs non-productive/out-of-frame)
    from pre-computed notebook output.
    """
    print("Loading Approach 1 (productive vs non-productive) from notebook CSV...")
    oof = pd.read_csv(OOF_CSV, dtype={'site': str})
    print(f"  Loaded {len(oof)} rows from {OOF_CSV}")

    result = oof.rename(columns={
        'mutcount_length_observed': 'oof_branch_length_observed',
        'mutcount_length_expected': 'oof_branch_length_expected',
        'mutation_acquired_observed': 'oof_observed_counts',
        'mutation_acquired_expected': 'oof_expected_counts',
        'mutation_acquired_adjusted_observed': 'oof_observed_counts_adjusted',
        'mutation_acquired_adjusted_expected': 'oof_expected_counts_adjusted',
        'rate_mutcount_adjusted_observed': 'oof_rate_observed',
        'rate_mutcount_adjusted_expected': 'oof_rate_expected',
        'ratio': 'oof_ratio',
        'log_ratio': 'oof_log_ratio',
        'log_selection_factor': 'dasm_log_selection_factor',
    })

    return result


def load_approach2_data():
    """
    Load Approach 2 data (productive vs Thrifty expected)
    from pre-computed notebook output.
    """
    print(f"Loading Approach 2 (productive vs Thrifty) from notebook CSV...")
    thrifty = pd.read_csv(THRIFTY_CSV, dtype={'site': str})
    print(f"  Loaded {len(thrifty)} rows from {THRIFTY_CSV}")

    # Add pseudocount-adjusted expected counts to match original script
    thrifty['expected_counts_adjusted'] = thrifty['expected_counts'] + PSEUDOCOUNT

    result = thrifty.rename(columns={
        'mutcount_length_observed': 'thrifty_branch_length',
        'observed_counts': 'thrifty_observed_counts',
        'observed_counts_adjusted': 'thrifty_observed_counts_adjusted',
        'expected_counts': 'thrifty_expected_counts',
        'expected_counts_adjusted': 'thrifty_expected_counts_adjusted',
        'ratio': 'thrifty_ratio',
        'log_ratio': 'thrifty_log_ratio',
        'log_selection_factor': 'dasm_log_selection_factor',
    })

    return result


def create_combined_table():
    """Create combined validation table with all approaches."""

    approach1_df = load_approach1_data()
    approach2_df = load_approach2_data()

    # Select columns for joining
    approach1_cols = [
        'v_family', 'site', 'parent_aa', 'child_aa', 'one_mutation_away',
        'dasm_log_selection_factor',
        'oof_branch_length_observed', 'oof_branch_length_expected',
        'oof_observed_counts', 'oof_expected_counts',
        'oof_observed_counts_adjusted', 'oof_expected_counts_adjusted',
        'oof_rate_observed', 'oof_rate_expected',
        'oof_ratio', 'oof_log_ratio',
    ]

    approach2_cols = [
        'v_family', 'site', 'parent_aa', 'child_aa',
        'dasm_log_selection_factor',
        'thrifty_branch_length',
        'thrifty_observed_counts', 'thrifty_expected_counts',
        'thrifty_observed_counts_adjusted', 'thrifty_expected_counts_adjusted',
        'thrifty_ratio', 'thrifty_log_ratio',
    ]

    # Merge: use DASM from approach 1 as primary, fill from approach 2
    print("\nMerging datasets...")
    combined = pd.merge(
        approach1_df[approach1_cols],
        approach2_df[approach2_cols],
        on=['v_family', 'site', 'parent_aa', 'child_aa'],
        how='outer',
        suffixes=('', '_thrifty_src'),
    )

    # Consolidate DASM: prefer approach 1 value, fill missing from approach 2
    if 'dasm_log_selection_factor_thrifty_src' in combined.columns:
        combined['dasm_log_selection_factor'] = combined['dasm_log_selection_factor'].fillna(
            combined['dasm_log_selection_factor_thrifty_src']
        )
        combined.drop(columns=['dasm_log_selection_factor_thrifty_src'], inplace=True)

    # Flag rows where expected counts pass the minimum threshold for reliable ratios
    combined['oof_passes_filter'] = combined['oof_expected_counts'] >= MIN_EXPECTED_COUNTS
    combined['thrifty_passes_filter'] = combined['thrifty_expected_counts'] >= MIN_EXPECTED_COUNTS

    # Final column order
    final_columns = [
        'v_family', 'site', 'parent_aa', 'child_aa', 'one_mutation_away',
        'dasm_log_selection_factor',
        'oof_observed_counts', 'oof_expected_counts',
        'oof_observed_counts_adjusted', 'oof_expected_counts_adjusted',
        'oof_rate_observed', 'oof_rate_expected',
        'oof_ratio', 'oof_log_ratio', 'oof_branch_length_observed', 'oof_branch_length_expected',
        'oof_passes_filter',
        'thrifty_observed_counts', 'thrifty_expected_counts',
        'thrifty_observed_counts_adjusted', 'thrifty_expected_counts_adjusted',
        'thrifty_ratio', 'thrifty_log_ratio', 'thrifty_branch_length',
        'thrifty_passes_filter',
    ]
    combined = combined[final_columns]

    print(f"\nCombined table: {len(combined)} rows, {len(combined.columns)} columns")
    print(f"  Rows with DASM data: {combined.dasm_log_selection_factor.notna().sum()}")
    print(f"  Rows with OOF data: {combined.oof_log_ratio.notna().sum()}")
    print(f"  Rows with Thrifty data: {combined.thrifty_log_ratio.notna().sum()}")
    print(f"  Rows with all three: {(combined.dasm_log_selection_factor.notna() & combined.oof_log_ratio.notna() & combined.thrifty_log_ratio.notna()).sum()}")

    return combined


def filter_for_entrenched_sites(combined_df):
    """Filter combined table to only include entrenched sites.

    Uses a left join from the entrenched sites list so that all entrenched
    substitutions appear in the output, even those without validation data
    in any approach (they will have NaN for validation columns).
    """
    _, entrenched_sites_aas, _, _, _, _ = load_entrenched_sites(NUMBERING_SCHEME)

    # Rename columns to match combined table
    entrenched_sites_aas_renamed = entrenched_sites_aas.rename(columns={
        'amino_acid': 'parent_aa',
        'target_amino_acid': 'child_aa'
    })

    # Left join: start from entrenched sites, bring in validation data where available
    filtered_df = pd.merge(
        entrenched_sites_aas_renamed[['site', 'v_family', 'parent_aa', 'child_aa']].drop_duplicates(),
        combined_df,
        on=['site', 'v_family', 'parent_aa', 'child_aa'],
        how='left'
    )

    print(f"\nEntrenched table: {len(filtered_df)} rows")
    print(f"  Rows with DASM data: {filtered_df.dasm_log_selection_factor.notna().sum()}")
    print(f"  Rows with OOF data: {filtered_df.oof_log_ratio.notna().sum()}")
    print(f"  Rows with Thrifty data: {filtered_df.thrifty_log_ratio.notna().sum()}")
    print(f"  Rows with all three: {(filtered_df.dasm_log_selection_factor.notna() & filtered_df.oof_log_ratio.notna() & filtered_df.thrifty_log_ratio.notna()).sum()}")
    print(f"  Rows with NO validation data: {(filtered_df.dasm_log_selection_factor.isna() & filtered_df.oof_log_ratio.isna() & filtered_df.thrifty_log_ratio.isna()).sum()}")

    return filtered_df


if __name__ == '__main__':
    # Create output directory
    output_dir = '_output'
    os.makedirs(output_dir, exist_ok=True)

    # Create and save unfiltered table
    combined_df = create_combined_table()
    unfiltered_path = f'{output_dir}/combined_validation_table_unfiltered.csv'
    combined_df.to_csv(unfiltered_path, index=False)
    print(f"\nSaved unfiltered to: {unfiltered_path}")

    # Create and save entrenched table
    entrenched_df = filter_for_entrenched_sites(combined_df)
    entrenched_path = f'{output_dir}/combined_validation_table_entrenched.csv'
    entrenched_df.to_csv(entrenched_path, index=False)
    print(f"Saved entrenched to: {entrenched_path}")

    # Create paper-ready version with selected columns, sorted and rounded
    entrenched_for_paper_path = f'{output_dir}/combined_validation_table_entrenched_for_paper.csv'
    entrenched_df.sort_values(['site', 'v_family'])[[
        'v_family', 'site', 'parent_aa', 'child_aa', 'dasm_log_selection_factor',
        'oof_observed_counts_adjusted', 'oof_expected_counts_adjusted',
        'oof_branch_length_observed', 'oof_branch_length_expected',
        'oof_log_ratio', 'oof_passes_filter',
        'thrifty_observed_counts_adjusted', 'thrifty_expected_counts_adjusted',
        'thrifty_branch_length', 'thrifty_log_ratio', 'thrifty_passes_filter',
    ]].round({
        'dasm_log_selection_factor': 2,
        'oof_branch_length_observed': 2,
        'oof_branch_length_expected': 2,
        'oof_log_ratio': 2,
        'thrifty_expected_counts_adjusted': 2,
        'thrifty_branch_length': 2,
        'thrifty_log_ratio': 2,
    }).to_csv(entrenched_for_paper_path, index=False)
    print(f"Saved entrenched (for paper) to: {entrenched_for_paper_path}")

    # Display sample
    print("\nColumn names:")
    print(combined_df.columns.tolist())
    print("\nSample of entrenched table (first 5 rows):")
    print(entrenched_df.head().to_string())
