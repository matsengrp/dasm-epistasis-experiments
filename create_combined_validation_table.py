"""
Create a combined validation table with results from both validation approaches
plus DASM selection factors, unfiltered.

Approach 1: Productive (v1rodriguez) vs Non-productive (tangshm/out-of-frame)
Approach 2: Productive (v1jaffe) vs Thrifty expected

This script uses the same loading patterns as the original notebooks to ensure consistency.
"""

import os
import numpy as np
import pandas as pd

from utils import (
    add_column_aa_one_mutation_away_from_aa,
    add_germline_information,
    load_and_process_dasm_data,
)
from dnsmex.neutral_mutability import CachedNeutralMutabilityDataset

# Configuration matching the notebooks
NUMBERING_SCHEME = 'chothia'

# DASM config (same in both notebooks)
DASM_MODEL_NAME = "dasm_4m-v1jaffeCC+v1tangCC-joint"
DASM_DATASET_NAME = "v1rodriguez"

# Approach 1 config (from rates_analysis_productive_non_productive.ipynb)
APPROACH1_OBSERVED_DATASET = 'v1rodriguez'
APPROACH1_EXPECTED_DATASET = 'tangshm'
APPROACH1_BRANCH_LENGTH_METHOD = 'synonymous_mutation_freq_branch'

# Approach 2 config (from rates_analysis_productive_w_thrifty.ipynb)
APPROACH2_DATASET = 'v1jaffe'
APPROACH2_BRANCH_LENGTH_METHOD = 'mutation_frequency'
APPROACH2_BRANCH_LENGTH_SCALE_FACTOR = 1.60

PSEUDOCOUNT = 0.5


def load_dasm_data():
    """
    Load DASM data using the same pattern as both notebooks.
    Returns aggregated selection factors per (v_family, site, parent_aa, child_aa).
    """
    print("Loading DASM data using load_and_process_dasm_data()...")
    _, _, aa_site_subs_selection_df_germline = load_and_process_dasm_data(
        model_name=DASM_MODEL_NAME,
        dataset_name=DASM_DATASET_NAME,
        numbering_scheme=NUMBERING_SCHEME
    )

    # Aggregate matching the entrenchment analysis in v_families_dasm.ipynb:
    # 1. depth == 2 (only PCPs 2 levels from naive)
    # 2. one_mutation_away == True
    # 3. is_germline_codon == True
    # 4. drop NaN log_selection_factor
    # 5. remove rare AAs (< 10 unique PCPs per v_family+site+parent_aa)
    # 6. median of log_selection_factor (not log of median)
    before_grouping = aa_site_subs_selection_df_germline[
        (aa_site_subs_selection_df_germline.depth == 2) &
        (aa_site_subs_selection_df_germline.one_mutation_away == True) &
        (aa_site_subs_selection_df_germline.is_germline_codon == True)
    ].copy()
    before_grouping.dropna(subset=['log_selection_factor'], inplace=True)

    # Remove rare AAs: keep only parent_aa with >= 10 unique PCPs at that site+family
    counts_of_aa_site_family = before_grouping[['v_family', 'site', 'parent_aa', 'pcp_index']].drop_duplicates().groupby(['v_family', 'site', 'parent_aa']).size().reset_index(name='count')
    counts_of_aa_site_family = counts_of_aa_site_family[counts_of_aa_site_family['count'] >= 10]
    before_grouping = pd.merge(before_grouping, counts_of_aa_site_family[['v_family', 'site', 'parent_aa']], on=['v_family', 'site', 'parent_aa'], how='inner')

    dasm_summarized = before_grouping.groupby(
        ['v_family', 'site', 'parent_aa', 'selection_factor_target_aa']
    ).log_selection_factor.median().reset_index()

    dasm_summarized = dasm_summarized.rename(columns={
        'selection_factor_target_aa': 'child_aa',
        'log_selection_factor': 'dasm_log_selection_factor'
    })

    print(f"  Loaded {len(dasm_summarized)} DASM substitutions")
    return dasm_summarized


def load_approach1_data():
    """
    Load Approach 1 data (productive vs non-productive/out-of-frame).
    Uses the exact pattern from rates_analysis_productive_non_productive.ipynb cells 20-21.
    """
    print("\nLoading Approach 1 (productive vs non-productive) data...")

    # File paths matching the notebook
    observed_path = f'_ignore/observed_counts_from_productive/observed_mutation_rates_{APPROACH1_OBSERVED_DATASET}_bl_{APPROACH1_BRANCH_LENGTH_METHOD}_{NUMBERING_SCHEME}_per_aa_mutation_rates.csv'
    expected_path = f'_ignore/expected_counts_from_non_productive/expected_mutation_rates_{APPROACH1_EXPECTED_DATASET}_bl_{APPROACH1_BRANCH_LENGTH_METHOD}_{NUMBERING_SCHEME}_per_aa_mutation_rates.csv'

    observed_aa_df = pd.read_csv(observed_path, dtype={'site': str})
    expected_aa_df = pd.read_csv(expected_path, dtype={'site': str})

    print(f"  Observed: {len(observed_aa_df)} rows from {observed_path}")
    print(f"  Expected: {len(expected_aa_df)} rows from {expected_path}")

    # Add one_mutation_away column exactly as in the notebook (cell 20)
    add_column_aa_one_mutation_away_from_aa(observed_aa_df, 'parent_aa', 'child_aa')
    add_column_aa_one_mutation_away_from_aa(expected_aa_df, 'parent_aa', 'child_aa')

    # Merge exactly as in the notebook (cell 21)
    rates_aa_summarized = pd.merge(
        observed_aa_df[['site', 'parent_aa', 'child_aa', 'rate_mutcount', 'rate_mutcount_adjusted',
                        'mutcount_length', 'mutation_acquired', 'mutation_acquired_adjusted',
                        'v_family', 'one_mutation_away']],
        expected_aa_df[['site', 'parent_aa', 'child_aa', 'rate_mutcount', 'rate_mutcount_adjusted',
                        'mutcount_length', 'mutation_acquired', 'mutation_acquired_adjusted',
                        'v_family', 'one_mutation_away']],
        on=['site', 'parent_aa', 'child_aa', 'v_family', 'one_mutation_away'],
        suffixes=('_observed', '_expected')
    )

    # Calculate ratio exactly as in the notebook
    rates_aa_summarized['ratio_oof'] = (
        rates_aa_summarized['rate_mutcount_adjusted_observed'] /
        rates_aa_summarized['rate_mutcount_adjusted_expected']
    )
    rates_aa_summarized['log_ratio_oof'] = np.log(rates_aa_summarized['ratio_oof'])

    # Rename columns for clarity in combined table
    result = rates_aa_summarized.rename(columns={
        'mutcount_length_observed': 'branch_length_oof_observed',
        'mutcount_length_expected': 'branch_length_oof_expected',
        'mutation_acquired_observed': 'observed_counts_oof',
        'mutation_acquired_expected': 'expected_counts_oof',
        'mutation_acquired_adjusted_observed': 'observed_counts_adj_oof',
        'mutation_acquired_adjusted_expected': 'expected_counts_adj_oof',
        'rate_mutcount_observed': 'rate_oof_observed',
        'rate_mutcount_expected': 'rate_oof_expected',
        'rate_mutcount_adjusted_observed': 'rate_adj_oof_observed',
        'rate_mutcount_adjusted_expected': 'rate_adj_oof_expected',
    })

    print(f"  Merged: {len(result)} rows")
    return result


def load_approach2_data():
    """
    Load Approach 2 data (productive vs Thrifty expected).
    Uses the exact pattern from rates_analysis_productive_w_thrifty.ipynb cells 4-7.
    """
    print("\nLoading Approach 2 (productive vs Thrifty) data...")

    # Load observed from CSV (cell 4 pattern)
    observed_path = f'_ignore/observed_counts_from_productive/observed_mutation_rates_{APPROACH2_DATASET}_bl_total_mutation_freq_branch_{NUMBERING_SCHEME}_per_aa_mutation_rates.csv'
    observed_aa_df = pd.read_csv(observed_path, dtype={'site': str})

    # Rename columns exactly as in notebook cell 4
    observed_aa_df = observed_aa_df.rename(columns={
        'mutation_acquired': 'observed_counts',
        'mutation_acquired_adjusted': 'observed_counts_adjusted',
        'mutcount_length': 'mutcount_length_observed'
    }).drop(columns=['rate_mutcount', 'rate_mutcount_adjusted'])

    print(f"  Observed: {len(observed_aa_df)} rows from {observed_path}")

    # Load expected from Thrifty cache (cell 6 pattern)
    print("  Loading Thrifty neutral model cache...")
    neutral_probs = CachedNeutralMutabilityDataset(
        dataset_nickname=APPROACH2_DATASET,
        branch_length_mode=APPROACH2_BRANCH_LENGTH_METHOD,
        branch_length_scale_factor=APPROACH2_BRANCH_LENGTH_SCALE_FACTOR,
        numbering_scheme=NUMBERING_SCHEME,
    )

    # Add germline information exactly as in notebook cell 6
    neutral_probs.aa_neutral_df = add_germline_information(
        neutral_probs.pcp_df,
        neutral_probs.aa_neutral_df,
        numbering_scheme=NUMBERING_SCHEME
    )

    # LEAF FILTERING exactly as in notebook cell 6
    pcp_indices_non_leaf = neutral_probs.pcp_df[~neutral_probs.pcp_df['child_is_leaf']].index

    # Aggregate expected counts exactly as in notebook cell 6
    expected_aa_df = neutral_probs.aa_neutral_df[
        (neutral_probs.aa_neutral_df.is_germline_codon == True) &
        (neutral_probs.aa_neutral_df.pcp_index.isin(pcp_indices_non_leaf))
    ].groupby(['site', 'current_aa', 'v_family', 'transition_aa']).substitution_probability.sum().reset_index()

    expected_aa_df = expected_aa_df.rename(columns={
        'current_aa': 'parent_aa',
        'transition_aa': 'child_aa',
        'substitution_probability': 'expected_counts'
    })

    print(f"  Expected (Thrifty): {len(expected_aa_df)} rows")

    # Merge exactly as in notebook cell 7
    merge_counts = pd.merge(
        expected_aa_df,
        observed_aa_df,
        on=['v_family', 'site', 'parent_aa', 'child_aa']
    )

    # Calculate ratio exactly as in notebook cell 7
    merge_counts['expected_counts_adjusted'] = merge_counts['expected_counts'] + PSEUDOCOUNT
    merge_counts['ratio_thrifty'] = (
        merge_counts['observed_counts_adjusted'] /
        merge_counts['expected_counts_adjusted']
    )
    merge_counts['log_ratio_thrifty'] = np.log(merge_counts['ratio_thrifty'])

    # Rename columns for clarity in combined table
    result = merge_counts.rename(columns={
        'mutcount_length_observed': 'branch_length_thrifty',
        'observed_counts': 'observed_counts_thrifty',
        'observed_counts_adjusted': 'observed_counts_adj_thrifty',
        'expected_counts': 'expected_counts_thrifty',
        'expected_counts_adjusted': 'expected_counts_adj_thrifty',
    })

    print(f"  Merged: {len(result)} rows")
    return result


def create_combined_table():
    """Create combined validation table with all approaches."""

    # Load all data
    dasm_df = load_dasm_data()
    approach1_df = load_approach1_data()
    approach2_df = load_approach2_data()

    # Select columns for joining
    approach1_cols = [
        'v_family', 'site', 'parent_aa', 'child_aa', 'one_mutation_away',
        'branch_length_oof_observed', 'branch_length_oof_expected',
        'observed_counts_adj_oof', 'expected_counts_adj_oof',
        'rate_adj_oof_observed', 'rate_adj_oof_expected',
        'ratio_oof', 'log_ratio_oof'
    ]

    approach2_cols = [
        'v_family', 'site', 'parent_aa', 'child_aa',
        'branch_length_thrifty',
        'observed_counts_adj_thrifty', 'expected_counts_adj_thrifty',
        'ratio_thrifty', 'log_ratio_thrifty'
    ]

    # Merge DASM with approach 1
    print("\nMerging datasets...")
    combined = pd.merge(
        dasm_df,
        approach1_df[approach1_cols],
        on=['v_family', 'site', 'parent_aa', 'child_aa'],
        how='outer'
    )

    # Merge with approach 2
    combined = pd.merge(
        combined,
        approach2_df[approach2_cols],
        on=['v_family', 'site', 'parent_aa', 'child_aa'],
        how='outer'
    )

    # Rename columns with prefix style and final names
    column_rename = {
        # OOF columns
        'observed_counts_adj_oof': 'oof_observed_counts',
        'expected_counts_adj_oof': 'oof_expected_counts',
        'rate_adj_oof_observed': 'oof_rate_observed',
        'rate_adj_oof_expected': 'oof_rate_expected',
        'ratio_oof': 'oof_ratio',
        'log_ratio_oof': 'oof_log_ratio',
        'branch_length_oof_observed': 'oof_branch_length_observed',
        'branch_length_oof_expected': 'oof_branch_length_expected',
        # Thrifty columns
        'observed_counts_adj_thrifty': 'thrifty_observed_counts',
        'expected_counts_adj_thrifty': 'thrifty_expected_counts',
        'ratio_thrifty': 'thrifty_ratio',
        'log_ratio_thrifty': 'thrifty_log_ratio',
        'branch_length_thrifty': 'thrifty_branch_length',
    }
    combined = combined.rename(columns=column_rename)

    # Final column order (counts before branch lengths)
    final_columns = [
        'v_family', 'site', 'parent_aa', 'child_aa', 'one_mutation_away',
        'dasm_log_selection_factor',
        'oof_observed_counts', 'oof_expected_counts', 'oof_rate_observed', 'oof_rate_expected',
        'oof_ratio', 'oof_log_ratio', 'oof_branch_length_observed', 'oof_branch_length_expected',
        'thrifty_observed_counts', 'thrifty_expected_counts',
        'thrifty_ratio', 'thrifty_log_ratio', 'thrifty_branch_length',
    ]
    combined = combined[final_columns]

    print(f"\nCombined table: {len(combined)} rows, {len(combined.columns)} columns")
    print(f"  Rows with DASM data: {combined.dasm_log_selection_factor.notna().sum()}")
    print(f"  Rows with OOF data: {combined.oof_log_ratio.notna().sum()}")
    print(f"  Rows with Thrifty data: {combined.thrifty_log_ratio.notna().sum()}")
    print(f"  Rows with all three: {(combined.dasm_log_selection_factor.notna() & combined.oof_log_ratio.notna() & combined.thrifty_log_ratio.notna()).sum()}")

    return combined


def filter_for_entrenched_sites(combined_df):
    """Filter combined table to only include entrenched sites."""
    from utils import load_entrenched_sites

    _, entrenched_sites_aas, _, _, _, _ = load_entrenched_sites(NUMBERING_SCHEME)

    # Rename columns to match combined table
    entrenched_sites_aas_renamed = entrenched_sites_aas.rename(columns={
        'amino_acid': 'parent_aa',
        'target_amino_acid': 'child_aa'
    })

    # Filter
    filtered_df = pd.merge(
        combined_df,
        entrenched_sites_aas_renamed[['site', 'v_family', 'parent_aa', 'child_aa']],
        on=['site', 'v_family', 'parent_aa', 'child_aa'],
        how='inner'
    )

    print(f"\nEntrenched table: {len(filtered_df)} rows")
    print(f"  Rows with DASM data: {filtered_df.dasm_log_selection_factor.notna().sum()}")
    print(f"  Rows with OOF data: {filtered_df.oof_log_ratio.notna().sum()}")
    print(f"  Rows with Thrifty data: {filtered_df.thrifty_log_ratio.notna().sum()}")
    print(f"  Rows with all three: {(filtered_df.dasm_log_selection_factor.notna() & filtered_df.oof_log_ratio.notna() & filtered_df.thrifty_log_ratio.notna()).sum()}")

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
        'oof_observed_counts', 'oof_expected_counts', 'oof_branch_length_observed', 'oof_branch_length_expected', 'oof_log_ratio',
        'thrifty_observed_counts', 'thrifty_expected_counts', 'thrifty_branch_length', 'thrifty_log_ratio',
    ]].round({
        'dasm_log_selection_factor': 2,
        'oof_branch_length_observed': 2,
        'oof_branch_length_expected': 2,
        'oof_log_ratio': 2,
        'thrifty_expected_counts': 2,
        'thrifty_log_ratio': 2,
    }).to_csv(entrenched_for_paper_path, index=False)
    print(f"Saved entrenched (for paper) to: {entrenched_for_paper_path}")

    # Display sample
    print("\nColumn names:")
    print(combined_df.columns.tolist())
    print("\nSample of entrenched table (first 5 rows):")
    print(entrenched_df.head().to_string())
