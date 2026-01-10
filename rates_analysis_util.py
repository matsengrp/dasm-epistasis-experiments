"""
Utility functions for analyzing mutation rates across different V gene families.

This module contains shared functions used by:
- rates_analysis_productive_non_productive.ipynb
- rates_analysis_productive_w_thrifty.ipynb
"""

import pandas as pd
import os
from netam.sequences import AA_STR_SORTED, CODONS, translate_codon


def compare_mutation_rates_on_different_vfamilies(site_sub_probs_df_germline, site, vfamilies=['IGHV1', 'IGHV3', 'IGHV4'], branch_length_method='synonymous_mutation_freq_branch'):
    '''
    Calculate mutation rates at a specific site across different V gene families.
    Mutation rates are calculated for codon

    This function analyzes mutation rates at a given site across different V families,
    using branch length as a measure of evolutionary time.

    Parameters:
    -----------
    site_sub_probs_df_germline : pd.DataFrame
        DataFrame containing site substitution probabilities with columns:
        - 'site': site position
        - 'sample_id': sample identifier
        - 'family': family identifier
        - 'pcp_index': phylogenetic branch identifier
        - 'branch_length': evolutionary distance (substitutions per site)
        - 'germline_amino_acid': original amino acid at germline
        - 'is_germline_aa': boolean indicating if current AA matches germline
        - 'mutation': boolean indicating if site is mutated from germline
        - 'mutations_per_branch': total mutations accumulated on this branch

    vfamilies : list
        List of V gene families to compare (default: ['IGHV1', 'IGHV3', 'IGHV4'])
    site : int
        Site position to analyze

    Returns:
    --------
    list
        List of dictionaries containing mutation rate data for each V family:
        - 'vfamily': V gene family name
        - 'site': site position
        - 'branch_length': total branch length
        - 'mutations_per_branch': total mutations per branch
        - 'mutation_acquired': number of mutations acquired
        - 'rate': mutations per unit branch length
        - 'rate_mutcount': mutations per mutation count

    Notes:
    ------
    - Only analyzes branches where the site has germline amino acid identity
    - Handles cases with zero branch length or mutation count
    - Calculates both branch-length normalized and mutation-count normalized rates
    '''
    if branch_length_method not in ['synonymous_mutation_freq_branch', 'nonsynonymous_mutation_freq_branch', 'total_mutation_freq_branch']:
        raise ValueError("branch_length_method must be one of 'synonymous_mutation_freq_branch', 'nonsynonymous_mutation_freq_branch', or 'total_mutation_freq_branch'")

    # filter only rows with the two relevant sites
    cur_df = site_sub_probs_df_germline[(site_sub_probs_df_germline['site'] == site) & (site_sub_probs_df_germline.v_family.isin(vfamilies))].copy()

    # filter only rows where germline identity is known
    cur_df = cur_df[(cur_df['germline_amino_acid'].notna())]

    vfamily_results = []
    vfamily_results_per_aa = []
    vfamily_results_per_codon = []

    for vfamily in vfamilies:
        # filter for the current V family and only branches that have the germline codon at site
        vfamily_df = cur_df[(cur_df['v_family'] == vfamily) & (cur_df.is_germline_codon == True)].copy()
        if vfamily_df.empty:
            print(f"No data for V family {vfamily} at site {site}")
            continue
        # Calculate mutation rates for the current V family
        length_mutcount = vfamily_df[branch_length_method].sum()

        mutation_acquired = len(vfamily_df[(vfamily_df['nucleotide_mutation_count'] > 0) & (vfamily_df['child_aa'] != vfamily_df['parent_aa'])]) # any nonsynonymous mutations

        if length_mutcount == 0:
            print(f"Zero branch length or mutation count for V family {vfamily} at site {site}")
            continue

        # Calculate mutation rates
        rate_mutcount = mutation_acquired / length_mutcount
        vfamily_results.append({
            'v_family': vfamily,
            'site': site,
            'mutcount_length': length_mutcount,
            'mutation_acquired': mutation_acquired,
            'rate_mutcount': rate_mutcount
        })


        # calculate aa specific rates
        for amino_acid in vfamily_df['germline_amino_acid'].unique():
            aa_df = vfamily_df[vfamily_df['parent_aa'] == amino_acid]
            aa_length_mutcount = aa_df[branch_length_method].sum()

            if aa_df.empty:
                continue

            for target_amino_acid in AA_STR_SORTED:
                if amino_acid != target_amino_acid:
                    aa_mutation_acquired = len(aa_df[(aa_df['nucleotide_mutation_count'] > 0) & (aa_df['child_aa'] == target_amino_acid)])

                    # Calculate mutation rates for the specific amino acid
                    #rate_aa = aa_mutation_acquired / aa_length
                    rate_aa_mutcount = aa_mutation_acquired / aa_length_mutcount

                    vfamily_results_per_aa.append({
                        'v_family': f"{vfamily}",
                        'site': site,
                        'parent_aa': amino_acid,
                        'child_aa': target_amino_acid,
                        #'branch_length': aa_length,
                        'mutcount_length': aa_length_mutcount,
                        'mutation_acquired': aa_mutation_acquired,
                        #'rate': rate_aa,
                        'rate_mutcount': rate_aa_mutcount
                    })

         # calculate codon specific rates
        for codon in vfamily_df['germline_codon'].unique():
            codon_df = vfamily_df[vfamily_df['parent_codon'] == codon]
            codon_length_mutcount = codon_df[branch_length_method].sum()

            if codon_df.empty:
                continue

            for target_codon in CODONS:
                if codon != target_codon:
                    if sum(1 for a, b in zip(codon, target_codon) if a != b) > 1:
                        continue  # skip codons that are more than one mutation away
                    codon_mutation_acquired = len(codon_df[(codon_df['nucleotide_mutation_count'] > 0) & (codon_df['child_codon'] == target_codon)])

                    # Calculate mutation rates for the specific amino acid
                    rate_codon_mutcount = codon_mutation_acquired / codon_length_mutcount

                    vfamily_results_per_codon.append({
                        'v_family': f"{vfamily}",
                        'site': site,
                        'parent_codon': codon,
                        'child_codon': target_codon,
                        'parent_aa': translate_codon(codon),
                        'child_aa': translate_codon(target_codon),
                        'mutcount_length': codon_length_mutcount,
                        'mutation_acquired': codon_mutation_acquired,
                        'rate_mutcount': rate_codon_mutcount
                    })


    return vfamily_results, vfamily_results_per_aa, vfamily_results_per_codon


def compare_mutation_rates_on_different_backgrounds_for_all_sites(site_sub_probs_df_germline, output_path, branch_length_method='synonymous_mutation_freq_branch', remove_leaves=True, pcp_df=None):
    """
    Calculate mutation rates for all sites across different V gene families.

    Args:
        site_sub_probs_df_germline (pd.DataFrame): DataFrame with site substitution probabilities and germline info
        output_path (str): Base output path to save the results CSV files (without extension)
        branch_length_method (str): Method to use for branch length calculation
        remove_leaves (bool): Whether to remove leaf nodes from the analysis (default: True)
        pcp_df (pd.DataFrame): Parent-child pair DataFrame with 'child_is_leaf' column (required if remove_leaves=True)

    Returns:
        tuple: (results_df, results_per_aa_df) - DataFrames with results for overall and per-AA analysis
    """
    # Remove leaf nodes if requested
    if remove_leaves:
        if pcp_df is None:
            print("Warning: remove_leaves=True but pcp_df is None. Skipping leaf removal.")
        else:
            print("Removing leaf nodes from the analysis")
            pcp_df_for_merge = pcp_df[['child_is_leaf']].copy()
            pcp_df_for_merge['pcp_index'] = pcp_df_for_merge.index
            site_sub_probs_df_germline = pd.merge(
                site_sub_probs_df_germline,
                pcp_df_for_merge,
                on='pcp_index',
                how='inner'
            )
            site_sub_probs_df_germline = site_sub_probs_df_germline[
                ~site_sub_probs_df_germline['child_is_leaf']
            ]

    # Run for all sites
    all_results = []
    all_results_per_aa = []
    all_results_per_codon = []
    site_list = site_sub_probs_df_germline.site.unique().tolist()

    print(f"Processing {len(site_list)} sites...")

    for i, site in enumerate(site_list):
        if i % 10 == 0:  # Progress indicator
            print(f"Processing site {i+1}/{len(site_list)}: site {site}")

        vfamily_results, vfamily_results_per_aa, vfamily_results_per_codon = compare_mutation_rates_on_different_vfamilies(
            site_sub_probs_df_germline, site, branch_length_method=branch_length_method
        )

        # Extend the lists with results from this site
        all_results.extend(vfamily_results)
        all_results_per_aa.extend(vfamily_results_per_aa)
        all_results_per_codon.extend(vfamily_results_per_codon)

    # Convert to DataFrames
    print("Converting results to DataFrames...")
    results_df = pd.DataFrame(all_results)
    results_per_aa_df = pd.DataFrame(all_results_per_aa)
    results_per_codon_df = pd.DataFrame(all_results_per_codon)

    # Generate output file paths
    base_path = output_path.replace('.csv', '')  # Remove .csv if present
    overall_output_path = f"{base_path}_overall_mutation_rates.csv"
    per_aa_output_path = f"{base_path}_per_aa_mutation_rates.csv"
    per_codon_output_path = f"{base_path}_per_codon_mutation_rates.csv"

    # Save overall results
    if not results_df.empty:
        print(f'Saving overall results to {overall_output_path}')
        results_df.to_csv(overall_output_path, index=False)
        print(f"Overall results: {len(results_df)} rows saved")
    else:
        print("Warning: No overall results to save")

    # Save per-amino acid results
    if not results_per_aa_df.empty:
        print(f'Saving per-amino acid results to {per_aa_output_path}')
        results_per_aa_df.to_csv(per_aa_output_path, index=False)
        print(f"Per-AA results: {len(results_per_aa_df)} rows saved")
    else:
        print("Warning: No per-amino acid results to save")

    # Save per-codon results
    if not results_per_codon_df.empty:
        print(f'Saving per-codon results to {per_codon_output_path}')
        results_per_codon_df.to_csv(per_codon_output_path, index=False)
        print(f"Per-codon results: {len(results_per_codon_df)} rows saved")
    else:
        print("Warning: No per-codon results to save")

    return results_df, results_per_aa_df, results_per_codon_df
