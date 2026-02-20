"""
Utility functions for analyzing mutation rates across different V gene families.

This module contains shared functions used by:
- rates_analysis_productive_non_productive.ipynb
- rates_analysis_productive_w_thrifty.ipynb
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.odr import ODR, Model, RealData
from netam.sequences import AA_STR_SORTED, CODONS, translate_codon

from utils import sort_antibody_sites


def orthogonal_regression(x, y):
    """
    Orthogonal/Deming regression - appropriate when both variables have measurement error.

    When comparing two models that both have measurement error/uncertainty, OLS produces
    attenuated slopes - biased toward 0. This is called "regression dilution" or
    "attenuation bias." If the data truly lies along y=x but both x and y have independent
    noise, OLS will give a slope < 1. The higher the noise, the more attenuation.

    Orthogonal regression (also called Deming regression or total least squares) minimizes
    perpendicular distances to the line, treating both variables as having error.

    Parameters:
    -----------
    x : array-like
        Independent variable values
    y : array-like
        Dependent variable values

    Returns:
    --------
    tuple : (slope, intercept)
        The slope and intercept of the orthogonal regression line
    """
    def linear(B, x):
        return B[0] * x + B[1]

    linear_model = Model(linear)
    data = RealData(x, y)
    odr = ODR(data, linear_model, beta0=[1., 0.])
    output = odr.run()
    return output.beta[0], output.beta[1]  # slope, intercept


def compare_mutation_rates_on_different_vfamilies(site_sub_probs_df_germline, site, vfamilies=['IGHV1', 'IGHV3', 'IGHV4'], branch_length_method='synonymous_mutation_freq_branch', pseudocount=0.5):
    '''
    Calculate mutation rates at a specific site across different V gene families.

    This function analyzes mutation rates at a given site across different V families,
    using branch length as a measure of evolutionary time. Results are calculated at
    three levels of granularity:
    1. Per V-family (overall mutation rate)
    2. Per amino acid substitution (parent_aa -> child_aa)
    3. Per codon substitution (parent_codon -> child_codon, single mutations only)

    Parameters:
    -----------
    site_sub_probs_df_germline : pd.DataFrame
        DataFrame containing site substitution probabilities with columns:
        - 'site': site position
        - 'v_family': V gene family
        - 'germline_amino_acid': original amino acid at germline
        - 'germline_codon': original codon at germline
        - 'is_germline_codon': boolean indicating if parent codon matches germline
        - 'parent_aa', 'child_aa': parent and child amino acids
        - 'parent_codon', 'child_codon': parent and child codons
        - 'nucleotide_mutation_count': number of nucleotide mutations
        - branch_length_method column: branch length measure

    site : str or int
        Site position to analyze

    vfamilies : list
        List of V gene families to compare (default: ['IGHV1', 'IGHV3', 'IGHV4'])

    branch_length_method : str
        Column name to use for branch length calculation

    pseudocount : float
        Pseudocount to add to mutation counts before calculating adjusted rates.
        This is applied at the count level: adjusted_rate = (count + pseudocount) / branch_length.
        Default is 0.5 (Laplace smoothing). Output includes both raw and adjusted values.

    Returns:
    --------
    tuple of lists
        (vfamily_results, vfamily_results_per_aa, vfamily_results_per_codon)
        Each list contains dictionaries with mutation rate data including:
        - mutation_acquired: raw count
        - mutation_acquired_adjusted: count + pseudocount
        - rate_mutcount: raw rate (count / branch_length)
        - rate_mutcount_adjusted: adjusted rate ((count + pseudocount) / branch_length)
    '''
    if branch_length_method not in ['synonymous_mutation_freq_branch', 'nonsynonymous_mutation_freq_branch', 'total_mutation_freq_branch']:
        raise ValueError("branch_length_method must be one of 'synonymous_mutation_freq_branch', 'nonsynonymous_mutation_freq_branch', or 'total_mutation_freq_branch'")

    # filter only rows with the relevant site and V families
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
            'mutation_acquired_adjusted': mutation_acquired + pseudocount,
            'rate_mutcount': rate_mutcount,
            'rate_mutcount_adjusted': (mutation_acquired + pseudocount) / length_mutcount
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
                    rate_aa_mutcount = aa_mutation_acquired / aa_length_mutcount

                    vfamily_results_per_aa.append({
                        'v_family': f"{vfamily}",
                        'site': site,
                        'parent_aa': amino_acid,
                        'child_aa': target_amino_acid,
                        'mutcount_length': aa_length_mutcount,
                        'mutation_acquired': aa_mutation_acquired,
                        'mutation_acquired_adjusted': aa_mutation_acquired + pseudocount,
                        'rate_mutcount': rate_aa_mutcount,
                        'rate_mutcount_adjusted': (aa_mutation_acquired + pseudocount) / aa_length_mutcount
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

                    # Calculate mutation rates for the specific codon
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
                        'mutation_acquired_adjusted': codon_mutation_acquired + pseudocount,
                        'rate_mutcount': rate_codon_mutcount,
                        'rate_mutcount_adjusted': (codon_mutation_acquired + pseudocount) / codon_length_mutcount
                    })


    return vfamily_results, vfamily_results_per_aa, vfamily_results_per_codon


def compare_mutation_rates_on_different_backgrounds_for_all_sites(site_sub_probs_df_germline, output_path, branch_length_method='synonymous_mutation_freq_branch', remove_leaves=True, pcp_df=None, pseudocount=0.5):
    """
    Calculate mutation rates for all sites across different V gene families.

    Args:
        site_sub_probs_df_germline (pd.DataFrame): DataFrame with site substitution probabilities and germline info
        output_path (str): Base output path to save the results CSV files (without extension)
        branch_length_method (str): Method to use for branch length calculation
        remove_leaves (bool): Whether to remove leaf nodes from the analysis (default: True)
        pcp_df (pd.DataFrame): Parent-child pair DataFrame with 'child_is_leaf' column (required if remove_leaves=True)
        pseudocount (float): Pseudocount to add to mutation counts (default: 0.5)

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
            site_sub_probs_df_germline, site, branch_length_method=branch_length_method, pseudocount=pseudocount
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


def add_mutation_counts_per_branch_for_branch_length(df):
    """
    Add columns for synonymous and nonsynonymous nucleotide mutation frequencies
    to use as alternative branch length measures.

    This function calculates per-branch mutation counts and frequencies that can be used
    as alternative branch length normalization methods when calculating mutation rates.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns:
        - 'pcp_index': identifier for parent-child pairs
        - 'parent_codon': parent codon sequence
        - 'child_codon': child codon sequence
        - 'parent_aa': parent amino acid
        - 'child_aa': child amino acid

    Returns:
    --------
    pd.DataFrame
        Input DataFrame with additional columns:
        - 'nucleotide_mutation_count': number of nucleotide differences per codon
        - 'synonymous_mutation_freq_branch': synonymous mutations / sequence length
        - 'nonsynonymous_mutation_freq_branch': nonsynonymous mutations / sequence length
        - 'total_mutation_freq_branch': total mutations / sequence length

    Notes:
    ------
    The mutation frequency columns provide alternative normalization methods:
    - 'synonymous_mutation_freq_branch': Useful when comparing productive and non-productive
      data, as synonymous mutations are under the same selection pressure in both.
    - 'total_mutation_freq_branch': General-purpose normalization.
    """
    # Calculate nucleotide sequence length (number of codons * 3)
    df['seq_nuc_length'] = df.groupby('pcp_index').transform('size') * 3

    # Count nucleotide differences between parent and child codons
    df['nucleotide_mutation_count'] = df.apply(
        lambda row: sum(c1 != c2 for c1, c2 in zip(row['child_codon'], row['parent_codon'])),
        axis=1
    )

    # Mark amino acid changes
    df['aa_mutation'] = df['parent_aa'] != df['child_aa']

    # Separate synonymous and nonsynonymous nucleotide mutations
    df['synonymous_nucleotide_mutation_count'] = np.where(
        df['aa_mutation'] == False,
        df['nucleotide_mutation_count'],
        0
    )
    df['nonsynonymous_nucleotide_mutation_count'] = np.where(
        df['aa_mutation'] == True,
        df['nucleotide_mutation_count'],
        0
    )

    # Sum mutations per branch (pcp_index)
    df['synonymous_mutations_per_branch'] = df.groupby('pcp_index')['synonymous_nucleotide_mutation_count'].transform('sum')
    df['nonsynonymous_mutations_per_branch'] = df.groupby('pcp_index')['nonsynonymous_nucleotide_mutation_count'].transform('sum')
    df['total_mutations_per_branch'] = df['nonsynonymous_mutations_per_branch'] + df['synonymous_mutations_per_branch']

    # Calculate mutation frequencies (mutations / sequence length)
    df['synonymous_mutation_freq_branch'] = df['synonymous_mutations_per_branch'] / df['seq_nuc_length']
    df['nonsynonymous_mutation_freq_branch'] = df['nonsynonymous_mutations_per_branch'] / df['seq_nuc_length']
    df['total_mutation_freq_branch'] = df['total_mutations_per_branch'] / df['seq_nuc_length']

    # Drop intermediate columns that are not needed downstream
    df.drop(columns=[
        'seq_nuc_length',
        'aa_mutation',
        'synonymous_nucleotide_mutation_count',
        'nonsynonymous_nucleotide_mutation_count',
        'synonymous_mutations_per_branch',
        'nonsynonymous_mutations_per_branch',
        'total_mutations_per_branch',
    ], inplace=True)

    return df


def plot_dasm_vs_rates_comparison(compare_dasm_rates, entrenched_sites_aas, site_color_map,
                                   savefig_prefix=None, title="Comparison of Observed/Expected Rates Ratio vs DASM Selection Factor", 
                                   title_extra='', figures_dir='figures/'):
    """
    Create a scatter plot comparing observed/expected rate ratios to DASM selection factors.

    Uses orthogonal regression (appropriate when both variables have measurement error)
    and reports Pearson correlation.

    Parameters:
    -----------
    compare_dasm_rates : pd.DataFrame
        DataFrame with columns:
        - 'log_ratio': log(observed_counts / expected_counts)
        - 'log_selection_factor': log of DASM selection factor
        - 'site', 'v_family', 'parent_aa', 'child_aa': for merging with entrenched sites

    entrenched_sites_aas : pd.DataFrame
        DataFrame with columns: 'site', 'v_family', 'amino_acid', 'target_amino_acid'

    site_color_map : dict
        Dictionary mapping site names to colors for consistent plotting

    savefig_prefix : str, optional
        If provided, save figure with this prefix

    title_extra : str, optional
        Additional text to add to the plot title

    figures_dir : str, optional
        Directory to save figures (default: 'figures/')

    Returns:
    --------
    None (displays plot and optionally saves to file)
    """
    # Calculate regression statistics
    x = compare_dasm_rates['log_ratio']
    y = compare_dasm_rates['log_selection_factor']

    # Remove any NaN values for regression calculation
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]
    n = len(x_clean)

    # Calculate Pearson correlation
    r_value, p_value = stats.pearsonr(x_clean, y_clean)

    # Calculate orthogonal regression
    slope_ortho, intercept_ortho = orthogonal_regression(x_clean.values, y_clean.values)

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot regular points in grey
    sns.scatterplot(data=compare_dasm_rates,
                    x='log_ratio', y='log_selection_factor',
                    color='grey', alpha=0.3, label='Other sites')

    # Filter entrenched data
    entrenched_compare_rates_dasm = pd.merge(
        entrenched_sites_aas.rename(columns={'amino_acid': 'parent_aa', 'target_amino_acid': 'child_aa'}),
        compare_dasm_rates,
        on=['site', 'v_family', 'parent_aa', 'child_aa'],
        how='inner'
    )
    print(f"Plotting {len(entrenched_compare_rates_dasm)} entrenched points")

    # Print which points were not found in the compare_dasm_rates
    not_found = pd.merge(
        entrenched_sites_aas.rename(columns={'amino_acid': 'parent_aa', 'target_amino_acid': 'child_aa'}),
        compare_dasm_rates,
        on=['site', 'v_family', 'parent_aa', 'child_aa'],
        how='outer',
        indicator=True
    )
    not_found = not_found[not_found['_merge'] == 'left_only']
    if len(not_found) > 0:
        print("The following entrenched points were not found in the comparison data:")
        print(not_found[['site', 'v_family', 'parent_aa', 'child_aa']])

    # Plot entrenched points in color
    entrenched_compare_rates_dasm['site'] = entrenched_compare_rates_dasm['site'].astype(str)
    sns.scatterplot(data=entrenched_compare_rates_dasm,
                    x='log_ratio', y='log_selection_factor',
                    s=90, hue='site', style='v_family', palette=site_color_map)

    # Add orthogonal regression line
    x_range = np.array([x_clean.min(), x_clean.max()])
    y_ortho = slope_ortho * x_range + intercept_ortho
    ax.plot(x_range, y_ortho, linestyle='-', color='blue', linewidth=2, label='Orthogonal regression')

    ax.axvline(0, color='black', linestyle=':', linewidth=1)
    ax.axhline(0, color='black', linestyle=':', linewidth=1)

    # Add legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.xlabel('Observed Rate / Expected Rate (log)')
    plt.ylabel('DASM Selection Factor (log)')

    # Format the equation for the title
    if intercept_ortho >= 0:
        equation_ortho = f'y = {slope_ortho:.3f}x + {intercept_ortho:.3f}'
    else:
        equation_ortho = f'y = {slope_ortho:.3f}x - {abs(intercept_ortho):.3f}'

    title = f'{title}\nOrthogonal regression: {equation_ortho}, R² = {r_value**2:.3f}\nn = {n} {title_extra}'
    plt.title(title)

    plt.tight_layout()
    plt.show()

    if savefig_prefix:
        fig.savefig(f'{figures_dir}/{savefig_prefix}validation_dasm_vs_rates_comparison.png', dpi=300, bbox_inches='tight')


def plot_rates_pairwise_analysis(compare_dasm_rates, pairwise_df_dict, site_color_map,
                                  savefig_prefix=None, title_extra='', figures_dir='figures/'):
    """
    Create pairwise comparison plots of observed/expected count ratios across V families.

    Uses orthogonal regression (appropriate when both variables have measurement error)
    and reports R² from Pearson correlation.

    Parameters:
    -----------
    compare_dasm_rates : pd.DataFrame
        DataFrame with columns:
        - 'v_family', 'site', 'parent_aa', 'child_aa'
        - 'observed_counts', 'expected_counts', 'log_ratio'

    pairwise_df_dict : dict
        Dictionary mapping comparison names (e.g., 'IGHV1_vs_IGHV3') to DataFrames
        with entrenched site information

    site_color_map : dict
        Dictionary mapping site names to colors for consistent plotting

    savefig_prefix : str, optional
        If provided, save figure with this prefix

    title_extra : str, optional
        Additional text to add to the plot title

    figures_dir : str, optional
        Directory to save figures (default: 'figures/')

    Returns:
    --------
    None (displays plot and optionally saves to file)
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 7))
    fig.subplots_adjust(hspace=0.5, wspace=0.3, right=0.85, top=0.8)
    axes = axes.flatten()
    ax_i = 0

    # Collect all unique sites across all subplots
    all_sites = set()

    plot_order = ['IGHV1_vs_IGHV3', 'IGHV1_vs_IGHV4', 'IGHV3_vs_IGHV4', 'within_IGHV1', 'within_IGHV3', 'within_IGHV4']

    for cur_pair_name in plot_order:
        cur_pairwise_df = pairwise_df_dict[cur_pair_name]
        # Get entrenched sites
        cur_pairwise_df = cur_pairwise_df[cur_pairwise_df.is_entrenched == True]

        # Create pairwise log ratio dataframe
        compare_dasm_rates1 = compare_dasm_rates[['v_family', 'site', 'parent_aa', 'child_aa', 'observed_counts', 'expected_counts', 'log_ratio']].copy()
        compare_dasm_rates1 = compare_dasm_rates1.rename(columns={
            'parent_aa': 'parent_aa_1_and_target_aa_2',
            'child_aa': 'parent_aa_2_and_target_aa_1'
        })

        compare_dasm_rates2 = compare_dasm_rates[['v_family', 'site', 'parent_aa', 'child_aa', 'observed_counts', 'expected_counts', 'log_ratio']].copy()
        compare_dasm_rates2 = compare_dasm_rates2.rename(columns={
            'parent_aa': 'parent_aa_2_and_target_aa_1',
            'child_aa': 'parent_aa_1_and_target_aa_2'
        })

        counts_pairwise = pd.merge(
            compare_dasm_rates1, compare_dasm_rates2,
            on=['site', 'parent_aa_1_and_target_aa_2', 'parent_aa_2_and_target_aa_1'],
            suffixes=('_1', '_2')
        )

        # Filter counts_pairwise by the v_families being compared
        if cur_pair_name.startswith('within_'):
            # Within family comparison - both v_families should be the same
            v_fam = cur_pair_name.replace('within_', '')
            counts_pairwise = counts_pairwise[
                (counts_pairwise['v_family_1'] == v_fam) &
                (counts_pairwise['v_family_2'] == v_fam)
            ]
        else:
            # Between family comparison (e.g., IGHV1_vs_IGHV3)
            v_fam1, v_fam2 = cur_pair_name.split('_vs_')
            counts_pairwise = counts_pairwise[
                (counts_pairwise['v_family_1'] == v_fam1) &
                (counts_pairwise['v_family_2'] == v_fam2)
            ]

        # Merge with current entrenched sites according to DASM analysis
        entrenched_merged_pairwise = pd.merge(
            counts_pairwise, cur_pairwise_df,
            on=['site', 'parent_aa_1_and_target_aa_2', 'parent_aa_2_and_target_aa_1', 'v_family_1', 'v_family_2'],
            how='inner'
        )

        # Count how many points were not found in the merge
        not_found = pd.merge(
            cur_pairwise_df, counts_pairwise,
            on=['site', 'parent_aa_1_and_target_aa_2', 'parent_aa_2_and_target_aa_1', 'v_family_1', 'v_family_2'],
            how='outer',
            indicator=True
        )
        not_found = not_found[not_found['_merge'] == 'left_only']

        # Set subplot title
        base_title = cur_pair_name.replace('_', ' ')
        n_plotted = len(entrenched_merged_pairwise)
        n_total = len(cur_pairwise_df)
        base_title += f"\n(plotting {n_plotted} out of {n_total} entrenched site+aa pairs)"
        axes[ax_i].set_title(base_title, fontsize=10)

        # Collect all unique sites
        all_sites.update(entrenched_merged_pairwise['site'].unique())

        # Plot
        sns.scatterplot(counts_pairwise, x='log_ratio_1', y='log_ratio_2', color='grey',
                       ax=axes[ax_i], alpha=0.3, label='Other sites')
        sns.scatterplot(entrenched_merged_pairwise, x='log_ratio_1', y='log_ratio_2',
                       hue='site', palette=site_color_map, ax=axes[ax_i], s=90)

        # Add y=x dashed line
        ax_min = min(counts_pairwise['log_ratio_1'].min(), counts_pairwise['log_ratio_2'].min())
        ax_max = max(counts_pairwise['log_ratio_1'].max(), counts_pairwise['log_ratio_2'].max())
        axes[ax_i].plot([ax_min, ax_max], [ax_min, ax_max], linestyle='--', color='black', linewidth=1, label='y=x')

        # Remove individual subplot legends
        if axes[ax_i].get_legend():
            axes[ax_i].get_legend().remove()

        axes[ax_i].axvline(0, color='black', linestyle=':', linewidth=1)
        axes[ax_i].axhline(0, color='black', linestyle=':', linewidth=1)
        ax_i += 1

    # Add title
    fig.suptitle(f'Germline-divergent sites comparison of Observed/Expected Counts Ratios\n{title_extra}')

    # Create legend handles for all unique sites
    legend_handles = []
    legend_labels = []

    # Add "Other sites" first
    legend_handles.append(plt.scatter([], [], color='grey', alpha=0.3))
    legend_labels.append('Other sites')

    # Add y=x line legend entry
    legend_handles.append(plt.Line2D([0], [0], color='black', linestyle='--', linewidth=1))
    legend_labels.append('y=x')

    # Add all unique sites in sorted order
    sorted_sites = sort_antibody_sites(list(all_sites))
    for site in sorted_sites:
        legend_handles.append(plt.scatter([], [], color=site_color_map[site], s=90))
        legend_labels.append(site)

    # Create unified legend on the right side
    fig.legend(legend_handles, legend_labels, loc='center left', bbox_to_anchor=(0.87, 0.5),
               frameon=True, title='Sites')

    if savefig_prefix:
        fig.savefig(f'{figures_dir}/{savefig_prefix}validation_rates_pairwise_comparison.pdf', dpi=800)

    fig.show()
