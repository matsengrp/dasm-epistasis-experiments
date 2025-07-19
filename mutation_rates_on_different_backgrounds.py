import itertools
import os
import numpy as np
import pandas as pd
import tqdm
from statsmodels.stats.rates import test_poisson_2indep
import argparse

from dnsmex import dxsm_data, dnsm_zoo
from dnsmex.dnsm_oe import write_sites_oe
from dnsmex.local import localify

figures_dir = localify("FIGURES_DIR")





def compare_mutation_rates_on_different_backgrounds(site_sub_probs_df_germline, site1, site2):
    '''
    Calculate the association between mutations at two sites, accounting for branch length.
    
    This function analyzes whether the presence of a mutation at site1 affects the rate
    of acquiring mutations at site2, using branch length as a measure of evolutionary time.
    
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
    
    site1 : int
        First site position (predictor variable)
    site2 : int  
        Second site position (outcome variable)
        
    Returns:
    --------
    dict
        Dictionary containing:
        - Site identifiers and summary statistics
        - Branch counts, lengths, and mutation rates for each group
        - Rate ratios and log rate ratios (both branch length and mutation count normalized)
        - Poisson test statistics and p-values (or NaN if insufficient data)
        
    Notes:
    ------
    - Only analyzes branches where site2 is still at germline (can acquire mutation)
    - Handles zero division explicitly: rates can be 0, inf, or finite values
    - Rate ratios can be 0, inf, or NaN depending on the data
    - Statistical tests return NaN if insufficient data for reliable inference
    - Removes branches where both sites are mutated on the same pcp, as the acquisition order cannot be resolved
    '''

    # filter only rows with the two relevant sites
    cur_df = site_sub_probs_df_germline[(site_sub_probs_df_germline['site'].isin([site1, site2]))].copy()
    # filter only rows where germline identity is known
    cur_df = cur_df[(cur_df['germline_amino_acid'].notna())]
    # merge data for both sites per pcp
    merged_df = pd.merge(cur_df[cur_df.site == site1], cur_df[cur_df.site == site2], on=['sample_id', 'family', 'pcp_index', 'branch_length', 'mutations_per_branch'], how='inner', suffixes=('_1', '_2'))
    # only keep pcps that have measurements for both sites
    merged_df.dropna(subset=['site_1', 'site_2'], inplace=True)
    # remove pcps/branches that already have mutation at site2, and thus cannot acquire it
    merged_df = merged_df[merged_df.is_germline_aa_2 == True]
    percent_of_unresolved_pcps = len(merged_df[((merged_df.mutation_1 == True) & (merged_df.mutation_2 == True))]) / len(merged_df) * 100 if not merged_df.empty else 0
    # remove branches that have a mutation at site1 and site2 on the same pcp - this means we cannot resolve who whether the mutation at site2 was acquired after site1 or not
    merged_df = merged_df[~((merged_df.mutation_1 == True) & (merged_df.mutation_2 == True))]

    # Calculate totals from actual data
    total_counts = len(merged_df[['family', 'sample_id', 'pcp_index']].drop_duplicates())
    total_branch_length = merged_df['branch_length'].sum()
    total_mutations_acquired = len(merged_df[(merged_df['mutation_2'] == True)])

    # Calculations:
    #                        total branch length     number of B mutation acquired  |
    # Has A mutation                                                                |  rate = mutations per total length
    # No A mutation                                                                 |  rate = mutations per total length

    # has A mutation - calculate actual values
    has_site1_mutation = merged_df[(merged_df.is_germline_aa_1 == False)]
    length_has_site1_mutation = has_site1_mutation['branch_length'].sum()
    mutcount_length_has_site1_mutation = has_site1_mutation['mutations_per_branch'].sum()
    has_site1_mutation_site2_acquired = len(has_site1_mutation[(has_site1_mutation['mutation_2'] == True)])

    # no A mutation - calculate actual values
    no_site1_mutation = merged_df[(merged_df.is_germline_aa_1 == True)]
    length_no_site1_mutation = no_site1_mutation['branch_length'].sum()
    mutcount_length_no_site1_mutation = no_site1_mutation['mutations_per_branch'].sum()
    no_site1_mutation_site2_acquired = len(no_site1_mutation[(no_site1_mutation['mutation_2'] == True)])

    # assertions to ensure data integrity
    assert total_counts == len(has_site1_mutation) + len(no_site1_mutation), \
        f"Mismatch in counts: {len(has_site1_mutation)} + {len(no_site1_mutation)} != {total_counts}"
    assert np.isclose(length_has_site1_mutation + length_no_site1_mutation, total_branch_length, atol=1e-5), \
        f"Mismatch in branch lengths: {length_has_site1_mutation} + {length_no_site1_mutation} != {total_branch_length}"    
    assert has_site1_mutation_site2_acquired + no_site1_mutation_site2_acquired == total_mutations_acquired, \
        f"Mismatch in mutations acquired: {has_site1_mutation_site2_acquired} + {no_site1_mutation_site2_acquired} != {total_mutations_acquired}"

    # Calculate mutation rates per branch length - handle zero division explicitly
    if length_has_site1_mutation > 0:
        rate_has_site1 = has_site1_mutation_site2_acquired / length_has_site1_mutation
    else: rate_has_site1 = 0.0  # No branch length means no mutations

    if length_no_site1_mutation > 0:
        rate_no_site1 = no_site1_mutation_site2_acquired / length_no_site1_mutation
    else: rate_no_site1 = 0.0  # No branch length means no mutations

    # Calculate rate ratio - handle edge cases
    if rate_no_site1 > 0 and np.isfinite(rate_no_site1):
        rate_ratio = rate_has_site1 / rate_no_site1
    elif rate_has_site1 > 0 and rate_no_site1 == 0: rate_ratio = np.inf  # Positive rate in numerator, zero in denominator
    else: rate_ratio = np.nan  # 0/0, inf/inf, or other indeterminate forms

    # Calculate log rate ratio
    if rate_ratio > 0 and np.isfinite(rate_ratio):
        log_rate_ratio = np.log(rate_ratio)
    elif rate_ratio == 0: log_rate_ratio = -np.inf
    elif rate_ratio == np.inf: log_rate_ratio = np.inf
    else: log_rate_ratio = np.nan

    # Calculate mutation rates per total mutation count - handle zero division explicitly
    if mutcount_length_has_site1_mutation > 0:
        rate_has_site1_mutcount = has_site1_mutation_site2_acquired / mutcount_length_has_site1_mutation
    else: rate_has_site1_mutcount = 0.0  # No mutations means rate is 0

    if mutcount_length_no_site1_mutation > 0:
        rate_no_site1_mutcount = no_site1_mutation_site2_acquired / mutcount_length_no_site1_mutation
    else: rate_no_site1_mutcount = 0.0  # No mutations means rate is 0

    # Calculate mutation count rate ratio - handle edge cases
    if rate_no_site1_mutcount > 0 and np.isfinite(rate_no_site1_mutcount):
        rate_ratio_mutcount = rate_has_site1_mutcount / rate_no_site1_mutcount
    elif rate_has_site1_mutcount > 0 and rate_no_site1_mutcount == 0: rate_ratio_mutcount = np.inf
    else: rate_ratio_mutcount = np.nan

    # Calculate log mutation count rate ratio
    if rate_ratio_mutcount > 0 and np.isfinite(rate_ratio_mutcount):
        log_rate_ratio_mutcount = np.log(rate_ratio_mutcount)
    elif rate_ratio_mutcount == 0: log_rate_ratio_mutcount = -np.inf
    elif rate_ratio_mutcount == np.inf: log_rate_ratio_mutcount = np.inf
    else: log_rate_ratio_mutcount = np.nan

       # Statsmodels Poisson test - only run if we have sufficient data
    min_count_threshold = 1  # Minimum count for reliable Poisson test
    min_exposure_threshold = 1e-10  # Minimum exposure for reliable Poisson test

    if (length_has_site1_mutation >= min_exposure_threshold and
        length_no_site1_mutation >= min_exposure_threshold):
        try:
            poisson_result = test_poisson_2indep(
                count1=has_site1_mutation_site2_acquired,
                exposure1=length_has_site1_mutation,
                count2=no_site1_mutation_site2_acquired,
                exposure2=length_no_site1_mutation,)
            poisson_statistic = poisson_result.statistic
            poisson_pvalue = poisson_result.pvalue
        except (ValueError, ZeroDivisionError):
            poisson_statistic = np.nan
            poisson_pvalue = np.nan
    else:
        poisson_statistic = np.nan
        poisson_pvalue = np.nan

    # Poisson test for mutation count version
    if (mutcount_length_has_site1_mutation >= min_exposure_threshold and
        mutcount_length_no_site1_mutation >= min_exposure_threshold):
        try:
            poisson_result_mutcount = test_poisson_2indep(
                count1=has_site1_mutation_site2_acquired,
                exposure1=mutcount_length_has_site1_mutation,
                count2=no_site1_mutation_site2_acquired,
                exposure2=mutcount_length_no_site1_mutation,)
            poisson_statistic_mutcount = poisson_result_mutcount.statistic
            poisson_pvalue_mutcount = poisson_result_mutcount.pvalue
        except (ValueError, ZeroDivisionError):
            poisson_statistic_mutcount = np.nan
            poisson_pvalue_mutcount = np.nan
    else:
        poisson_statistic_mutcount = np.nan
        poisson_pvalue_mutcount = np.nan

    # Create results dictionary with actual observed values
    results = {
        'site1': site1,
        'site2': site2,
        'total_branch_counts': total_counts,
        'total_branch_length': total_branch_length,
        'total_mutations_acquired': total_mutations_acquired,
        'percent_of_removed_unresolved_pcps': percent_of_unresolved_pcps,
        
        # Has site1 mutation group
        'has_site1_mutation_branch_count': len(has_site1_mutation),
        'has_site1_mutation_branch_length': length_has_site1_mutation,
        'has_site1_mutation_site2_acquired': has_site1_mutation_site2_acquired,
        'has_site1_mutation_rate': rate_has_site1,
        
        # No site1 mutation group
        'no_site1_mutation_branch_count': len(no_site1_mutation),
        'no_site1_mutation_branch_length': length_no_site1_mutation,
        'no_site1_mutation_site2_acquired': no_site1_mutation_site2_acquired,
        'no_site1_mutation_rate': rate_no_site1,
        
        # Rate comparison
        'rate_ratio': rate_ratio,
        'log_rate_ratio': log_rate_ratio,
        'rate_ratio_mutcount': rate_ratio_mutcount,
        'log_rate_ratio_mutcount': log_rate_ratio_mutcount,
        
        # Statistical test results
        'poisson_statistic': poisson_statistic,
        'poisson_pvalue': poisson_pvalue,
        'poisson_statistic_mutcount': poisson_statistic_mutcount,
        'poisson_pvalue_mutcount': poisson_pvalue_mutcount,
    }
    
    return results

def compare_mutation_rates_on_different_backgrounds_for_all_sites(site_sub_probs_df_germline, output_path):
    """
    Calculate mutation association for all pairs of sites in the provided list.
    
    Args:
        site_sub_probs_df_germline (pd.DataFrame): DataFrame with site substitution probabilities and germline info
        output_path (str): Output path to save the results CSV file        
    Returns:
        pd.DataFrame: DataFrame with results for each site pair
    """
    # prep site_sub_probs_df_germline with mutation counts per branch

    # run for all site combinations
    all_results = []
    combinations = list(itertools.permutations(site_sub_probs_df_germline.site.unique().tolist(), 2))
    for site1, site2 in tqdm.tqdm(combinations):
        result = compare_mutation_rates_on_different_backgrounds(site_sub_probs_df_germline, site1, site2)
        all_results.append(result)

    results_df = pd.DataFrame(all_results)

    print(f'Saving to {output_path}')
    results_df.to_csv(output_path, index=False)






def load_or_create_datafile(dataset_name='v1rodriguez'):
    model_name = "dnsm_1m-v1jaffe+v1tang-joint"

    crepe_prefix = localify(f"DNSM_TRAINED_MODELS_DIR/{model_name}")
    test_output_prefix = localify(f"DNSM_TEST_OUTPUT_DIR/{model_name}-ON-{dataset_name}")

    try:
        site_sub_probs_df = pd.read_csv(f"{test_output_prefix}-site_sub_probs_df.csv")
        pcp_df = pd.read_csv(f"{test_output_prefix}-pcp_df.csv")


    except FileNotFoundError:
        if not os.path.exists(f"{test_output_prefix}.branch_lengths_csv"):
            dnsm_zoo.write_branch_lengths("/fh/fast/matsen_e/shared/bcr-mut-sel/dnsm/dnsm-experiments-1/dnsm-train/trained_models/dnsm_1m-v1jaffe+v1tang-joint", dataset_name, f"{test_output_prefix}.branch_lengths_csv")


        complete_plotter, plotter_dict = write_sites_oe(
            crepe_prefix=crepe_prefix,
            dataset_name=dataset_name,
            branch_length_path=f"{test_output_prefix}.branch_lengths_csv",
            csv_output_path=f"{test_output_prefix}-sites-oe.csv",
            fig_out_path=f"{figures_dir}/sites-oe-V1,3,4.svg",
            min_log_prob=-4,
            replace_title=True,
        )
        complete_plotter = complete_plotter["heavy"]


        pcp_df = complete_plotter.pcp_df.copy()
        site_sub_probs_df = complete_plotter.site_sub_probs_df.copy()
        for column_name in ['selection_factor', 'neutral_prob', 'prob']:
            site_sub_probs_df[column_name] = site_sub_probs_df[column_name].apply(lambda x: x.item() if hasattr(x, 'item') else x)


        pcp_df.to_csv(f"{test_output_prefix}-pcp_df.csv", index=False)
        site_sub_probs_df.to_csv(f"{test_output_prefix}-site_sub_probs_df.csv", index=False)


    # add family annotations to aa and site substitution dfs
    pcp_df['j_family'] = pcp_df['j_gene'].str.split('*').str[0]
    pcp_df_for_merge = pcp_df.copy()
    pcp_df_for_merge = pcp_df_for_merge[['v_gene', 'j_gene', 'v_family', 'j_family', 'sample_id', 'family', 'distance', 'branch_length', 'parent_name', 'child_name']]
    pcp_df_for_merge['pcp_index'] = pcp_df_for_merge.index

    site_sub_probs_df = pd.merge(site_sub_probs_df, pcp_df_for_merge, on='pcp_index', how='inner')

    site_sub_probs_df['log_selection_factor'] = np.log(site_sub_probs_df['selection_factor'])


    ## add germline information
    germline_codons_df = pd.read_csv(localify(f"DATA_DIR/germline_codons.csv"))
    germline_codons_df['site'] = germline_codons_df['site'].astype(float)
    site_sub_probs_df_germline = pd.merge(site_sub_probs_df, germline_codons_df.rename(columns={'codon':'germline_codon', 'amino_acid':'germline_amino_acid'}).drop(columns=['v_family']), on=['v_gene', 'site'], how='left')
    site_sub_probs_df_germline['is_germline_aa'] = site_sub_probs_df_germline.parent_aa == site_sub_probs_df_germline.germline_amino_acid
    site_sub_probs_df_germline['is_germline_codon'] = site_sub_probs_df_germline.parent_codon == site_sub_probs_df_germline.germline_codon

    # count mutations per branch as an alternative to branch length
    site_sub_probs_df_germline['mutations_per_branch'] = site_sub_probs_df_germline.groupby('pcp_index')['mutation'].transform('sum')
    return site_sub_probs_df_germline




def main():
    parser = argparse.ArgumentParser(description='Calculate mutation associations')
    parser.add_argument('inputs', nargs='+', help='Input files or parameters')
    
    args = parser.parse_args()
    
    print(f"Received {len(args.inputs)} inputs: {args.inputs}")
    
    all_site_sub_probs_df_germline = []
    for input_file in args.inputs:
         site_sub_probs_df_germline = load_or_create_datafile(dataset_name=input_file)
         all_site_sub_probs_df_germline.append(site_sub_probs_df_germline)
    
    site_sub_probs_df_germline = pd.concat(all_site_sub_probs_df_germline, ignore_index=True)

    # Create output filename from input datasets
    input_names = "_".join(args.inputs)
    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    output_path = localify(f"DATA_DIR/epistasis/mutation_rates_on_different_backgrounds/mutation_rates_on_different_backgrounds_{input_names}_{timestamp}.csv")
    

    compare_mutation_rates_on_different_backgrounds_for_all_sites(site_sub_probs_df_germline, output_path)


if __name__ == "__main__":
    main()