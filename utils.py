from tqdm import tqdm
import pandas as pd
import numpy as np
import os

from netam.sequences import translate_sequence, AA_STR_SORTED
from netam.codon_table import single_mutant_aa_indices

from Bio.Data import CodonTable

from dnsmex.local import localify
from dnsmex import dasm_oe, dnsm_oe, dasm_zoo, dnsm_zoo

GERMLINE_PATH_DICTIONARY = {'imgt':'germline/germline_codons_imgt.csv', 'chothia':'germline/germline_codons_chothia.csv'}


def get_cdr_definitions(numbering_scheme='imgt', chain='heavy'):
    """
    Get CDR region definitions based on numbering scheme and chain type.

    Parameters:
    -----------
    numbering_scheme : str
        Either 'imgt' or 'chothia' (default: 'imgt')
    chain : str
        Either 'heavy' or 'light' (default: 'heavy')
        Only affects Chothia numbering scheme

    Returns:
    --------
    list of tuples : List of (start, end) tuples for each CDR region
        [(cdr1_start, cdr1_end), (cdr2_start, cdr2_end), (cdr3_start, cdr3_end)]

    Notes:
    ------
    CDR boundaries:
    - IMGT (same for heavy and light): CDR1 (27-38), CDR2 (56-65), CDR3 (105-117)
    - Chothia Heavy: CDR1 (26-32), CDR2 (52-56), CDR3 (95-102)
    - Chothia Light: CDR1 (24-34), CDR2 (50-56), CDR3 (89-97)

    Example:
    --------
    >>> cdr_regions = get_cdr_definitions('imgt')
    """
    if numbering_scheme == 'imgt':
        # IMGT CDR boundaries (same for heavy and light chains)
        return [(27, 38), (56, 65), (105, 117)]
    elif numbering_scheme == 'chothia':
        if chain == 'heavy':
            # Chothia heavy chain CDR boundaries
            return [(26, 32), (52, 56), (95, 102)]
        elif chain == 'light':
            # Chothia light chain CDR boundaries
            return [(24, 34), (50, 56), (89, 97)]
        else:
            raise ValueError(f"Invalid chain: {chain}. Must be 'heavy' or 'light'.")
    else:
        raise ValueError(f"Invalid numbering_scheme: {numbering_scheme}. Must be 'imgt' or 'chothia'.")


def is_in_cdr(site, numbering_scheme='imgt', chain='heavy'):
    """
    Check if a site is in any CDR region.

    Parameters:
    -----------
    site : int or str
        Site number to check (can be int like 30 or str like '30A' for Chothia insertions)
    numbering_scheme : str
        Either 'imgt' or 'chothia' (default: 'imgt')
    chain : str
        Either 'heavy' or 'light' (default: 'heavy')

    Returns:
    --------
    bool : True if site is in a CDR region, False otherwise

    Example:
    --------
    >>> is_in_cdr(30, 'imgt')  # CDR1 in IMGT
    True
    >>> is_in_cdr('30A', 'chothia')  # CDR1 in Chothia with insertion
    True
    >>> is_in_cdr(50, 'imgt')  # Framework region
    False
    """
    # Extract base number from site (handles both int and str with insertion codes)
    site_str = str(site)

    # Handle IMGT decimal notation (e.g., '30.1', '30.2')
    if '.' in site_str:
        base_num = int(site_str.split('.')[0])
    else:
        # Handle Chothia letter notation (e.g., '30A', '30B') or plain numbers
        num_part = ''.join(c for c in site_str if c.isdigit())
        base_num = int(num_part) if num_part else 0

    cdr_regions = get_cdr_definitions(numbering_scheme, chain)
    return any(start <= base_num <= end for start, end in cdr_regions)


def create_codon_aa_mutation_dict():
    """Create a dictionary mapping (parent_codon, target_aa) -> boolean"""

    standard_table = CodonTable.unambiguous_dna_by_id[1]
    all_codons = list(standard_table.forward_table.keys())

    mutation_dict = {}

    for parent_codon in all_codons:
        # Get all amino acids reachable by single mutation
        reachable_aa_indices = single_mutant_aa_indices(parent_codon)
        reachable_aas = {AA_STR_SORTED[i] for i in reachable_aa_indices}

        # For each possible amino acid, set True/False
        for aa in AA_STR_SORTED:
            mutation_dict[(parent_codon, aa)] = aa in reachable_aas

    return mutation_dict


def add_column_aa_one_mutation_away_from_codon(
    df, parent_codon_col="parent_codon", target_aa_col="selection_factor_target_aa"
):
    """
    Add a column to the DataFrame indicating if the target amino acid is one mutation away from the parent codon.

    Parameters:
    - df: DataFrame with columns 'parent_codon' and 'selection_factor_target_aa'.

    Returns:
    - DataFrame with an additional column 'one_mutation_away'.
    """

    # Create the dictionary
    CODON_AA_MUTATION_DICT = create_codon_aa_mutation_dict()

    # Vectorized lookup using map - much faster than apply
    print("Adding one_mutation_away column (vectorized)...")
    df["one_mutation_away"] = df[[parent_codon_col, target_aa_col]].apply(
        lambda row: (row[parent_codon_col], row[target_aa_col]), axis=1
    ).map(CODON_AA_MUTATION_DICT)


def create_aa_aa_mutation_dict():
    """Create a dictionary mapping (parent_aa, target_aa) -> boolean"""

    standard_table = CodonTable.unambiguous_dna_by_id[1]
    all_codons = list(standard_table.forward_table.keys())

    mutation_dict = {}

    # First, collect all reachable amino acids for each amino acid
    aa_reachable = {}

    for parent_codon in all_codons:
        parent_aa = standard_table.forward_table[parent_codon]

        # Get all amino acids reachable by single mutation from this codon
        reachable_aa_indices = single_mutant_aa_indices(parent_codon)
        reachable_aas = {AA_STR_SORTED[i] for i in reachable_aa_indices}

        # Add to the set of reachable amino acids for this parent amino acid
        if parent_aa not in aa_reachable:
            aa_reachable[parent_aa] = set()
        aa_reachable[parent_aa].update(reachable_aas)

    # Create the final dictionary mapping (parent_aa, target_aa) -> boolean
    for parent_aa in AA_STR_SORTED:
        reachable_aas = aa_reachable.get(parent_aa, set())

        for target_aa in AA_STR_SORTED:
            mutation_dict[(parent_aa, target_aa)] = target_aa in reachable_aas

    return mutation_dict


def add_column_aa_one_mutation_away_from_aa(
    df, parent_aa_col="parent_aa", target_aa_col="selection_factor_target_aa"
):
    """
    Add a column to the DataFrame indicating if the target amino acid is one mutation away from the parent codon.

    Parameters:
    - df: DataFrame with columns 'parent_aa' and 'selection_factor_target_aa'.

    Returns:
    - DataFrame with an additional column 'one_mutation_away'.
    """

    # Create the dictionary
    AA_AA_MUTATION_DICT = create_aa_aa_mutation_dict()

    # Vectorized lookup using map - much faster than apply
    print("Adding one_mutation_away column (vectorized)...")
    df["one_mutation_away"] = df[[parent_aa_col, target_aa_col]].apply(
        lambda row: (row[parent_aa_col], row[target_aa_col]), axis=1
    ).map(AA_AA_MUTATION_DICT)


def sort_antibody_sites(sites):
    """
    Sort antibody numbering sites properly, handling both IMGT and Chothia insertion codes.

    IMGT uses decimal notation: 30, 30.1, 30.2, 31
    Chothia uses letter notation: 30, 30A, 30B, 31

    Args:
        sites: Iterable of site identifiers (can be int, float, or str)

    Returns:
        List of sites sorted in proper antibody numbering order

    Examples:
        >>> sort_antibody_sites(['1', '1A', '2', '3'])
        ['1', '1A', '2', '3']
        >>> sort_antibody_sites([1, 1.1, 1.2, 2])
        [1, 1.1, 1.2, 2]
        >>> sort_antibody_sites(['30', '30A', '30B', '31', '31A'])
        ['30', '30A', '30B', '31', '31A']
    """
    def site_key(site):
        site_str = str(site)

        # Handle IMGT decimal notation (e.g., '30.1', '30.2')
        if '.' in site_str:
            parts = site_str.split('.')
            base_num = int(parts[0])
            decimal_part = int(parts[1]) if len(parts) > 1 else 0
            # Return: (base number, decimal part, empty string for letter)
            # This ensures 30 < 30.1 < 30.2 < 31
            return (base_num, decimal_part, '')

        # Handle Chothia letter notation (e.g., '30A', '30B')
        # Split into number and letter parts
        num_part = ''
        letter_part = ''
        for char in site_str:
            if char.isdigit():
                num_part += char
            else:
                letter_part += char

        base_num = int(num_part) if num_part else 0
        # Return: (base number, 0 for no decimal, letter part)
        # Empty letter sorts before any letter: 30 < 30A < 30B < 31
        return (base_num, 0, letter_part)

    return sorted(sites, key=site_key)


def add_cdr_shading(ax, sorted_sites, cdr_regions=None, numbering_scheme='imgt',
                    chain='heavy', color='red', alpha=0.1):
    """
    Add CDR background shading to a matplotlib axis with string site labels.

    This function handles the conversion from site numbers to x-axis positions,
    supporting both IMGT (decimal notation) and Chothia (letter notation) schemes.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis to add shading to
    sorted_sites : list
        List of site identifiers in the order they appear on the x-axis
        (e.g., ['1', '2', '27', '27A', '28', ...])
    cdr_regions : list of tuples, optional
        List of (start, end) site numbers for CDR regions
        If None, uses default regions based on numbering_scheme and chain
    numbering_scheme : str
        Either 'imgt' or 'chothia' (default: 'imgt')
    chain : str
        Either 'heavy' or 'light' (default: 'heavy')
        Only affects Chothia numbering scheme
    color : str
        Color for the shading (default: 'red')
    alpha : float
        Transparency of the shading (default: 0.1)

    Returns:
    --------
    None (modifies ax in place)

    Notes:
    ------
    Default CDR boundaries:
    - IMGT (same for heavy and light): CDR1 (27-38), CDR2 (56-65), CDR3 (105-117)
    - Chothia Heavy: CDR1 (26-32), CDR2 (52-56), CDR3 (95-102)
    - Chothia Light: CDR1 (24-34), CDR2 (50-56), CDR3 (89-97)

    Example:
    --------
    >>> fig, ax = plt.subplots()
    >>> sorted_sites = sort_antibody_sites(df['site'].unique())
    >>> ax.scatter(range(len(sorted_sites)), values)
    >>> ax.set_xticks(range(len(sorted_sites)))
    >>> ax.set_xticklabels(sorted_sites, rotation=90)
    >>> add_cdr_shading(ax, sorted_sites, numbering_scheme='chothia', chain='heavy')
    """
    if cdr_regions is None:
        cdr_regions = get_cdr_definitions(numbering_scheme, chain)

    # Create a mapping from site to x-position (index in sorted sites)
    site_to_position = {site: i for i, site in enumerate(sorted_sites)}

    # Add shading for each CDR region
    for cdr_start, cdr_end in cdr_regions:
        # Find sites within this CDR region
        # Extract numeric part of site for comparison
        cdr_sites = []
        for site in sorted_sites:
            site_str = str(site)
            # Extract base number (handles both '30.1' and '30A' formats)
            if '.' in site_str:
                base_num = int(site_str.split('.')[0])
            else:
                # Remove any letters to get base number
                num_part = ''.join(c for c in site_str if c.isdigit())
                base_num = int(num_part) if num_part else 0

            if cdr_start <= base_num <= cdr_end:
                cdr_sites.append(site)

        # Shade the region if sites were found
        if cdr_sites:
            start_pos = site_to_position[cdr_sites[0]]
            end_pos = site_to_position[cdr_sites[-1]]
            ax.axvspan(start_pos, end_pos, color=color, alpha=alpha)


def load_and_process_dasm_data(
    model_name,
    dataset_name,
    numbering_scheme='imgt',
    figures_dir=None,
    force_recompute=False
):
    """
    Load or compute DASM data with germline annotations and additional features.

    This function handles the complete workflow of:
    1. Loading cached CSV files if available (or computing if not)
    2. Running write_sites_oe if data needs to be computed
    3. Adding germline annotations
    4. Adding tree/clonal family information
    5. Adding one mutation away annotations
    6. Computing log selection factors

    Parameters:
    -----------
    model_name : str
        Name of the trained DASM model (e.g., "dasm_4m-v1jaffeCC+v1tangCC-joint")
    dataset_name : str
        Name of the dataset to analyze (e.g., "v1rodriguez_short")
    numbering_scheme : str
        Either 'imgt' or 'chothia'
    figures_dir : str, optional
        Directory for output figures. If None, no figure path is generated
    force_recompute : bool
        If True, recompute even if cached files exist

    Returns:
    --------
    tuple of (site_sub_probs_df, pcp_df, aa_site_subs_selection_df_germline)
        - site_sub_probs_df: Site-level substitution probabilities
        - pcp_df: Parent-child pair information
        - aa_site_subs_selection_df_germline: Amino acid substitution selection factors
          with germline annotations, tree info, and derived features

    Example:
    --------
    >>> site_df, pcp_df, aa_df = load_and_process_dasm_data(
    ...     model_name="dasm_4m-v1jaffeCC+v1tangCC-joint",
    ...     dataset_name="v1rodriguez_short",
    ...     numbering_scheme='chothia',
    ...     figures_dir="figures"
    ... )
    """
    crepe_prefix = localify(f"DASM_TRAINED_MODELS_DIR/{model_name}")
    test_output_prefix = localify(f"DASM_TEST_OUTPUT_DIR/{model_name}-ON-{dataset_name}")

    # Define file paths
    site_sub_probs_path = f"{test_output_prefix}-{numbering_scheme}-site_sub_probs_df.csv"
    pcp_path = f"{test_output_prefix}-{numbering_scheme}-pcp_df.csv"
    aa_site_subs_path = f"{test_output_prefix}-{numbering_scheme}-aa_site_subs_selection_df.csv"

    # Try to load cached data
    if not force_recompute:
        try:
            site_sub_probs_df = pd.read_csv(site_sub_probs_path, dtype={'site': str})
            pcp_df = pd.read_csv(pcp_path)
            aa_site_subs_selection_df = pd.read_csv(aa_site_subs_path, dtype={'site': str})

        except FileNotFoundError:
            force_recompute = True

    # Compute if needed
    if force_recompute:
        # Ensure branch lengths exist
        branch_length_path = f"{test_output_prefix}.branch_lengths_csv"
        if not os.path.exists(branch_length_path):
            dasm_zoo.write_branch_lengths(
                localify(f"DASM_TRAINED_MODELS_DIR/{model_name}"),
                dataset_name,
                branch_length_path
            )

        # Prepare figure output path
        fig_out_path = None
        if figures_dir is not None:
            fig_out_path = f"{figures_dir}/sites-oe-V1,3,4.svg"

        # Run write_sites_oe
        complete_plotter, plotter_dict = dasm_oe.write_sites_oe(
            crepe_prefix=crepe_prefix,
            dataset_name=dataset_name,
            branch_length_path=branch_length_path,
            csv_output_path=f"{test_output_prefix}-sites-oe.csv",
            fig_out_path=fig_out_path,
            min_log_prob=-4,
            replace_title=True,
            numbering_scheme=numbering_scheme,
        )
        complete_plotter = complete_plotter["heavy"]

        # Extract dataframes
        pcp_df = complete_plotter.pcp_df.copy()
        site_sub_probs_df = complete_plotter.site_sub_probs_df.copy()
        aa_site_subs_selection_df = complete_plotter.aa_site_subs_selection_df.copy()

        # Save to CSV
        pcp_df.to_csv(pcp_path, index=False)
        site_sub_probs_df.to_csv(site_sub_probs_path, index=False)
        aa_site_subs_selection_df.to_csv(aa_site_subs_path, index=False)

    # Add tree/clonal family information
    pcp_df_for_merge = pcp_df.copy()
    pcp_df_for_merge = pcp_df_for_merge[['sample_id', 'family', 'depth', 'distance']]
    pcp_df_for_merge['pcp_index'] = pcp_df_for_merge.index

    aa_site_subs_selection_df = pd.merge(
        aa_site_subs_selection_df,
        pcp_df_for_merge,
        on='pcp_index',
        how='inner'
    )

    # Add log selection factor
    aa_site_subs_selection_df['log_selection_factor'] = np.log(
        aa_site_subs_selection_df.selection_factor
    )

    # Add germline information
    aa_site_subs_selection_df_germline = add_germline_information(
        pcp_df,
        aa_site_subs_selection_df,
        numbering_scheme=numbering_scheme
    )

    # Add one mutation away annotation
    add_column_aa_one_mutation_away_from_codon(aa_site_subs_selection_df_germline)

    # Add cdr annotation
    aa_site_subs_selection_df_germline['is_cdr'] = aa_site_subs_selection_df_germline['site'].apply(is_in_cdr, numbering_scheme=numbering_scheme)

    return site_sub_probs_df, pcp_df, aa_site_subs_selection_df_germline


def load_and_process_dnsm_data(
    model_name,
    dataset_name,
    numbering_scheme='imgt',
    figures_dir=None,
    force_recompute=False
):
    """
    Load or compute DNSM data with germline annotations and additional features.

    This function handles the complete workflow for DNSM models:
    1. Loading cached CSV files if available (or computing if not)
    2. Running write_sites_oe if data needs to be computed
    3. Converting tensor values to scalars
    4. Adding germline annotations
    5. Adding tree/clonal family information
    6. Computing log selection factors

    Parameters:
    -----------
    model_name : str
        Name of the trained DNSM model (e.g., "dnsm_1m-v1jaffe+v1tang-joint")
    dataset_name : str
        Name of the dataset to analyze (e.g., "v1rodriguez")
    numbering_scheme : str
        Either 'imgt' or 'chothia'
    figures_dir : str, optional
        Directory for output figures. If None, no figure path is generated
    force_recompute : bool
        If True, recompute even if cached files exist

    Returns:
    --------
    tuple of (site_sub_probs_df, pcp_df)
        - site_sub_probs_df: Site-level substitution probabilities with germline annotations
        - pcp_df: Parent-child pair information

    Example:
    --------
    >>> site_df, pcp_df = load_and_process_dnsm_data(
    ...     model_name="dnsm_1m-v1jaffe+v1tang-joint",
    ...     dataset_name="v1rodriguez",
    ...     numbering_scheme='imgt',
    ...     figures_dir="figures"
    ... )
    """
    crepe_prefix = localify(f"DNSM_TRAINED_MODELS_DIR/{model_name}")
    test_output_prefix = localify(f"DNSM_TEST_OUTPUT_DIR/{model_name}-ON-{dataset_name}")

    # Define file paths
    site_sub_probs_path = f"{test_output_prefix}-{numbering_scheme}-site_sub_probs_df.csv"
    pcp_path = f"{test_output_prefix}-{numbering_scheme}-pcp_df.csv"

    # Try to load cached data
    if not force_recompute:
        try:
            site_sub_probs_df = pd.read_csv(site_sub_probs_path, dtype={'site': str})
            pcp_df = pd.read_csv(pcp_path)

        except FileNotFoundError:
            force_recompute = True

    # Compute if needed
    if force_recompute:
        # Ensure branch lengths exist
        branch_length_path = f"{test_output_prefix}.branch_lengths_csv"
        if not os.path.exists(branch_length_path):
            dnsm_zoo.write_branch_lengths(
                localify(f"DNSM_TRAINED_MODELS_DIR/{model_name}"),
                dataset_name,
                branch_length_path
            )

        # Prepare figure output path
        fig_out_path = None
        if figures_dir is not None:
            fig_out_path = f"{figures_dir}/sites-oe-V1,3,4.svg"

        # Run write_sites_oe for DNSM
        complete_plotter, plotter_dict = dnsm_oe.write_sites_oe(
            crepe_prefix=crepe_prefix,
            dataset_name=dataset_name,
            branch_length_path=branch_length_path,
            csv_output_path=f"{test_output_prefix}-sites-oe.csv",
            fig_out_path=fig_out_path,
            min_log_prob=-4,
            replace_title=True,
            numbering_scheme=numbering_scheme,
        )
        complete_plotter = complete_plotter["heavy"]

        # Extract dataframes
        pcp_df = complete_plotter.pcp_df.copy()
        site_sub_probs_df = complete_plotter.site_sub_probs_df.copy()

        # Convert tensor values to scalars (DNSM-specific)
        for column_name in ['selection_factor', 'neutral_prob', 'prob']:
            if column_name in site_sub_probs_df.columns:
                site_sub_probs_df[column_name] = site_sub_probs_df[column_name].apply(
                    lambda x: x.item() if hasattr(x, 'item') else x
                )

        # Save to CSV
        pcp_df.to_csv(pcp_path, index=False)
        site_sub_probs_df.to_csv(site_sub_probs_path, index=False)


    # Add log selection factor
    site_sub_probs_df['log_selection_factor'] = np.log(
        site_sub_probs_df['selection_factor']
    )

    # Add j_family annotation
    pcp_df['j_family'] = pcp_df['j_gene'].str.split('*').str[0]

    # Add tree/clonal family information
    pcp_df_for_merge = pcp_df.copy()
    pcp_df_for_merge = pcp_df_for_merge[[
        'j_gene', 'j_family',
        'sample_id', 'family', 'depth', 'distance'
    ]]
    pcp_df_for_merge['pcp_index'] = pcp_df_for_merge.index

    site_sub_probs_df = pd.merge(
        site_sub_probs_df,
        pcp_df_for_merge,
        on='pcp_index',
        how='inner'
    )

    # Add germline information
    site_sub_probs_df_germline = add_germline_information(
        pcp_df,
        site_sub_probs_df,
        numbering_scheme=numbering_scheme
    )

    # Add cdr annotation
    site_sub_probs_df_germline['is_cdr'] = site_sub_probs_df_germline['site'].apply(is_in_cdr, numbering_scheme=numbering_scheme)


    return site_sub_probs_df_germline, pcp_df


def add_germline_information(
    pcp_df, site_df, numbering_scheme='imgt'
):
    """
    Add germline information and family annotations to a site-level DataFrame.

    Parameters:
    -----------
    pcp_df : pd.DataFrame
        DataFrame with phylogenetic information containing columns:
        - v_gene, j_gene, v_family, sample_id, family, depth, distance
        - Index will be used as pcp_index

    site_df : pd.DataFrame
        DataFrame with site-level data containing at minimum:
        - pcp_index, site (IMGT aligned), parent_aa
        - Optionally: parent_codon

    germline_codons_path : str
        Path to germline_codons.csv file

    Returns:
    --------
    pd.DataFrame
        Enhanced DataFrame with germline information and family annotations
    """
    try:
        germline_codons_path = GERMLINE_PATH_DICTIONARY[numbering_scheme]
    except KeyError:
        raise ValueError(f"Invalid numbering_scheme: {numbering_scheme}. Must be 'imgt' or 'chothia'.")

    # Create a copy to avoid modifying original DataFrames
    site_df_enhanced = site_df.copy()
    pcp_df_work = pcp_df.copy()

    # Prepare merge DataFrame with family annotations
    pcp_df_for_merge = pcp_df_work[
        [
            "v_gene",
            "j_gene",
            "v_family",
        ]
    ].copy()
    pcp_df_for_merge["pcp_index"] = pcp_df_for_merge.index

    # Merge family information
    site_df_enhanced = pd.merge(
        site_df_enhanced, pcp_df_for_merge, on="pcp_index", how="inner"
    )

    germline_codons_df = pd.read_csv(germline_codons_path)

    # Ensure site column compatibility - keep as string to support both IMGT and Chothia
    # IMGT uses numeric sites, Chothia uses insertion codes like '52A'
    germline_codons_df["site"] = germline_codons_df["site"].astype(str)
    site_df_enhanced["site"] = site_df_enhanced["site"].astype(str)

    # Prepare germline columns for merge
    germline_for_merge = germline_codons_df.rename(
        columns={"codon": "germline_codon", "amino_acid": "germline_amino_acid"}
    )

    # Remove v_family from germline if it exists (to avoid conflicts)
    if "v_family" in germline_for_merge.columns:
        germline_for_merge = germline_for_merge.drop(columns=["v_family"])

    # Merge germline information
    site_df_with_germline = pd.merge(
        site_df_enhanced, germline_for_merge, on=["v_gene", "site"], how="left"
    )

    # Add germline comparison columns
    site_df_with_germline["is_germline_aa"] = (
        site_df_with_germline.parent_aa == site_df_with_germline.germline_amino_acid
    ).where(
        ~(site_df_with_germline.germline_amino_acid.isna())
    )

    # Add germline codon comparison if parent_codon exists
    if "parent_codon" in site_df_with_germline.columns:
        site_df_with_germline["is_germline_codon"] = (
            site_df_with_germline.parent_codon == site_df_with_germline.germline_codon
        ).where(
            ~(site_df_with_germline.germline_codon.isna())
        )

    return site_df_with_germline
