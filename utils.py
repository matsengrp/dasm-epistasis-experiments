from tqdm import tqdm
import pandas as pd
import numpy as np

from netam.sequences import translate_sequence, AA_STR_SORTED
from netam.codon_table import single_mutant_aa_indices

from Bio.Data import CodonTable


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

    # Add the one_mutation_away column to the dataframe with progress bar
    tqdm.pandas(desc="Adding one_mutation_away column")
    df["one_mutation_away"] = df.progress_apply(
        lambda row: CODON_AA_MUTATION_DICT[(row[parent_codon_col], row[target_aa_col])],
        axis=1,
    )


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

    # Add the one_mutation_away column to the dataframe with progress bar
    tqdm.pandas(desc="Adding one_mutation_away column")
    df["one_mutation_away"] = df.progress_apply(
        lambda row: AA_AA_MUTATION_DICT[(row[parent_aa_col], row[target_aa_col])],
        axis=1,
    )


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


def add_germline_information(
    pcp_df, site_df, germline_codons_path="germline/germline_codons.csv"
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

    # Ensure site column compatibility
    germline_codons_df["site"] = germline_codons_df["site"].astype(float)
    site_df_enhanced["site"] = site_df_enhanced["site"].astype(float)

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
    )

    # Add germline codon comparison if parent_codon exists
    if "parent_codon" in site_df_with_germline.columns:
        site_df_with_germline["is_germline_codon"] = (
            site_df_with_germline.parent_codon == site_df_with_germline.germline_codon
        )

    return site_df_with_germline
