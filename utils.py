
from tqdm import tqdm
import pandas as pd

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



def add_column_aa_one_mutation_away_from_codon(df, parent_codon_col='parent_codon', target_aa_col='selection_factor_target_aa'):
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
    df['one_mutation_away'] = df.progress_apply(lambda row: CODON_AA_MUTATION_DICT[(row[parent_codon_col], row[target_aa_col])], axis=1)


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

def add_column_aa_one_mutation_away_from_aa(df, parent_aa_col='parent_aa', target_aa_col='selection_factor_target_aa'):
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
    df['one_mutation_away'] = df.progress_apply(lambda row: AA_AA_MUTATION_DICT[(row[parent_aa_col], row[target_aa_col])], axis=1)