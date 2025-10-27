import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio.Seq import Seq

from dnsmex.local import localify
from dnsmex.neutral_mutability import MutabilityContainer

from netam.oe_plot import annotate_sites_df
from netam.models import DEFAULT_NEUTRAL_MODEL




### create expected mutation counts using neutral model
### we will inherit fro the MutabilityContainer class from dnsmex, but will use branch length per sepcific pcp index to calcaulte the observed mutation counts instead of a fixed branch length of 0.1.
### has alternate branch length calculation methods, see doc below. We will likely use from_snynoymous_mutations

class NeutralMutationProbability(MutabilityContainer):
    """
    Variant of MutabilityContainer that uses per-sequence branch lengths from pcp_df.

    Unlike MutabilityContainer which uses a fixed branch length (default 0.1) for all sequences,
    this class reads branch lengths from the 'branch_length' column in pcp_df for each pcp_index,
    or calculates them from observed mutations.

    Branch length calculation methods:
    - 'from_tree': Use pre-existing 'branch_length' column (requires pcp_df to have this column)
    - 'from_synonymous_mutations': Calculate from observed synonymous mutation frequency
    - 'from_nonsynonymous_mutations': Calculate from observed nonsynonymous mutation frequency
    - 'from_total_mutations': Calculate from total observed mutation frequency

    Requirements:
    - For 'from_tree': pcp_df must have a 'branch_length' column
    - For mutation-based methods: pcp_df must have parent/child codon information

    All other functionality is identical to MutabilityContainer.
    """

    def __init__(
        self,
        dataset,
        neutral_model_name=DEFAULT_NEUTRAL_MODEL,
        subset_size=None,
        branch_length_method='from_tree',
    ):
        """
        Initialize NeutralMutationProbability.
        
        Parameters:
        -----------
        dataset : str
            Dataset name
        neutral_model_name : str
            Name of the neutral model to use
        subset_size : int, optional
            Number of sequences to use (for testing)
        branch_length_method : str
            Method for calculating branch lengths:
            - 'from_tree': Use existing branch_length column
            - 'from_synonymous_mutations': Calculate from synonymous mutations
            - 'from_nonsynonymous_mutations': Calculate from nonsynonymous mutations
            - 'from_total_mutations': Calculate from total mutations
        """
        # Validate branch_length_method
        valid_methods = ['from_tree', 'from_synonymous_mutations', 
                        'from_nonsynonymous_mutations', 'from_total_mutations']
        if branch_length_method not in valid_methods:
            raise ValueError(
                f"branch_length_method must be one of {valid_methods}, "
                f"got '{branch_length_method}'"
            )
        
        self.branch_length_method = branch_length_method
        
        # Call parent __init__ but don't pass branch_length parameter
        # (we'll override the method that uses it)
        super().__init__(
            dataset=dataset,
            neutral_model_name=neutral_model_name,
            branch_length=None,  # Not used in this subclass
            subset_size=subset_size,
        )

        # Validate requirements based on method
        if branch_length_method == 'from_tree':
            if 'branch_length' not in self.pcp_df.columns:
                raise ValueError(
                    "pcp_df must have a 'branch_length' column when using "
                    "branch_length_method='from_tree'. "
                    "Use simulation.add_branch_lengths_to_trees() or similar to add branch lengths."
                )

    def _calculate_branch_lengths_from_mutations(self, nuc_neutral_df):
        """
        Calculate branch lengths from observed mutations between parent and child codons.
        
        Returns a dataframe with columns: ['pcp_index', 'branch_length']
        """
        # Get parent/child codon information
        temp_site_data = nuc_neutral_df[['pcp_index', 'site',
                                         'parent_codon', 'parent_aa', 
                                         'child_codon', 'child_aa']].drop_duplicates()
        
        # Calculate sequence nucleotide length per pcp_index
        temp_site_data['seq_nuc_length'] = temp_site_data.groupby('pcp_index').transform('size') * 3
        
        # Count nucleotide mutations per codon
        temp_site_data['nucleotide_mutation_count'] = temp_site_data.apply(
            lambda row: sum(c1 != c2 for c1, c2 in zip(row['child_codon'], row['parent_codon'])), 
            axis=1
        )
        
        # Determine if mutation is synonymous or nonsynonymous
        temp_site_data['mutation'] = temp_site_data['parent_aa'] != temp_site_data['child_aa']
        
        # Count synonymous mutations
        temp_site_data['synonymous_nucleotide_mutation_count'] = np.where(
            temp_site_data['mutation'] == False,
            temp_site_data['nucleotide_mutation_count'],
            0
        )
        
        # Count nonsynonymous mutations
        temp_site_data['nonsynonymous_nucleotide_mutation_count'] = np.where(
            temp_site_data['mutation'] == True,
            temp_site_data['nucleotide_mutation_count'],
            0
        )
        
        # Aggregate per branch (pcp_index)
        temp_site_data['synonymous_mutations_per_branch'] = temp_site_data.groupby('pcp_index')[
            'synonymous_nucleotide_mutation_count'].transform('sum')
        temp_site_data['nonsynonymous_mutations_per_branch'] = temp_site_data.groupby('pcp_index')[
            'nonsynonymous_nucleotide_mutation_count'].transform('sum')
        temp_site_data['total_mutations_per_branch'] = (
            temp_site_data['nonsynonymous_mutations_per_branch'] + 
            temp_site_data['synonymous_mutations_per_branch']
        )
        
        # Calculate mutation frequencies (mutations per nucleotide)
        temp_site_data['synonymous_mutation_freq_branch'] = (
            temp_site_data['synonymous_mutations_per_branch'] / temp_site_data['seq_nuc_length']
        )
        temp_site_data['nonsynonymous_mutation_freq_branch'] = (
            temp_site_data['nonsynonymous_mutations_per_branch'] / temp_site_data['seq_nuc_length']
        )
        temp_site_data['total_mutation_freq_branch'] = (
            temp_site_data['total_mutations_per_branch'] / temp_site_data['seq_nuc_length']
        )
        
        # Get per-branch summary
        branch_data = temp_site_data[[
            'pcp_index', 'seq_nuc_length',
            'synonymous_mutations_per_branch', 'nonsynonymous_mutations_per_branch',
            'total_mutations_per_branch', 'synonymous_mutation_freq_branch',
            'nonsynonymous_mutation_freq_branch', 'total_mutation_freq_branch'
        ]].drop_duplicates()
        
        # Select the appropriate frequency based on method
        if self.branch_length_method == 'from_synonymous_mutations':
            branch_data['branch_length'] = branch_data['synonymous_mutation_freq_branch']
        elif self.branch_length_method == 'from_nonsynonymous_mutations':
            branch_data['branch_length'] = branch_data['nonsynonymous_mutation_freq_branch']
        elif self.branch_length_method == 'from_total_mutations':
            branch_data['branch_length'] = branch_data['total_mutation_freq_branch']
        
        return branch_data[['pcp_index', 'branch_length']]

    def create_nucleotide_neutral_rates_df(self):
        """
        Create a dataframe with one row per sequence, site, and nucleotide.
        Uses per-sequence branch lengths from pcp_df['branch_length'] or calculates them
        from observed mutations, depending on branch_length_method.
        """
        # Create the raw nucleotide data (same as parent)
        nuc_neutral_df = self._create_raw_nucleotide_data()

        nuc_neutral_df = nuc_neutral_df[~nuc_neutral_df.current_codon.isna()]

        # amino acid site annotation - create amino acid site wise annotations for every (pcp_index, site)
        temp_annotation_df = (
            nuc_neutral_df[["pcp_index", "site"]].drop_duplicates().copy()
        )
        temp_annotation_df["unannotated_nuc_site"] = temp_annotation_df["site"].copy()
        temp_annotation_df["aa_site"] = temp_annotation_df.site // 3
        annotated_temp_annotation_df = annotate_sites_df(
            temp_annotation_df[["pcp_index", "aa_site"]].drop_duplicates(),
            self.pcp_df,
            numbering_dict=self.numbering,
            add_codons_aas=True,
        )
        # merge back the annotated sites
        temp_annotation_df = pd.merge(
            annotated_temp_annotation_df,
            temp_annotation_df[["pcp_index", "aa_site", "unannotated_nuc_site"]],
            on=["pcp_index", "aa_site"],
            how="inner",
        )
        temp_annotation_df.drop(columns=["aa_site"], inplace=True)
        temp_annotation_df = temp_annotation_df[
            [
                "pcp_index",
                "site",
                "unannotated_nuc_site",
                "is_cdr",
                "parent_codon",
                "parent_aa",
                "child_codon",
                "child_aa",
            ]
        ].rename(columns={"unannotated_nuc_site": "nuc_site"})
        nuc_neutral_df = pd.merge(
            nuc_neutral_df.rename(columns={"site": "nuc_site"}),
            temp_annotation_df,
            on=["pcp_index", "nuc_site"],
            how="inner",
        )

        # Get branch lengths based on the specified method
        if self.branch_length_method == 'from_tree':
            # Use existing branch lengths from pcp_df
            branch_length_df = self.pcp_df[['branch_length']].reset_index().rename(
                columns={'index': 'pcp_index'}
            )
        else:
            # Calculate branch lengths from observed mutations
            branch_length_df = self._calculate_branch_lengths_from_mutations(nuc_neutral_df)
        
        # Merge branch lengths
        nuc_neutral_df = pd.merge(
            nuc_neutral_df,
            branch_length_df,
            on='pcp_index',
            how='left'
        )

        # Use per-sequence branch length to calculate probability
        nuc_neutral_df["probability"] = nuc_neutral_df.apply(
            lambda row: 1.0 - np.exp(-row['branch_length'] * row['rate']),
            axis=1
        )

        # get substitution probability of mutation to specific nucleotide
        nuc_neutral_df["substitution_probability"] = (
            nuc_neutral_df["probability"] * nuc_neutral_df["csp"]
        )

        # add amino acids transition information
        nuc_neutral_df["current_aa"] = nuc_neutral_df.apply(
            lambda row: (
                Seq(row["current_codon"]).translate()
                if row["current_codon"] is not None
                else None
            ),
            axis=1,
        )
        nuc_neutral_df["transition_aa"] = nuc_neutral_df.apply(
            lambda row: (
                Seq(row["transition_codon"]).translate()
                if row["transition_codon"] is not None
                else None
            ),
            axis=1,
        )

        # reorder columns (include branch_length in output for transparency)
        nuc_neutral_df = nuc_neutral_df[
            [
                "pcp_index",
                "nuc_site",
                "current_nucleotide",
                "transition_nucleotide",
                "rate",
                "csp",
                "branch_length",
                "probability",
                "substitution_probability",
                "codon_position",
                "current_codon",
                "transition_codon",
                "current_aa",
                "transition_aa",
                "site",
                "is_cdr",
                "parent_codon",
                "parent_aa",
                "child_codon",
                "child_aa",
            ]
        ]

        self.nuc_neutral_df = nuc_neutral_df


class CachedNeutralMutationProbability:
    """
    Cached version of NeutralMutationProbability that saves/loads dataframes to/from gzip-compressed CSV files.
    
    Usage:
    - If cache files exist, loads them instantly
    - If cache files don't exist, creates them using NeutralMutationProbability and saves to cache
    - Provides the same interface as NeutralMutationProbability
    """

    def __init__(
        self,
        dataset,
        neutral_model_name=DEFAULT_NEUTRAL_MODEL,
        subset_size=None,
        branch_length_method='from_tree',
        cache_dir=localify(f"ANALYSIS_CACHE/"),
    ):
        """
        Initialize CachedNeutralMutationProbability.
        
        Parameters:
        -----------
        dataset : str
            Dataset name
        neutral_model_name : str
            Name of the neutral model to use
        subset_size : int, optional
            Number of sequences to use (for testing)
        branch_length_method : str
            Method for calculating branch lengths:
            - 'from_tree': Use existing branch_length column
            - 'from_synonymous_mutations': Calculate from synonymous mutations
            - 'from_nonsynonymous_mutations': Calculate from nonsynonymous mutations
            - 'from_total_mutations': Calculate from total mutations
        cache_dir : str
            Directory to store cache files
        """
        self.neutral_model_name = neutral_model_name
        self.dataset = dataset
        self.subset_size = subset_size
        self.branch_length_method = branch_length_method
        self.cache_dir = cache_dir

        # Create cache directory if it doesn't exist
        os.makedirs(cache_dir, exist_ok=True)

        # Generate cache file paths with subset and branch_length_method indication
        cache_parts = [f"neutral_mutation_prob", neutral_model_name, dataset]
        
        if subset_size:
            cache_parts.append(f"subset{subset_size}")
        
        # Add branch length method to cache name (short form for readability)
        method_abbrev = {
            'from_tree': 'tree',
            'from_synonymous_mutations': 'syn',
            'from_nonsynonymous_mutations': 'nonsyn',
            'from_total_mutations': 'total'
        }
        cache_parts.append(f"bl_{method_abbrev[branch_length_method]}")
        
        cache_base = f"{cache_dir}/{'_'.join(cache_parts)}"

        self.cache_files = {
            "nucleotide": f"{cache_base}_nucleotide.csv.gz",
            "amino_acid": f"{cache_base}_amino_acid.csv.gz",
            "amino_acid_to_any": f"{cache_base}_amino_acid_to_any.csv.gz",
            "codon": f"{cache_base}_codon.csv.gz",
            "codon_to_any": f"{cache_base}_codon_to_any.csv.gz",
            "pcp_df": f"{cache_base}_pcp_df.csv.gz",
        }

        # Try to load from cache, otherwise create new data
        if self._cache_exists():
            self._load_from_cache()
        else:
            self._create_and_cache()

    def _cache_exists(self):
        """Check if all required cache files exist"""
        return all(os.path.exists(filepath) for filepath in self.cache_files.values())

    def _load_from_cache(self):
        """Load dataframes from gzip-compressed cache files"""
        subset_msg = f" (subset {self.subset_size})" if self.subset_size else ""
        method_msg = f" [method: {self.branch_length_method}]"
        print(f"Loading NeutralMutationProbability data from gzip cache{subset_msg}{method_msg}...")

        try:
            self.nuc_neutral_df = pd.read_csv(
                self.cache_files["nucleotide"], compression="gzip"
            )
            self.aa_neutral_df = pd.read_csv(
                self.cache_files["amino_acid"], compression="gzip"
            )
            self.aa_to_any_neutral_df = pd.read_csv(
                self.cache_files["amino_acid_to_any"], compression="gzip"
            )
            self.codon_neutral_df = pd.read_csv(
                self.cache_files["codon"], compression="gzip"
            )
            self.codon_to_any_neutral_df = pd.read_csv(
                self.cache_files["codon_to_any"], compression="gzip"
            )
            self.pcp_df = pd.read_csv(self.cache_files["pcp_df"], compression="gzip")

            # Validate that DataFrames are not empty
            if any(
                len(df) == 0
                for df in [
                    self.nuc_neutral_df,
                    self.aa_neutral_df,
                    self.aa_to_any_neutral_df,
                    self.codon_neutral_df,
                    self.codon_to_any_neutral_df,
                    self.pcp_df,
                ]
            ):
                raise ValueError("One or more cached DataFrames are empty")

            print(f"✓ Loaded from gzip cache:")
            print(f"  - Nucleotide DataFrame: {len(self.nuc_neutral_df):,} rows")
            print(f"  - Amino Acid DataFrame: {len(self.aa_neutral_df):,} rows")
            print(
                f"  - Amino Acid to Any DataFrame: {len(self.aa_to_any_neutral_df):,} rows"
            )
            print(f"  - Codon DataFrame: {len(self.codon_neutral_df):,} rows")
            print(
                f"  - Codon to Any DataFrame: {len(self.codon_to_any_neutral_df):,} rows"
            )
            print(f"  - PCP DataFrame: {len(self.pcp_df):,} rows")

        except Exception as e:
            print(f"Error loading from gzip cache: {e}")
            print("Creating new data...")
            self._create_and_cache()

    def _create_and_cache(self):
        """Create new NeutralMutationProbability and save to cache"""
        subset_msg = f" (subset {self.subset_size})" if self.subset_size else ""
        method_msg = f" [method: {self.branch_length_method}]"
        print(f"Creating new NeutralMutationProbability data{subset_msg}{method_msg}...")

        # Create the container using NeutralMutationProbability
        container = NeutralMutationProbability(
            self.dataset, 
            self.neutral_model_name, 
            subset_size=self.subset_size,
            branch_length_method=self.branch_length_method
        )

        # Copy the dataframes
        self.nuc_neutral_df = container.nuc_neutral_df.copy()
        self.aa_neutral_df = container.aa_neutral_df.copy()
        self.aa_to_any_neutral_df = container.aa_to_any_neutral_df.copy()
        self.codon_neutral_df = container.codon_neutral_df.copy()
        self.codon_to_any_neutral_df = container.codon_to_any_neutral_df.copy()
        self.pcp_df = container.pcp_df.copy()

        # Copy other attributes that might be needed
        self.rates = container.rates
        self.csp_logits = container.csp_logits
        self.csp_rates = container.csp_rates
        self.numbering = container.numbering
        if hasattr(container, "numbering_type"):
            self.numbering_type = container.numbering_type

        # Save to cache
        self._save_to_cache()

    def _save_to_cache(self):
        """Save dataframes to gzip-compressed cache files"""
        print("Saving data to gzip cache...")

        try:
            self.nuc_neutral_df.to_csv(
                self.cache_files["nucleotide"], index=False, compression="gzip"
            )
            self.aa_neutral_df.to_csv(
                self.cache_files["amino_acid"], index=False, compression="gzip"
            )
            self.aa_to_any_neutral_df.to_csv(
                self.cache_files["amino_acid_to_any"], index=False, compression="gzip"
            )
            self.codon_neutral_df.to_csv(
                self.cache_files["codon"], index=False, compression="gzip"
            )
            self.codon_to_any_neutral_df.to_csv(
                self.cache_files["codon_to_any"], index=False, compression="gzip"
            )
            self.pcp_df.to_csv(
                self.cache_files["pcp_df"], index=False, compression="gzip"
            )

            print(f"✓ Saved to gzip cache:")
            print(f"  - {self.cache_files['nucleotide']}")
            print(f"  - {self.cache_files['amino_acid']}")
            print(f"  - {self.cache_files['amino_acid_to_any']}")
            print(f"  - {self.cache_files['codon']}")
            print(f"  - {self.cache_files['codon_to_any']}")
            print(f"  - {self.cache_files['pcp_df']}")

        except Exception as e:
            print(f"Error saving to gzip cache: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate neutral mutation probabilities for a dataset')
    parser.add_argument('dataset', type=str, help='Dataset name (e.g., v1tang, v1rodriguez)')
    parser.add_argument('--branch-length-method', type=str, default='from_total_mutations',
                        choices=['from_tree', 'from_synonymous_mutations',
                                'from_nonsynonymous_mutations', 'from_total_mutations'],
                        help='Method for calculating branch lengths (default: from_total_mutations)')

    args = parser.parse_args()

    neutral_probabilties = CachedNeutralMutationProbability(
        args.dataset,
        branch_length_method=args.branch_length_method
    )
