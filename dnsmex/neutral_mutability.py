import pandas as pd
import numpy as np
from Bio.Seq import Seq
import os
import torch
from torch.nn import functional as F
from tqdm import tqdm

from netam import pretrained
from netam.models import DEFAULT_NEUTRAL_MODEL
from netam.dasm import DASMDataset
from netam.sequences import BASES, CODONS
from dnsmex.dxsm_data import pcp_df_of_nickname, get_anarci_dict
from netam.oe_plot import annotate_sites_df, get_numbering_dict
from dnsmex.local import localify

# Pre-compute codon-to-amino-acid lookup for fast translation
CODON_TO_AA = {codon: str(Seq(codon).translate()) for codon in CODONS}


def nt_mutation_frequency(parent, child):
    """Return the fraction of nucleotide sites that differ between parent and child.

    Excludes positions where either parent or child has 'N'.
    """
    return sum(
        1 for p, c in zip(parent, child) if p != c and p != "N" and c != "N"
    ) / len(parent)


def synonymous_mutation_frequency(parent, child):
    """Return the fraction of synonymous mutations between parent and child sequences.

    A synonymous mutation is a nucleotide change that does not change the amino acid.
    Excludes codons containing 'N' in either parent or child.

    Args:
        parent: Parent nucleotide sequence (must be multiple of 3)
        child: Child nucleotide sequence (must be multiple of 3)

    Returns:
        Fraction of nucleotide positions with synonymous mutations
    """
    if len(parent) % 3 != 0 or len(child) % 3 != 0:
        raise ValueError("Sequences must be multiples of 3 for codon analysis")
    if len(parent) != len(child):
        raise ValueError("Parent and child must be same length")

    synonymous_count = 0
    total_positions = 0

    # Iterate over codons
    for i in range(0, len(parent), 3):
        parent_codon = parent[i : i + 3]
        child_codon = child[i : i + 3]

        # Skip codons with N
        if "N" in parent_codon or "N" in child_codon:
            continue

        # Translate codons to amino acids
        parent_aa = str(Seq(parent_codon).translate())
        child_aa = str(Seq(child_codon).translate())

        # Check each position in the codon
        for j in range(3):
            total_positions += 1
            if parent_codon[j] != child_codon[j]:
                # This position has a mutation
                if parent_aa == child_aa:
                    # The mutation is synonymous
                    synonymous_count += 1

    if total_positions == 0:
        return 0.0

    return synonymous_count / total_positions


class NeutralMutabilityDataset(DASMDataset):
    """DASMDataset extended with neutral mutation probability DataFrame methods.

    Provides analysis-friendly DataFrames with IMGT/Chothia numbering for:
    - Nucleotide-level neutral mutation probabilities
    - Codon-level neutral mutation probabilities
    - Amino acid-level neutral mutation probabilities

    This class uses DASMDataset's computation of neutral mutation
    probabilities via molevol.neutral_codon_probs_of_seq(), which correctly
    includes (1-p) terms for non-mutating positions in each codon.

    Distinct from DASMOEPlotter which analyzes selection factors from trained models.
    """

    def __init__(
        self,
        *args,
        numbering_dict=None,
        pcp_df=None,
        numbering_scheme=None,
        dataset_nickname=None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.numbering_dict = numbering_dict
        self.pcp_df = pcp_df  # Store for annotations
        self.numbering_scheme = numbering_scheme
        self.dataset_nickname = dataset_nickname

    @classmethod
    def of_pcp_df_with_numbering(
        cls,
        dataset_nickname: str,
        branch_length_mode: str = "constant",
        branch_length: float = 0.1,
        branch_length_scale_factor: float = 1.0,
        numbering_scheme: str = "imgt",
        multihit_model=None,
        subset_size=None,
    ):
        """Create NeutralMutabilityDataset from dataset nickname with IMGT numbering.

        Args:
            dataset_nickname: Dataset name (e.g., 'v1rodriguez')
            branch_length_mode: One of 'constant', 'mutation_frequency', or 'synonymous_frequency'
            branch_length: Used if mode='constant' (default: 0.1)
            branch_length_scale_factor: Multiplicative factor for all branch lengths
            numbering_scheme: Numbering scheme to use (default: 'imgt')
            multihit_model: Optional multihit model (default: None)
            subset_size: Optional subset size for testing (default: None)
        """
        # Validate branch_length_mode
        valid_modes = ["constant", "mutation_frequency", "synonymous_frequency"]
        if branch_length_mode not in valid_modes:
            raise ValueError(
                f"branch_length_mode must be one of {valid_modes}, got '{branch_length_mode}'"
            )
        # Always infer model_known_token_count from neutral model
        neutral_crepe = pretrained.load(DEFAULT_NEUTRAL_MODEL, device="cpu")
        model_known_token_count = neutral_crepe.encoder.site_count
        print(
            f"Inferred model_known_token_count={model_known_token_count} from {DEFAULT_NEUTRAL_MODEL}"
        )

        # Load pcp_df with rates and CSPs already computed!
        pcp_df = pcp_df_of_nickname(
            dataset_nickname,
            add_shm_outputs=True,  # This adds nt_rates_heavy, nt_csps_heavy columns
            sample_count=subset_size,
        )

        # Use parent class constructor - handles everything!
        # Note: We use branch_length_scale_factor directly here, not multiplied by 5
        # The factor of 5 is only used during training to calibrate mutation frequency-based branch lengths
        dataset = super(NeutralMutabilityDataset, cls).of_pcp_df(
            pcp_df,
            model_known_token_count=model_known_token_count,
            branch_length_multiplier=branch_length_scale_factor,
            multihit_model=multihit_model,
        )

        # Calculate branch lengths based on mode
        # Note: pcp_df uses "_heavy" suffix for column names at this point
        if branch_length_mode == "constant":
            # Use fixed branch length for all sequences
            dataset.branch_lengths = torch.full(
                (len(pcp_df),), branch_length * branch_length_scale_factor
            )
            dataset.update_neutral_probs()  # Recompute with new branch lengths
        elif branch_length_mode == "mutation_frequency":
            # Use total nucleotide mutation frequency
            branch_lengths = np.array(
                [
                    nt_mutation_frequency(parent, child) * branch_length_scale_factor
                    for parent, child in zip(
                        pcp_df["parent_heavy"], pcp_df["child_heavy"]
                    )
                ]
            )
            dataset.branch_lengths = torch.tensor(branch_lengths)
            dataset.update_neutral_probs()  # Recompute with new branch lengths
        elif branch_length_mode == "synonymous_frequency":
            # Use synonymous mutation frequency only
            # Add small epsilon to avoid zero branch lengths (required by DXSMDataset)
            epsilon = 1e-6
            branch_lengths = np.array(
                [
                    max(synonymous_mutation_frequency(parent, child), epsilon)
                    * branch_length_scale_factor
                    for parent, child in zip(
                        pcp_df["parent_heavy"], pcp_df["child_heavy"]
                    )
                ]
            )
            dataset.branch_lengths = torch.tensor(branch_lengths)
            dataset.update_neutral_probs()  # Recompute with new branch lengths

        # Add numbering dict and pcp_df for DataFrame methods
        anarci_path = get_anarci_dict(numbering_scheme)[dataset_nickname]

        # Create simplified column names for get_numbering_dict compatibility
        pcp_df_heavy = pcp_df.copy()
        for colname in pcp_df.columns:
            if colname.endswith("_heavy"):
                shortened_colname = colname[: -len("_heavy")]
                pcp_df_heavy[shortened_colname] = pcp_df_heavy[colname]

        dataset.numbering_dict, _ = get_numbering_dict(
            anarci_path["heavy"], pcp_df_heavy, verbose=True, checks=numbering_scheme
        )
        dataset.pcp_df = pcp_df_heavy
        dataset.numbering_scheme = numbering_scheme
        dataset.dataset_nickname = dataset_nickname

        return dataset

    def get_nucleotide_probabilities_df(self):
        """Create DataFrame with nucleotide-level mutation probabilities.

        Returns DataFrame with same structure as MutabilityContainer.nuc_neutral_df.
        Note: This computes marginal nucleotide mutation propensities, NOT the
        probability of a single-nucleotide change in a codon.
        """
        n_sequences = self.nt_ratess.shape[0]
        n_sites = self.nt_ratess.shape[1]
        rows = []

        for seq_idx in tqdm(
            range(n_sequences), desc="Processing sequences to nucleotide-level rates"
        ):
            parent_seq = self.pcp_df.iloc[seq_idx]["parent"]
            seq_length = len(parent_seq)

            for site_idx in range(n_sites):
                rate = self.nt_ratess[seq_idx, site_idx].item()

                # Only process sites that exist in the sequence
                if site_idx < seq_length:
                    # Get codon information for this site
                    codon_start = (site_idx // 3) * 3
                    codon_position = site_idx % 3

                    # Make sure we don't go beyond sequence length for codon
                    if codon_start + 2 < seq_length:
                        current_codon = parent_seq[codon_start : codon_start + 3]
                    else:
                        current_codon = None

                    current_nucleotide = parent_seq[site_idx]
                else:
                    # Sites beyond sequence length - model padding
                    current_codon = None
                    codon_position = None
                    current_nucleotide = None

                for nuc_idx, transition_nucleotide in enumerate(BASES):
                    csp_rate = self.nt_cspss[seq_idx, site_idx, nuc_idx].item()

                    # Determine the transition being tested
                    if current_nucleotide is not None and current_codon is not None:
                        # Calculate the resulting codon after transition
                        codon_list = list(current_codon)
                        codon_list[codon_position] = transition_nucleotide
                        transition_codon = "".join(codon_list)
                    else:
                        transition_codon = None

                    rows.append(
                        {
                            "pcp_index": seq_idx,
                            "site": site_idx,
                            "transition_nucleotide": transition_nucleotide,
                            "rate": rate,
                            "csp": csp_rate,
                            "current_codon": current_codon,
                            "codon_position": codon_position,
                            "current_nucleotide": current_nucleotide,
                            "transition_codon": transition_codon,
                        }
                    )

        nuc_neutral_df = pd.DataFrame(rows)
        nuc_neutral_df = nuc_neutral_df[~nuc_neutral_df.current_codon.isna()]

        # amino acid site annotation
        temp_annotation_df = (
            nuc_neutral_df[["pcp_index", "site"]].drop_duplicates().copy()
        )
        temp_annotation_df["unannotated_nuc_site"] = temp_annotation_df["site"].copy()
        temp_annotation_df["aa_site"] = temp_annotation_df.site // 3
        annotated_temp_annotation_df = annotate_sites_df(
            temp_annotation_df[["pcp_index", "aa_site"]].drop_duplicates(),
            self.pcp_df,
            numbering_dict=self.numbering_dict,
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

        # change rate to probability with per-sequence branch length (vectorized)
        # Create a mapping from pcp_index to branch_length for fast lookup
        branch_length_map = {
            idx: self.branch_lengths[idx].item()
            for idx in nuc_neutral_df["pcp_index"].unique()
        }
        nuc_neutral_df["probability"] = 1.0 - np.exp(
            -nuc_neutral_df["pcp_index"].map(branch_length_map) * nuc_neutral_df["rate"]
        )

        # get substitution probability of mutation to specific nucleotide (vectorized)
        # Special case: when current_nucleotide == transition_nucleotide,
        # substitution_probability should be 1 - probability (probability of NOT mutating)
        is_same_nucleotide = (
            nuc_neutral_df["current_nucleotide"]
            == nuc_neutral_df["transition_nucleotide"]
        )
        nuc_neutral_df["substitution_probability"] = np.where(
            is_same_nucleotide,
            1.0 - nuc_neutral_df["probability"],  # probability of NOT mutating
            nuc_neutral_df["probability"] * nuc_neutral_df["csp"],  # probability * csp
        )

        # add amino acids transition information using fast lookup table
        nuc_neutral_df["current_aa"] = nuc_neutral_df["current_codon"].map(CODON_TO_AA)
        nuc_neutral_df["transition_aa"] = nuc_neutral_df["transition_codon"].map(
            CODON_TO_AA
        )

        # reorder columns to match original structure
        nuc_neutral_df = nuc_neutral_df[
            [
                "pcp_index",
                "nuc_site",
                "current_nucleotide",
                "transition_nucleotide",
                "rate",
                "csp",
                "probability",
                "substitution_probability",
                "codon_position",
                "current_codon",
                "transition_codon",
                "current_aa",
                "transition_aa",
                "site",
                "parent_codon",
                "parent_aa",
                "child_codon",
                "child_aa",
            ]
        ]

        return nuc_neutral_df

    def get_codon_probabilities_df(self):
        """Create DataFrame with codon-level mutation probabilities.

        Uses self.log_neutral_codon_probss which is correctly computed with (1-p) terms!
        Returns DataFrame with same structure as MutabilityContainer.codon_neutral_df.

        Results are cached in _codon_probabilities_df for performance.
        """
        # Return cached result if available
        if hasattr(self, "_codon_probabilities_df"):
            return self._codon_probabilities_df

        # Convert log probabilities to linear space
        codon_probs = torch.exp(
            self.log_neutral_codon_probss
        )  # (n_seq, max_codon_len, 64)

        n_sequences = codon_probs.shape[0]
        max_codon_len = codon_probs.shape[1]
        rows = []

        for seq_idx in tqdm(range(n_sequences), desc="Processing codon probabilities"):
            parent_seq = self.pcp_df.iloc[seq_idx]["parent"]
            seq_length = len(parent_seq)
            n_codons = seq_length // 3

            for codon_idx in range(min(n_codons, max_codon_len)):
                codon_start = codon_idx * 3
                if codon_start + 2 < seq_length:
                    current_codon_str = parent_seq[codon_start : codon_start + 3]

                    # Skip if codon contains N
                    if "N" in current_codon_str:
                        continue

                    current_aa = CODON_TO_AA[current_codon_str]

                    # Iterate through all possible target codons
                    for target_codon_idx, target_codon_str in enumerate(CODONS):
                        prob = codon_probs[seq_idx, codon_idx, target_codon_idx].item()

                        # Only include if probability is non-zero
                        if prob > 0:
                            target_aa = CODON_TO_AA[target_codon_str]

                            # Only include single-nucleotide mutations
                            hamming_dist = sum(
                                c1 != c2
                                for c1, c2 in zip(current_codon_str, target_codon_str)
                            )
                            if hamming_dist == 1:
                                rows.append(
                                    {
                                        "pcp_index": seq_idx,
                                        "aa_site": codon_idx,
                                        "current_codon": current_codon_str,
                                        "transition_codon": target_codon_str,
                                        "current_aa": current_aa,
                                        "transition_aa": target_aa,
                                        "substitution_probability": prob,
                                    }
                                )

        codon_neutral_df = pd.DataFrame(rows)

        # Add annotations using annotate_sites_df
        temp_annotation_df = annotate_sites_df(
            codon_neutral_df[["pcp_index", "aa_site"]].drop_duplicates(),
            self.pcp_df,
            numbering_dict=self.numbering_dict,
            add_codons_aas=True,
        )

        # Merge annotations
        codon_neutral_df = pd.merge(
            codon_neutral_df,
            temp_annotation_df,
            on=["pcp_index", "aa_site"],
            how="inner",
        )

        # Reorder columns to match original structure
        codon_neutral_df = codon_neutral_df[
            [
                "pcp_index",
                "site",
                "current_codon",
                "transition_codon",
                "current_aa",
                "transition_aa",
                "parent_codon",
                "parent_aa",
                "child_codon",
                "child_aa",
                "substitution_probability",
            ]
        ]

        # Cache the result for future calls
        self._codon_probabilities_df = codon_neutral_df

        return codon_neutral_df

    def get_aa_probabilities_df(self):
        """Create DataFrame with amino acid-level mutation probabilities.

        Aggregates codon probabilities to AA level.
        Returns DataFrame with same structure as MutabilityContainer.aa_neutral_df.
        """
        # Get codon probabilities first
        codon_df = self.get_codon_probabilities_df()

        # Aggregate by amino acid transitions
        aa_neutral_df = (
            codon_df.groupby(
                [
                    "pcp_index",
                    "site",
                    "current_aa",
                    "transition_aa",
                    "parent_codon",
                    "parent_aa",
                    "child_codon",
                    "child_aa",
                ],
                sort=False,
            )
            .agg({"substitution_probability": "sum"})
            .reset_index()
        )

        return aa_neutral_df

    def get_aa_to_any_probabilities_df(self):
        """Create DataFrame with non-synonymous amino acid mutation probabilities.

        Aggregates all non-synonymous mutations for each amino acid.
        Returns DataFrame with same structure as MutabilityContainer.aa_to_any_neutral_df.
        """
        # Get AA probabilities first
        aa_df = self.get_aa_probabilities_df()

        # Remove synonymous mutations
        temp_aa_neutral_df = aa_df[aa_df.current_aa != aa_df.transition_aa].copy()

        # Capture original pcp_index and site order
        order_df = aa_df[["pcp_index", "site"]].drop_duplicates()
        order_df["_order"] = range(len(order_df))

        # Aggregate to any non-synonymous mutation
        aa_to_any_neutral_df = (
            temp_aa_neutral_df.groupby(
                [
                    "pcp_index",
                    "site",
                    "current_aa",
                    "parent_codon",
                    "parent_aa",
                    "child_codon",
                    "child_aa",
                ]
            )
            .agg({"substitution_probability": "sum"})
            .reset_index()
        )

        # Sort by original pcp_index and site order
        aa_to_any_neutral_df = (
            aa_to_any_neutral_df.merge(order_df, on=["pcp_index", "site"], how="left")
            .sort_values("_order")
            .drop(columns="_order")
            .reset_index(drop=True)
        )

        return aa_to_any_neutral_df

    def get_codon_to_any_probabilities_df(self):
        """Create DataFrame with non-synonymous codon mutation probabilities.

        Aggregates all non-synonymous mutations for each codon.
        Returns DataFrame with same structure as MutabilityContainer.codon_to_any_neutral_df.
        """
        # Get codon probabilities first
        codon_df = self.get_codon_probabilities_df()

        # Remove synonymous mutations
        temp_codon_neutral_df = codon_df[
            codon_df.current_aa != codon_df.transition_aa
        ].copy()

        # Capture original pcp_index and site order
        order_df = codon_df[["pcp_index", "site"]].drop_duplicates()
        order_df["_order"] = range(len(order_df))

        # Aggregate to any non-synonymous mutation
        codon_to_any_neutral_df = (
            temp_codon_neutral_df.groupby(
                [
                    "pcp_index",
                    "site",
                    "current_codon",
                    "current_aa",
                    "parent_codon",
                    "parent_aa",
                    "child_codon",
                    "child_aa",
                ]
            )
            .agg({"substitution_probability": "sum"})
            .reset_index()
        )

        # Sort by original pcp_index and site order
        codon_to_any_neutral_df = (
            codon_to_any_neutral_df.merge(
                order_df, on=["pcp_index", "site"], how="left"
            )
            .sort_values("_order")
            .drop(columns="_order")
            .reset_index(drop=True)
        )

        return codon_to_any_neutral_df


class CachedNeutralMutabilityDataset:
    """Cached version of NeutralMutabilityDataset that saves/loads dataframes to/from disk.

    Usage:
    - If cache files exist, loads them instantly
    - If cache files don't exist, creates them using NeutralMutabilityDataset and saves to cache
    - Provides the same interface as NeutralMutabilityDataset
    """

    def __init__(
        self,
        dataset_nickname,
        branch_length_mode="constant",
        branch_length=0.1,
        branch_length_scale_factor=1.0,
        numbering_scheme="imgt",
        multihit_model=None,
        subset_size=None,
        cache_dir=localify(f"ANALYSIS_CACHE/"),
        skip_nucleotide=False,
    ):
        self.dataset_nickname = dataset_nickname
        self.branch_length_mode = branch_length_mode
        self.branch_length = branch_length
        self.branch_length_scale_factor = branch_length_scale_factor
        self.numbering_scheme = numbering_scheme
        self.multihit_model = multihit_model
        self.subset_size = subset_size
        self.cache_dir = cache_dir
        self.skip_nucleotide = skip_nucleotide

        # Create cache directory if it doesn't exist
        os.makedirs(cache_dir, exist_ok=True)

        # Generate cache file paths
        # Include branch_length and scale_factor in cache name for constant mode
        if branch_length_mode == "constant":
            mode_spec = (
                f"{branch_length_mode}_bl{branch_length}_sf{branch_length_scale_factor}"
            )
        else:
            mode_spec = f"{branch_length_mode}_sf{branch_length_scale_factor}"

        # Add multihit model to cache name if used
        if multihit_model is not None:
            multihit_spec = f"_multihit_{type(multihit_model).__name__}"
        else:
            multihit_spec = ""

        if subset_size:
            cache_base = f"{cache_dir}/neutral_mutability_{dataset_nickname}_{numbering_scheme}_{mode_spec}{multihit_spec}_subset{subset_size}"
        else:
            cache_base = f"{cache_dir}/neutral_mutability_{dataset_nickname}_{numbering_scheme}_{mode_spec}{multihit_spec}"

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
        """Check if all required cache files exist."""
        return all(
            os.path.exists(filepath)
            for key, filepath in self.cache_files.items()
            if not (self.skip_nucleotide and key == "nucleotide")
        )

    def _load_from_cache(self):
        """Load dataframes from gzip-compressed cache files."""
        subset_msg = f" (subset {self.subset_size})" if self.subset_size else ""
        print(f"Loading NeutralMutabilityDataset data from gzip cache{subset_msg}...")

        try:
            # Explicitly specify site column as string to handle alphanumeric sites
            if self.skip_nucleotide:
                self.nuc_neutral_df = None
            else:
                self.nuc_neutral_df = pd.read_csv(
                    self.cache_files["nucleotide"],
                    compression="gzip",
                    dtype={"site": str},
                )
            self.aa_neutral_df = pd.read_csv(
                self.cache_files["amino_acid"], compression="gzip", dtype={"site": str}
            )
            self.aa_to_any_neutral_df = pd.read_csv(
                self.cache_files["amino_acid_to_any"],
                compression="gzip",
                dtype={"site": str},
            )
            self.codon_neutral_df = pd.read_csv(
                self.cache_files["codon"], compression="gzip", dtype={"site": str}
            )
            self.codon_to_any_neutral_df = pd.read_csv(
                self.cache_files["codon_to_any"],
                compression="gzip",
                dtype={"site": str},
            )
            self.pcp_df = pd.read_csv(self.cache_files["pcp_df"], compression="gzip")

            # Validate that DataFrames are not empty
            dfs_to_check = [
                self.aa_neutral_df,
                self.aa_to_any_neutral_df,
                self.codon_neutral_df,
                self.codon_to_any_neutral_df,
                self.pcp_df,
            ]
            if not self.skip_nucleotide:
                dfs_to_check.insert(0, self.nuc_neutral_df)
            if any(len(df) == 0 for df in dfs_to_check):
                raise ValueError("One or more cached DataFrames are empty")

            print(f"✓ Loaded from gzip cache:")
            if self.skip_nucleotide:
                print(f"  - Nucleotide DataFrame: skipped")
            else:
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

        except (pd.errors.EmptyDataError, ValueError) as e:
            # Expected errors from invalid/corrupt cache - safe to recreate
            print(f"⚠ Cache validation failed: {e}")
            print("Recreating data from source...")
            self._create_and_cache()

        except (OSError, IOError, PermissionError) as e:
            # File system errors - user needs to know!
            raise RuntimeError(
                f"Cannot access cache files in {self.cache_dir}. "
                f"Check file permissions and disk space. Error: {e}"
            ) from e

        except MemoryError as e:
            # Out of memory - definitely don't try to recreate!
            raise RuntimeError(
                f"Insufficient memory to load cache files. "
                f"Consider using subset_size parameter to reduce memory usage. "
                f"Error: {e}"
            ) from e

        except Exception as e:
            # Unexpected errors should fail loudly with actionable message
            raise RuntimeError(
                f"Unexpected error loading cache from {self.cache_dir}. "
                f"Cache may be corrupted or incompatible. "
                f"To fix: delete cache files and retry. "
                f"Cache files: {list(self.cache_files.values())}. "
                f"Error: {type(e).__name__}: {e}"
            ) from e

    def _create_and_cache(self):
        """Create new NeutralMutabilityDataset and save to cache."""
        subset_msg = f" (subset {self.subset_size})" if self.subset_size else ""
        print(f"Creating new NeutralMutabilityDataset data{subset_msg}...")

        # Create the dataset
        dataset = NeutralMutabilityDataset.of_pcp_df_with_numbering(
            dataset_nickname=self.dataset_nickname,
            branch_length_mode=self.branch_length_mode,
            branch_length=self.branch_length,
            branch_length_scale_factor=self.branch_length_scale_factor,
            numbering_scheme=self.numbering_scheme,
            multihit_model=self.multihit_model,
            subset_size=self.subset_size,
        )

        # Generate the dataframes
        if self.skip_nucleotide:
            self.nuc_neutral_df = None
        else:
            self.nuc_neutral_df = dataset.get_nucleotide_probabilities_df()
        self.aa_neutral_df = dataset.get_aa_probabilities_df()
        self.aa_to_any_neutral_df = dataset.get_aa_to_any_probabilities_df()
        self.codon_neutral_df = dataset.get_codon_probabilities_df()
        self.codon_to_any_neutral_df = dataset.get_codon_to_any_probabilities_df()
        self.pcp_df = dataset.pcp_df.copy()

        # Save to cache
        self._save_to_cache()

    def _save_to_cache(self):
        """Save dataframes to gzip-compressed cache files."""
        print("Saving data to gzip cache...")

        try:
            if self.nuc_neutral_df is not None:
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

        except (OSError, IOError, PermissionError) as e:
            # File system errors - user needs to fix permissions/disk space
            raise RuntimeError(
                f"Cannot write cache files to {self.cache_dir}. "
                f"Check directory permissions, disk space, and write access. "
                f"Error: {e}"
            ) from e

        except MemoryError as e:
            # Out of memory during serialization
            raise RuntimeError(
                f"Insufficient memory to save cache files. "
                f"DataFrames may be too large. Consider using subset_size parameter. "
                f"Error: {e}"
            ) from e

        except Exception as e:
            # Unexpected errors should fail loudly
            raise RuntimeError(
                f"Unexpected error saving cache to {self.cache_dir}. "
                f"Attempted to save: {list(self.cache_files.keys())}. "
                f"Error: {type(e).__name__}: {e}"
            ) from e
