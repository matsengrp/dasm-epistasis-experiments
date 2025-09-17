#!/usr/bin/env python3
"""
Solvent Accessibility Analysis Script

This script calculates SASA and RSA for PDB structures with V/D/J gene metadata integration.
Based on the Jupyter notebook solvent_accessibility_analysis.ipynb.

Usage:
    python run_sasa_analysis.py --help
    python run_sasa_analysis.py --scheme opig-imgt --max-files 50 --output results.csv
"""

import os
import sys
import random
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

# BioPython imports
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.DSSP import DSSP

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)


class AntibodyChainSelector(Select):
    """Select specific chains from PDB structure."""
    
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids if isinstance(chain_ids, list) else [chain_ids]
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


class ChainSelector(Select):
    """Select only specific chains from PDB structure."""

    def __init__(self, keep_chains):
        self.keep_chains = keep_chains if isinstance(keep_chains, list) else [keep_chains]

    def accept_chain(self, chain):
        return chain.id in self.keep_chains


class AntigenRemover(Select):
    """Remove antigen chains, keeping only antibody chains."""

    def __init__(self, antibody_chains, antigen_chains=None):
        self.antibody_chains = antibody_chains if isinstance(antibody_chains, list) else [antibody_chains]
        self.antigen_chains = antigen_chains if isinstance(antigen_chains, list) else ([antigen_chains] if antigen_chains else [])

    def accept_chain(self, chain):
        # If we have explicit antigen chain info, exclude those
        if self.antigen_chains and chain.id in self.antigen_chains:
            return False
        # Otherwise, keep only antibody chains (original behavior as fallback)
        return chain.id in self.antibody_chains


class StructureFilter:
    """Unified filtering pipeline for PDB structures."""

    def __init__(self, metadata_extractor, verbose=False):
        self.metadata = metadata_extractor
        self.verbose = verbose
        self.filter_stats = {}

    def apply_filters(self, pdb_directory, organism_filter=None, max_files=None):
        """Apply all filters in sequence with consolidated reporting."""

        # Stage 1: Get all available PDB files
        all_files = list(Path(pdb_directory).glob("*.pdb"))
        self._record_stage("available_files", len(all_files), all_files)

        if not self.metadata:
            # No metadata = no filtering possible
            final_files = all_files[:max_files] if max_files else all_files
            self._record_stage("final_files", len(final_files), final_files, note="No metadata filtering")
            self._print_filter_summary()
            return final_files

        # Stage 2: Metadata-based filtering
        valid_ids = self._apply_metadata_filters(organism_filter)
        files_with_metadata = self._filter_by_pdb_ids(all_files, valid_ids)
        self._record_stage("metadata_filtered", len(files_with_metadata), files_with_metadata)

        # Stage 3: Chain validation (require heavy + light + antigen)
        files_with_complete_chains = self._filter_by_complete_chains(files_with_metadata)
        self._record_stage("complete_chains", len(files_with_complete_chains), files_with_complete_chains)

        # Stage 4: Apply max_files limit
        if max_files:
            final_files = files_with_complete_chains[:max_files]
            self._record_stage("max_files_limited", len(final_files), final_files)
        else:
            final_files = files_with_complete_chains

        self._print_filter_summary()
        return final_files

    def _apply_metadata_filters(self, organism_filter):
        """Apply organism and V/J consistency filters."""
        if organism_filter:
            # Get organism-filtered PDBs
            organism_ids = self._get_organism_filtered_ids(organism_filter)
            # Get V/J consistent PDBs
            consistent_ids, kept, removed = self.metadata.get_pdbs_with_consistent_vj_genes()
            # Intersect both filters
            valid_ids = set(organism_ids) & set(consistent_ids)

            if self.verbose:
                print(f"Organism filter: {len(organism_ids)} PDBs match '{organism_filter}'")
                print(f"V/J consistency: {kept} PDBs kept, {removed} removed")
                print(f"Combined filters: {len(valid_ids)} PDBs")
        else:
            # Only V/J consistency filter
            valid_ids, kept, removed = self.metadata.get_pdbs_with_consistent_vj_genes()
            if self.verbose:
                print(f"V/J consistency filter: {kept} PDBs kept, {removed} removed")

        return valid_ids

    def _get_organism_filtered_ids(self, organism_filter):
        """Get PDB IDs matching organism filter."""
        if not self.metadata.abid_df.empty:
            filtered_df = self.metadata.abid_df[
                self.metadata.abid_df['organism'].str.contains(organism_filter, case=False, na=False)
            ]
            return set(filtered_df['pdbid'].unique())
        return set()

    def _filter_by_pdb_ids(self, pdb_files, valid_ids):
        """Keep only files with PDB IDs in the valid set."""
        valid_ids_upper = {str(pid).upper() for pid in valid_ids}
        return [f for f in pdb_files if f.stem.upper() in valid_ids_upper]

    def _filter_by_complete_chains(self, pdb_files):
        """Keep only files with heavy, light, and antigen chain information."""
        valid_files = []
        rejected_count = {"no_heavy": 0, "no_light": 0, "no_antigen": 0}

        for pdb_file in pdb_files:
            pdb_id = pdb_file.stem
            chain_info = self.metadata.get_chain_info(pdb_id)

            has_heavy = bool(chain_info.get('heavy_chain'))
            has_light = bool(chain_info.get('light_chain'))
            has_antigen = bool(chain_info.get('antigen_chains'))

            if has_heavy and has_light and has_antigen:
                valid_files.append(pdb_file)
            else:
                # Track rejection reasons
                if not has_heavy:
                    rejected_count["no_heavy"] += 1
                if not has_light:
                    rejected_count["no_light"] += 1
                if not has_antigen:
                    rejected_count["no_antigen"] += 1

                if self.verbose:
                    missing = []
                    if not has_heavy: missing.append('heavy')
                    if not has_light: missing.append('light')
                    if not has_antigen: missing.append('antigen')
                    print(f"Skipping {pdb_id}: missing {', '.join(missing)} chain(s)")

        if self.verbose and any(rejected_count.values()):
            total_rejected = len(pdb_files) - len(valid_files)
            print(f"Chain filtering: {len(valid_files)} kept, {total_rejected} rejected")
            for reason, count in rejected_count.items():
                if count > 0:
                    print(f"  {reason}: {count} files")

        return valid_files

    def _record_stage(self, stage_name, count, files, note=None):
        """Record statistics for a filtering stage."""
        self.filter_stats[stage_name] = {
            'count': count,
            'files': files,
            'note': note
        }

    def _print_filter_summary(self):
        """Print consolidated filtering summary."""
        if not self.verbose:
            return

        print("\n=== Filtering Pipeline Summary ===")
        for stage, data in self.filter_stats.items():
            count = data['count']
            note = f" ({data['note']})" if data.get('note') else ""
            print(f"{stage.replace('_', ' ').title()}: {count} files{note}")

        # Show sample of final files
        if self.filter_stats:
            final_stage = list(self.filter_stats.keys())[-1]
            final_files = self.filter_stats[final_stage]['files'][:5]
            if final_files:
                print(f"Sample final files: {[f.stem for f in final_files]}")


class ProteinMetadataExtractor:
    """Extract protein metadata from SAbDab summary tables."""

    def __init__(self, verbose=False):
        self.abid_info_path = "/fh/fast/matsen_e/shared/sabdab_pb/sabdab_summary_2024-01-26_abid_info.tsv"
        self.chain_info_path = "/fh/fast/matsen_e/shared/sabdab_pb/sabdab_summary_all_2024-01-26.tsv"
        self.abid_df = None
        self.chain_df = None
        self.verbose = verbose
        self.pdb_lookup = {}
        self.load_sabdab_data()
    
    def load_sabdab_data(self):
        """Load the SAbDab summary tables with V/D/J gene and chain information."""
        # Load abid_info table (V/D/J genes, organism)
        if os.path.exists(self.abid_info_path):
            self.abid_df = pd.read_table(self.abid_info_path)
            # Select relevant columns: pdbid, organism, va, ja, vb, jb
            required_cols = ["pdbid", "organism", "va", "ja", "vb", "jb"]
            available_cols = [col for col in required_cols if col in self.abid_df.columns]
            self.abid_df = self.abid_df[available_cols]
            if self.verbose:
                print(f"Loaded abid info with {len(self.abid_df)} entries")
        else:
            if self.verbose:
                print(f"Warning: Abid info file not found at {self.abid_info_path}")
            self.abid_df = pd.DataFrame()

        # Load chain_info table (chain identifications)
        if os.path.exists(self.chain_info_path):
            self.chain_df = pd.read_table(self.chain_info_path)
            # Select relevant columns: pdb, Hchain, Lchain, antigen_chain
            required_cols = ["pdb", "Hchain", "Lchain", "antigen_chain"]
            available_cols = [col for col in required_cols if col in self.chain_df.columns]
            self.chain_df = self.chain_df[available_cols]
            if self.verbose:
                print(f"Loaded chain info with {len(self.chain_df)} entries")
        else:
            if self.verbose:
                print(f"Warning: Chain info file not found at {self.chain_info_path}")
            self.chain_df = pd.DataFrame()

        # Create lookup dictionary for faster access
        self.pdb_lookup = {}

        # Process abid_df
        if not self.abid_df.empty:
            for _, row in self.abid_df.iterrows():
                pdb_id = str(row['pdbid']).strip()
                if pdb_id not in self.pdb_lookup:
                    self.pdb_lookup[pdb_id] = {}
                self.pdb_lookup[pdb_id].update(row.to_dict())
                # Also add case variations
                self.pdb_lookup[pdb_id.upper()] = self.pdb_lookup[pdb_id]
                self.pdb_lookup[pdb_id.lower()] = self.pdb_lookup[pdb_id]

        # Process chain_df
        if not self.chain_df.empty:
            for _, row in self.chain_df.iterrows():
                pdb_id = str(row['pdb']).strip()
                if pdb_id not in self.pdb_lookup:
                    self.pdb_lookup[pdb_id] = {}
                self.pdb_lookup[pdb_id].update(row.to_dict())
                # Also add case variations
                self.pdb_lookup[pdb_id.upper()] = self.pdb_lookup[pdb_id]
                self.pdb_lookup[pdb_id.lower()] = self.pdb_lookup[pdb_id]
    
    def get_protein_metadata(self, pdb_id):
        """Get protein metadata for a given PDB ID."""
        if not self.pdb_lookup:
            return {'pdb_id': pdb_id}

        # Try to find the entry in our lookup dictionary
        pdb_id_clean = str(pdb_id).strip()
        entry = None

        # Try different case variations
        for key in [pdb_id_clean.upper(), pdb_id_clean.lower(), pdb_id_clean]:
            if key in self.pdb_lookup:
                entry = self.pdb_lookup[key]
                break

        if entry is None:
            return {'pdb_id': pdb_id}

        metadata = {
            'pdb_id': pdb_id,
            'organism': entry.get('organism', ''),
            'v_gene_light': entry.get('va', ''),   # Light chain V gene
            'j_gene_light': entry.get('ja', ''),   # Light chain J gene
            'v_gene_heavy': entry.get('vb', ''),   # Heavy chain V gene
            'j_gene_heavy': entry.get('jb', ''),   # Heavy chain J gene
            'heavy_chain_id': entry.get('Hchain', ''),
            'light_chain_id': entry.get('Lchain', ''),
            'antigen_chains': entry.get('antigen_chain', '')
        }

        return metadata
    
    def get_pdbs_with_consistent_vj_genes(self):
        """Get PDB IDs where all entries have the same V and J genes."""
        if self.abid_df is None or len(self.abid_df) == 0:
            return set(), 0, 0

        # Group by PDB ID and check for consistency
        consistent_pdbs = set()
        total_pdbs = self.abid_df['pdbid'].nunique()

        for pdb_id, group in self.abid_df.groupby('pdbid'):
            # Check if all V and J genes are the same within this PDB
            va_unique = group['va'].dropna().unique() if 'va' in group.columns else []
            ja_unique = group['ja'].dropna().unique() if 'ja' in group.columns else []
            vb_unique = group['vb'].dropna().unique() if 'vb' in group.columns else []
            jb_unique = group['jb'].dropna().unique() if 'jb' in group.columns else []

            # PDB is consistent if each gene type has at most one unique value
            is_consistent = (
                len(va_unique) <= 1 and
                len(ja_unique) <= 1 and
                len(vb_unique) <= 1 and
                len(jb_unique) <= 1
            )

            if is_consistent:
                consistent_pdbs.add(pdb_id)

        removed_count = total_pdbs - len(consistent_pdbs)

        if self.verbose:
            print(f"V/J gene consistency filter: {len(consistent_pdbs)} PDBs kept, {removed_count} PDBs removed (out of {total_pdbs} total)")

        return consistent_pdbs, len(consistent_pdbs), removed_count


    def get_chain_info(self, pdb_id):
        """Get chain identification information for a PDB ID."""
        metadata = self.get_protein_metadata(pdb_id)

        # Parse antigen chains (can be multiple, separated by ' | ')
        antigen_chain_str = metadata.get('antigen_chains', '')
        antigen_chains = []
        if antigen_chain_str and pd.notna(antigen_chain_str):
            antigen_chains = [chain.strip() for chain in str(antigen_chain_str).split('|')]

        # Get antibody chains
        antibody_chains = []
        if metadata.get('heavy_chain_id'):
            antibody_chains.append(metadata['heavy_chain_id'])
        if metadata.get('light_chain_id'):
            antibody_chains.append(metadata['light_chain_id'])

        return {
            'antibody_chains': antibody_chains,
            'antigen_chains': antigen_chains,
            'heavy_chain': metadata.get('heavy_chain_id', ''),
            'light_chain': metadata.get('light_chain_id', '')
        }


class SASACalculator:
    """Calculate SASA and RSA for PDB structures with protein metadata integration."""
    
    def __init__(self, temp_dir="_temp_sasa", metadata_extractor=None):
        import os
        # Make temp directory unique per process to avoid conflicts
        unique_temp_dir = f"{temp_dir}_{os.getpid()}"
        self.temp_dir = Path(unique_temp_dir)
        self.temp_dir.mkdir(exist_ok=True)
        self.metadata_extractor = metadata_extractor
        
        # Wilke reference values for calculating absolute ASA from relative ASA
        self.wilke_reference = {
            'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0,
            'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0,
            'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0,
            'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0,
            'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0
        }

    def _create_nan_result(self, pdb_id, scenario, metadata):
        """Create a DataFrame with NaN values when SASA calculation cannot be performed."""
        return pd.DataFrame([{
            'chain_id': np.nan,
            'residue_number': np.nan,
            'insertion_code': np.nan,
            'amino_acid': np.nan,
            'sasa': np.nan,
            'rsa': np.nan,
            'ca_coordinates': np.nan,
            'scenario': scenario,
            # Add protein metadata
            'pdb_id': pdb_id,
            'organism': metadata.get('organism', ''),
            'v_gene_light': metadata.get('v_gene_light', ''),
            'j_gene_light': metadata.get('j_gene_light', ''),
            'v_gene_heavy': metadata.get('v_gene_heavy', ''),
            'j_gene_heavy': metadata.get('j_gene_heavy', ''),
            'heavy_chain_id': metadata.get('heavy_chain_id', ''),
            'light_chain_id': metadata.get('light_chain_id', ''),
            'antigen_chains': metadata.get('antigen_chains', '')
        }])
    
    def _calculate_sasa(self, pdb_path, scenario="heavy_light_antigen", site_range=None):
        """Internal method to calculate SASA with protein metadata.

        Parameters:
        -----------
        scenario : str
            One of: "heavy_light_antigen", "heavy_light", "heavy_only"
        """
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure('pdb', pdb_path)

        # Extract PDB ID and get metadata
        pdb_id = Path(pdb_path).stem
        metadata = {}
        chain_info = {}
        if self.metadata_extractor:
            metadata = self.metadata_extractor.get_protein_metadata(pdb_id)
            chain_info = self.metadata_extractor.get_chain_info(pdb_id)

            # Always use automatic chain detection from metadata

        # Determine which chains to keep based on scenario
        heavy_chain = chain_info.get('heavy_chain', '')
        light_chain = chain_info.get('light_chain', '')
        antigen_chains = chain_info.get('antigen_chains', [])

        # Determine which chains to keep based on scenario
        if scenario == "heavy_light_antigen":
            keep_chains = [heavy_chain, light_chain] + antigen_chains
            target_chains = [heavy_chain, light_chain]  # Only analyze antibody residues
        elif scenario == "heavy_light":
            keep_chains = [heavy_chain, light_chain]
            target_chains = [heavy_chain, light_chain]
        elif scenario == "heavy_only":
            keep_chains = [heavy_chain]
            target_chains = [heavy_chain]
        else:
            raise ValueError(f"Invalid scenario: {scenario}")

        # Remove empty chain IDs
        keep_chains = [c for c in keep_chains if c]
        target_chains = [c for c in target_chains if c]

        if not target_chains:
            print(f"Warning: No valid target chains for {pdb_id} in scenario {scenario}")
            return self._create_nan_result(pdb_id, scenario, metadata)
        
        # Create temporary PDB file
        temp_id = random.randint(10000, 99999)
        temp_pdb_path = self.temp_dir / f"temp_{temp_id}.pdb"
        
        try:
            io = PDBIO()
            io.set_structure(structure)

            # Use strict chain filtering - keep only the chains specified for this scenario
            io.save(str(temp_pdb_path), ChainSelector(keep_chains))
            
            # Reload the filtered structure
            filtered_structure = parser.get_structure('filtered', str(temp_pdb_path))
            filtered_model = filtered_structure[0]
            
            # Calculate DSSP with warning suppression
            import warnings
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="parse error at line 1")
                dssp = DSSP(filtered_model, str(temp_pdb_path), dssp="mkdssp", acc_array="Wilke")
            
            # Extract data for antibody chains only
            sasa_data = []
            
            for key in dssp.keys():
                chain_id, res_id = key

                # Only process target chains (antibody chains for analysis)
                if chain_id not in target_chains:
                    continue
                
                # Extract residue information
                residue_num = res_id[1]
                insertion_code = res_id[2].strip()
                
                # Apply site range filter if specified
                if site_range and (residue_num < site_range[0] or residue_num > site_range[1]):
                    continue
                
                # Get DSSP data
                dssp_data = dssp[key]
                aa = dssp_data[1]  # Amino acid
                rel_asa = dssp_data[3]  # Relative accessible surface area (Wilke normalized, 0-1)

                # Calculate absolute ASA from relative ASA
                max_asa = self.wilke_reference[aa]  # Will raise KeyError for unknown amino acids
                abs_asa = rel_asa * max_asa
                
                # Get CA coordinates if available
                ca_coord = None
                try:
                    residue = filtered_model[chain_id][res_id]
                    if "CA" in residue:
                        ca_coord = tuple(residue["CA"].coord)
                except:
                    pass
                
                # Build result row with both SASA data and metadata
                row_data = {
                    'chain_id': chain_id,
                    'residue_number': residue_num,
                    'insertion_code': insertion_code,
                    'amino_acid': aa,
                    'sasa': abs_asa,  # Absolute ASA in Å²
                    'rsa': rel_asa,   # Relative ASA (0-1 scale)
                    'ca_coordinates': ca_coord,
                    'scenario': scenario,
                    # Add protein metadata
                    'pdb_id': pdb_id,
                    'organism': metadata.get('organism', ''),
                    'v_gene_light': metadata.get('v_gene_light', ''),
                    'j_gene_light': metadata.get('j_gene_light', ''),
                    'v_gene_heavy': metadata.get('v_gene_heavy', ''),
                    'j_gene_heavy': metadata.get('j_gene_heavy', ''),
                    'heavy_chain_id': metadata.get('heavy_chain_id', ''),
                    'light_chain_id': metadata.get('light_chain_id', ''),
                    'antigen_chains': metadata.get('antigen_chains', '')
                }
                
                sasa_data.append(row_data)
            
            return pd.DataFrame(sasa_data)
            
        finally:
            # Clean up temporary file
            if temp_pdb_path.exists():
                temp_pdb_path.unlink()
    
    def calculate_multi_scenario_sasa(self, pdb_path, site_range=None):
        """Calculate SASA for all three scenarios and compare."""

        scenarios = ["heavy_light_antigen", "heavy_light", "heavy_only"]
        scenario_results = {}

        # Calculate SASA for each scenario
        for scenario in scenarios:
            try:
                results = self._calculate_sasa(pdb_path, scenario=scenario,
                                             site_range=site_range)
                if not results.empty:
                    scenario_results[scenario] = results
            except Exception as e:
                # Continue with other scenarios if one fails
                print(f"Warning: Failed to calculate SASA for scenario {scenario}: {e}")
                continue

        if not scenario_results:
            return pd.DataFrame()

        # Merge all scenarios
        metadata_cols = ['pdb_id', 'organism', 'v_gene_light', 'j_gene_light', 'v_gene_heavy', 'j_gene_heavy',
                        'heavy_chain_id', 'light_chain_id', 'antigen_chains']
        base_cols = ['chain_id', 'residue_number', 'insertion_code', 'amino_acid', 'ca_coordinates']

        # Start with first available scenario
        first_scenario = list(scenario_results.keys())[0]
        comparison = scenario_results[first_scenario][base_cols + metadata_cols + ['sasa', 'rsa']].copy()
        comparison = comparison.rename(columns={
            'sasa': f'sasa_{first_scenario}',
            'rsa': f'rsa_{first_scenario}'
        })

        # Merge other scenarios
        for scenario in list(scenario_results.keys())[1:]:
            scenario_data = scenario_results[scenario][['chain_id', 'residue_number', 'insertion_code', 'sasa', 'rsa']]
            scenario_data = scenario_data.rename(columns={
                'sasa': f'sasa_{scenario}',
                'rsa': f'rsa_{scenario}'
            })
            comparison = pd.merge(
                comparison, scenario_data,
                on=['chain_id', 'residue_number', 'insertion_code'],
                how='outer'
            )

        # Calculate differences and relative changes
        if 'sasa_heavy_light_antigen' in comparison.columns and 'sasa_heavy_light' in comparison.columns:
            comparison['sasa_antigen_effect'] = comparison['sasa_heavy_light_antigen'] - comparison['sasa_heavy_light']
            comparison['rsa_antigen_effect'] = comparison['rsa_heavy_light_antigen'] - comparison['rsa_heavy_light']

            comparison['sasa_antigen_relative'] = (
                comparison['sasa_antigen_effect'] / comparison['sasa_heavy_light']
            ).replace([np.inf, -np.inf], np.nan)

            comparison['rsa_antigen_relative'] = (
                comparison['rsa_antigen_effect'] / comparison['rsa_heavy_light']
            ).replace([np.inf, -np.inf], np.nan)

        if 'sasa_heavy_light' in comparison.columns and 'sasa_heavy_only' in comparison.columns:
            comparison['sasa_light_effect'] = comparison['sasa_heavy_light'] - comparison['sasa_heavy_only']
            comparison['rsa_light_effect'] = comparison['rsa_heavy_light'] - comparison['rsa_heavy_only']

            comparison['sasa_light_relative'] = (
                comparison['sasa_light_effect'] / comparison['sasa_heavy_only']
            ).replace([np.inf, -np.inf], np.nan)

            comparison['rsa_light_relative'] = (
                comparison['rsa_light_effect'] / comparison['rsa_heavy_only']
            ).replace([np.inf, -np.inf], np.nan)

        return comparison
    
    def cleanup(self):
        """Remove temporary directory."""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)


def process_pdb_directory(pdb_directory, site_range=None,
                         output_file=None, max_files=None, include_metadata=True,
                         metadata_extractor=None, organism_filter=None, verbose=False):
    """
    Process all PDB files in a directory with simplified unified filtering pipeline.

    Parameters:
    -----------
    pdb_directory : str or Path
        Directory containing PDB files
    site_range : tuple or None
        Range of residue numbers to include (min, max) or None for all
    output_file : str or None
        Output CSV file path or None to return DataFrame only
    max_files : int or None
        Maximum number of files to process (for testing)
    include_metadata : bool
        Whether to include V/D/J gene and organism metadata
    metadata_extractor : ProteinMetadataExtractor or None
        Metadata extractor instance, or None to create one
    organism_filter : str or None
        Filter to process only structures from specific organism
    verbose : bool
        Whether to print verbose output

    Returns:
    --------
    pd.DataFrame
        Combined results for all PDB files with metadata
    """
    # Initialize metadata extractor if requested
    if include_metadata and metadata_extractor is None:
        metadata_extractor = ProteinMetadataExtractor(verbose=verbose)

    # Apply unified filtering pipeline
    filter_pipeline = StructureFilter(metadata_extractor if include_metadata else None, verbose=verbose)
    pdb_files = filter_pipeline.apply_filters(pdb_directory, organism_filter, max_files)

    if not pdb_files:
        if organism_filter:
            print(f"No PDB files found matching organism '{organism_filter}' in {pdb_directory}")
        else:
            print(f"No valid PDB files found in {pdb_directory}")
        return pd.DataFrame()

    calculator = SASACalculator(metadata_extractor=metadata_extractor if include_metadata else None)
    all_results = []
    failed_files = []

    print(f"Processing {len(pdb_files)} validated PDB files...")

    for pdb_file in tqdm(pdb_files, desc="Processing PDBs"):
        try:
            results = calculator.calculate_multi_scenario_sasa(str(pdb_file), site_range)
            if not results.empty:
                all_results.append(results)
            else:
                failed_files.append(pdb_file.name)

        except Exception as e:
            if verbose:
                print(f"Error processing {pdb_file.name}: {str(e)}")
            failed_files.append(pdb_file.name)
            continue
    
    # Cleanup temporary files
    try:
        calculator.cleanup()
    except Exception as cleanup_error:
        if verbose:
            print(f"Warning: Cleanup failed: {cleanup_error}")

    if not all_results:
        print("No files were successfully processed")
        return pd.DataFrame()

    # Combine all results
    combined_results = pd.concat(all_results, ignore_index=True)

    # Save intermediate results to prevent total data loss
    if output_file and len(all_results) > 0:
        temp_output = f"{output_file}.partial"
        try:
            combined_results.to_csv(temp_output, index=False)
            if verbose:
                print(f"Intermediate results saved to {temp_output}")
        except Exception as save_error:
            if verbose:
                print(f"Warning: Could not save intermediate results: {save_error}")
    
    # Reorder columns for better readability
    if include_metadata:
        column_order = [
            'pdb_id', 'organism', 'v_gene_light', 'j_gene_light', 'v_gene_heavy', 'j_gene_heavy',
            'heavy_chain_id', 'light_chain_id', 'antigen_chains', 'chain_id',
            'residue_number', 'insertion_code', 'amino_acid',
            # Raw SASA values
            'sasa_heavy_light_antigen', 'sasa_heavy_light', 'sasa_heavy_only',
            'rsa_heavy_light_antigen', 'rsa_heavy_light', 'rsa_heavy_only',
            # Effect calculations
            'sasa_antigen_effect', 'rsa_antigen_effect', 'sasa_antigen_relative', 'rsa_antigen_relative',
            'sasa_light_effect', 'rsa_light_effect', 'sasa_light_relative', 'rsa_light_relative',
            'ca_coordinates'
        ]
    else:
        column_order = [
            'pdb_id', 'chain_id', 'residue_number', 'insertion_code', 'amino_acid',
            'sasa_heavy_light_antigen', 'sasa_heavy_light', 'sasa_heavy_only',
            'rsa_heavy_light_antigen', 'rsa_heavy_light', 'rsa_heavy_only',
            'sasa_antigen_effect', 'rsa_antigen_effect', 'sasa_antigen_relative', 'rsa_antigen_relative',
            'sasa_light_effect', 'rsa_light_effect', 'sasa_light_relative', 'rsa_light_relative',
            'ca_coordinates'
        ]
    
    # Only include columns that exist in the DataFrame
    available_columns = [col for col in column_order if col in combined_results.columns]
    combined_results = combined_results.reindex(columns=available_columns)
    
    # Print summary
    print(f"\n=== Processing Results ===")
    print(f"Successfully processed: {len(all_results)} files")
    print(f"Failed to process: {len(failed_files)} files")
    print(f"Total residues analyzed: {len(combined_results)}")
    print(f"Unique PDB structures: {combined_results['pdb_id'].nunique()}")
    
    if include_metadata:
        if 'v_gene_heavy' in combined_results.columns:
            v_heavy_count = combined_results['v_gene_heavy'].notna().sum()
            j_heavy_count = combined_results['j_gene_heavy'].notna().sum()
            print(f"Residues with heavy chain V gene data: {v_heavy_count}")
            print(f"Residues with heavy chain J gene data: {j_heavy_count}")
        if 'v_gene_light' in combined_results.columns:
            v_light_count = combined_results['v_gene_light'].notna().sum()
            j_light_count = combined_results['j_gene_light'].notna().sum()
            print(f"Residues with light chain V gene data: {v_light_count}")
            print(f"Residues with light chain J gene data: {j_light_count}")
        if 'organism' in combined_results.columns:
            print(f"Unique organisms: {combined_results['organism'].nunique()}")
    
    # Save to file if requested
    if output_file:
        combined_results.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    
    if verbose and failed_files:
        print(f"\nFailed files: {failed_files[:10]}{'...' if len(failed_files) > 10 else ''}")
    
    return combined_results


def validate_output_file(output_path):
    """Test if the output file can be created and written to."""
    try:
        # Check if parent directory exists and is writable
        output_file = Path(output_path)
        parent_dir = output_file.parent

        if not parent_dir.exists():
            print(f"Error: Output directory does not exist: {parent_dir}")
            return False

        if not os.access(parent_dir, os.W_OK):
            print(f"Error: No write permission for directory: {parent_dir}")
            return False

        # Test writing to the file
        test_content = "pdb_id,test\ntest_pdb,123\n"
        with open(output_path, 'w') as f:
            f.write(test_content)

        # Clean up test file
        os.remove(output_path)

        print(f"Output file validation successful: {output_path}")
        return True

    except Exception as e:
        print(f"Error: Cannot write to output file {output_path}: {e}")
        return False


def print_sample_results(results, num_samples=5):
    """Print sample results with metadata."""
    if results.empty:
        print("No results to display")
        return
    
    print(f"\n=== Sample Results ({num_samples} rows) ===")
    
    # Select relevant columns for display
    display_cols = ['pdb_id', 'organism', 'v_gene_heavy', 'j_gene_heavy', 'chain_id',
                   'residue_number', 'amino_acid', 'sasa_antigen_effect', 'sasa_light_effect']
    available_cols = [col for col in display_cols if col in results.columns]
    
    sample_data = results[available_cols].head(num_samples)
    print(sample_data.to_string(index=False))


def main():
    """Main function for command line execution."""
    parser = argparse.ArgumentParser(
        description="Analyze solvent accessibility (SASA/RSA) for PDB structures with V/D/J gene metadata.\n"
                   "Calculates SASA in three scenarios: heavy+light+antigen, heavy+light, and heavy-only,\n"
                   "then compares to determine antigen and light chain effects.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_sasa_analysis.py --scheme opig-imgt --output results.csv
  python run_sasa_analysis.py --scheme rcsb --max-files 10 --verbose
  python run_sasa_analysis.py --pdb-dir /path/to/pdbs --no-metadata --output basic_results.csv

Output columns include:
  - sasa_heavy_light_antigen, sasa_heavy_light, sasa_heavy_only: Raw SASA values
  - sasa_antigen_effect: SASA change due to antigen binding (heavy+light+antigen - heavy+light)
  - sasa_light_effect: SASA change due to light chain (heavy+light - heavy_only)
  - *_relative columns: Relative changes as fractions
        """
    )
    
    parser.add_argument('--scheme', default='opig-imgt', 
                       choices=['opig-imgt', 'opig-chothia', 'rcsb'],
                       help='PDB numbering scheme to use (default: opig-imgt)')
    
    parser.add_argument('--pdb-dir', type=str, 
                       help='Custom PDB directory path (overrides --scheme)')
    
    parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output CSV file path')
    
    # Antibody chains are now automatically detected from metadata
    
    parser.add_argument('--site-range', nargs=2, type=int, metavar=('MIN', 'MAX'),
                       help='Residue number range to analyze (e.g., --site-range 1 150)')
    
    parser.add_argument('--max-files', type=int,
                       help='Maximum number of PDB files to process (for testing)')
    
    parser.add_argument('--no-metadata', action='store_true',
                       help='Disable V/D/J gene metadata extraction')
    
    # Metadata file paths are now hardcoded in ProteinMetadataExtractor
    
    parser.add_argument('--organism', type=str,
                       help='Filter analysis to specific organism (e.g., "Homo sapiens", "Mus musculus")')

    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Print verbose output')
    
    args = parser.parse_args()
    
    # Set up environment
    os.environ.setdefault('LIBCIFPP_DATA_DIR', '/home/nharel/miniforge3/envs/netam_env/share/libcifpp')

    # Validate output file early to catch issues before processing
    if not validate_output_file(args.output):
        print("Fix output file issues before running analysis.")
        sys.exit(1)

    # Determine PDB directory
    if args.pdb_dir:
        pdb_directory = args.pdb_dir
    else:
        base_dir = "/fh/fast/matsen_e/shared/bcr-mut-sel/sabdab/pdb-db"
        pdb_directory = f"{base_dir}/pdb/{args.scheme}"
    
    if not os.path.exists(pdb_directory):
        print(f"Error: PDB directory not found: {pdb_directory}")
        sys.exit(1)
    
    # Convert site range to tuple if provided
    site_range = tuple(args.site_range) if args.site_range else None
    
    # Set up metadata extraction
    include_metadata = not args.no_metadata
    metadata_extractor = None
    if include_metadata:
        metadata_extractor = ProteinMetadataExtractor(verbose=args.verbose)
    
    # Print configuration
    print(f"=== Multi-Scenario SASA Analysis Configuration ===")
    print(f"PDB directory: {pdb_directory}")
    print(f"Filtering: Unified pipeline (V/J consistency + complete chains)")
    print(f"Chain requirement: heavy + light + antigen chains required")
    print(f"Site range: {site_range if site_range else 'All residues'}")
    print(f"Max files: {args.max_files if args.max_files else 'All files'}")
    print(f"Include metadata: {include_metadata}")
    print(f"Organism filter: {args.organism if args.organism else 'All organisms'}")
    print(f"Output file: {args.output}")
    print(f"Scenarios: heavy+light+antigen, heavy+light, heavy-only")
    print()
    
    # Run analysis
    try:
        results = process_pdb_directory(
            pdb_directory=pdb_directory,
            site_range=site_range,
            output_file=args.output,
            max_files=args.max_files,
            include_metadata=include_metadata,
            metadata_extractor=metadata_extractor,
            organism_filter=args.organism,
            verbose=args.verbose
        )
        
        if not results.empty:
            print_sample_results(results, num_samples=10)
            
            # Basic statistics
            print(f"\n=== Basic Statistics ===")
            if 'sasa_antigen_effect' in results.columns:
                print(f"Mean SASA antigen effect: {results['sasa_antigen_effect'].mean():.3f} ± {results['sasa_antigen_effect'].std():.3f} A^2")
                print(f"Mean RSA antigen effect: {results['rsa_antigen_effect'].mean():.4f} ± {results['rsa_antigen_effect'].std():.4f}")
            if 'sasa_light_effect' in results.columns:
                print(f"Mean SASA light chain effect: {results['sasa_light_effect'].mean():.3f} ± {results['sasa_light_effect'].std():.3f} A^2")
                print(f"Mean RSA light chain effect: {results['rsa_light_effect'].mean():.4f} ± {results['rsa_light_effect'].std():.4f}")
            
            print(f"\nAnalysis complete! Results saved to {args.output}")
        else:
            print("No results generated!")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()