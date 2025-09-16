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
from Bio import PDB
from Bio.SeqUtils import seq1

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)


class AntibodyChainSelector(Select):
    """Select specific chains from PDB structure."""
    
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids if isinstance(chain_ids, list) else [chain_ids]
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


class AntigenRemover(Select):
    """Remove antigen chains, keeping only antibody chains."""
    
    def __init__(self, antibody_chains):
        self.antibody_chains = antibody_chains if isinstance(antibody_chains, list) else [antibody_chains]
    
    def accept_chain(self, chain):
        return chain.id in self.antibody_chains


class ProteinMetadataExtractor:
    """Extract protein metadata from SAbDab summary table."""
    
    def __init__(self, sabdab_summary_path="/fh/fast/matsen_e/shared/bcr-mut-sel/sabdab/sabdab_summary_for_dnsm.tsv", verbose=False):
        self.sabdab_summary_path = sabdab_summary_path
        self.sabdab_df = None
        self.verbose = verbose
        self.pdb_lookup = {}
        self.load_sabdab_summary()
    
    def load_sabdab_summary(self):
        """Load the SAbDab summary table with V/D/J gene information."""
        if not os.path.exists(self.sabdab_summary_path):
            if self.verbose:
                print(f"Warning: SAbDab summary not found at {self.sabdab_summary_path}")
            self.sabdab_df = pd.DataFrame()
            return
        
        self.sabdab_df = pd.read_table(self.sabdab_summary_path)
        # Select relevant columns as done in dms_viz.py
        self.sabdab_df = self.sabdab_df[
            ["organism", "pdbid", "abid", "vb", "jb", "chainseq_b"]
        ]
        
        # Create lookup dictionary for faster access
        self.pdb_lookup = {}
        for _, row in self.sabdab_df.iterrows():
            pdb_id = str(row['pdbid']).strip()
            self.pdb_lookup[pdb_id.upper()] = row
            self.pdb_lookup[pdb_id.lower()] = row
        
        if self.verbose:
            print(f"Loaded SAbDab summary with {len(self.sabdab_df)} entries")
        
        return self.sabdab_df
    
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
            'abid': entry.get('abid', ''),
            'v_gene': entry.get('vb', ''),  # V gene family (heavy chain)
            'j_gene': entry.get('jb', ''),  # J gene family (heavy chain)
            'chain_sequence': entry.get('chainseq_b', '')
        }
        
        # Extract chain information from abid if available
        if pd.notna(entry.get('abid', '')):
            abid = str(entry.get('abid', ''))
            if len(abid) >= 2:
                metadata['heavy_chain_id'] = abid[-2]  # Second to last character
                metadata['light_chain_id'] = abid[-1]  # Last character
        
        return metadata
    
    def get_available_pdb_ids(self):
        """Get list of PDB IDs available in the summary table."""
        if self.sabdab_df is None or len(self.sabdab_df) == 0:
            return []
        return sorted(self.sabdab_df['pdbid'].unique())


class SASACalculator:
    """Calculate SASA and RSA for PDB structures with protein metadata integration."""
    
    def __init__(self, temp_dir="_temp_sasa", metadata_extractor=None):
        self.temp_dir = Path(temp_dir)
        self.temp_dir.mkdir(exist_ok=True)
        self.metadata_extractor = metadata_extractor
        
        # Standard amino acid volumes for RSA calculation (Angstrom squared)
        self.aa_max_asa = {
            'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0,
            'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0,
            'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0,
            'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0,
            'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0
        }
    
    def _calculate_sasa(self, pdb_path, antibody_chains, with_antigen=True, site_range=None):
        """Internal method to calculate SASA with protein metadata."""
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure('pdb', pdb_path)
        model = structure[0]
        
        # Extract PDB ID and get metadata
        pdb_id = Path(pdb_path).stem
        metadata = {}
        if self.metadata_extractor:
            metadata = self.metadata_extractor.get_protein_metadata(pdb_id)
        
        # Create temporary PDB file
        temp_id = random.randint(10000, 99999)
        temp_pdb_path = self.temp_dir / f"temp_{temp_id}.pdb"
        
        try:
            io = PDBIO()
            io.set_structure(structure)
            
            if with_antigen:
                # Keep all chains
                io.save(str(temp_pdb_path))
            else:
                # Keep only antibody chains
                io.save(str(temp_pdb_path), AntigenRemover(antibody_chains))
            
            # Reload the filtered structure
            filtered_structure = parser.get_structure('filtered', str(temp_pdb_path))
            filtered_model = filtered_structure[0]
            
            # Calculate DSSP
            dssp = DSSP(filtered_model, str(temp_pdb_path), dssp="mkdssp", acc_array="Wilke", file_type="PDB")
            
            # Extract data for antibody chains only
            sasa_data = []
            
            for key in dssp.keys():
                chain_id, res_id = key
                
                # Only process antibody chains
                if chain_id not in antibody_chains:
                    continue
                
                # Extract residue information
                residue_num = res_id[1]
                insertion_code = res_id[2].strip()
                
                # Apply site range filter if specified
                if site_range and (residue_num < site_range[0] or residue_num > site_range[1]):
                    continue
                
                # Get DSSP data
                aa = dssp[key][1]  # Amino acid
                sasa = dssp[key][3]  # Accessible surface area
                
                # Calculate RSA
                max_asa = self.aa_max_asa.get(aa, 200.0)  # Default if unknown AA
                rsa = sasa / max_asa if max_asa > 0 else 0.0
                
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
                    'sasa': sasa,
                    'rsa': rsa,
                    'ca_coordinates': ca_coord,
                    'with_antigen': with_antigen,
                    # Add protein metadata
                    'pdb_id': pdb_id,
                    'organism': metadata.get('organism', ''),
                    'abid': metadata.get('abid', ''),
                    'v_gene': metadata.get('v_gene', ''),
                    'j_gene': metadata.get('j_gene', ''),
                    'heavy_chain_id': metadata.get('heavy_chain_id', ''),
                    'light_chain_id': metadata.get('light_chain_id', '')
                }
                
                sasa_data.append(row_data)
            
            return pd.DataFrame(sasa_data)
            
        finally:
            # Clean up temporary file
            if temp_pdb_path.exists():
                temp_pdb_path.unlink()
    
    def compare_sasa(self, pdb_path, antibody_chains, site_range=None):
        """Calculate SASA both with and without antigen and compare."""
        
        sasa_with = self._calculate_sasa(pdb_path, antibody_chains, with_antigen=True, site_range=site_range)
        sasa_without = self._calculate_sasa(pdb_path, antibody_chains, with_antigen=False, site_range=site_range)
        
        # Merge the results - keep all metadata columns from with_antigen
        metadata_cols = ['pdb_id', 'organism', 'abid', 'v_gene', 'j_gene', 'heavy_chain_id', 'light_chain_id']
        base_cols = ['chain_id', 'residue_number', 'insertion_code', 'amino_acid', 'ca_coordinates']
        
        comparison = pd.merge(
            sasa_with[base_cols + metadata_cols + ['sasa', 'rsa']],
            sasa_without[['chain_id', 'residue_number', 'insertion_code', 'sasa', 'rsa']],
            on=['chain_id', 'residue_number', 'insertion_code'],
            suffixes=('_with_antigen', '_without_antigen'),
            how='outer'
        )
        
        # Calculate differences
        comparison['sasa_difference'] = comparison['sasa_with_antigen'] - comparison['sasa_without_antigen']
        comparison['rsa_difference'] = comparison['rsa_with_antigen'] - comparison['rsa_without_antigen']
        
        # Calculate relative change
        comparison['sasa_relative_change'] = (
            comparison['sasa_difference'] / comparison['sasa_without_antigen']
        ).replace([np.inf, -np.inf], np.nan)
        
        comparison['rsa_relative_change'] = (
            comparison['rsa_difference'] / comparison['rsa_without_antigen']
        ).replace([np.inf, -np.inf], np.nan)
        
        return comparison
    
    def cleanup(self):
        """Remove temporary directory."""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)


def process_pdb_directory(pdb_directory, antibody_chains=['H', 'L'], site_range=None, 
                         output_file=None, max_files=None, include_metadata=True, 
                         metadata_extractor=None, verbose=False):
    """
    Process all PDB files in a directory with protein metadata integration.
    
    Parameters:
    -----------
    pdb_directory : str or Path
        Directory containing PDB files
    antibody_chains : list
        List of chain IDs for antibody chains (default: ['H', 'L'])
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
    verbose : bool
        Whether to print verbose output
        
    Returns:
    --------
    pd.DataFrame
        Combined results for all PDB files with metadata
    """
    pdb_dir = Path(pdb_directory)
    pdb_files = list(pdb_dir.glob("*.pdb"))
    
    if max_files:
        pdb_files = pdb_files[:max_files]
    
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir}")
        return pd.DataFrame()
    
    # Initialize metadata extractor if requested
    if include_metadata and metadata_extractor is None:
        metadata_extractor = ProteinMetadataExtractor(verbose=verbose)
        if verbose and metadata_extractor.sabdab_df is not None:
            print(f"Loaded metadata for {len(metadata_extractor.sabdab_df)} structures")
    
    calculator = SASACalculator(metadata_extractor=metadata_extractor if include_metadata else None)
    all_results = []
    failed_files = []
    
    print(f"Processing {len(pdb_files)} PDB files...")
    
    for pdb_file in tqdm(pdb_files, desc="Processing PDBs"):
        try:
            pdb_id = pdb_file.stem
            results = calculator.compare_sasa(str(pdb_file), antibody_chains, site_range)
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
    calculator.cleanup()
    
    if not all_results:
        print("No files were successfully processed")
        return pd.DataFrame()
    
    # Combine all results
    combined_results = pd.concat(all_results, ignore_index=True)
    
    # Reorder columns for better readability
    if include_metadata:
        column_order = [
            'pdb_id', 'organism', 'abid', 'v_gene', 'j_gene', 
            'heavy_chain_id', 'light_chain_id', 'chain_id', 
            'residue_number', 'insertion_code', 'amino_acid',
            'sasa_with_antigen', 'sasa_without_antigen', 'sasa_difference', 'sasa_relative_change',
            'rsa_with_antigen', 'rsa_without_antigen', 'rsa_difference', 'rsa_relative_change',
            'ca_coordinates'
        ]
    else:
        column_order = [
            'pdb_id', 'chain_id', 'residue_number', 'insertion_code', 'amino_acid',
            'sasa_with_antigen', 'sasa_without_antigen', 'sasa_difference', 'sasa_relative_change',
            'rsa_with_antigen', 'rsa_without_antigen', 'rsa_difference', 'rsa_relative_change',
            'ca_coordinates'
        ]
    
    # Only include columns that exist in the DataFrame
    available_columns = [col for col in column_order if col in combined_results.columns]
    combined_results = combined_results.reindex(columns=available_columns)
    
    # Print summary
    print(f"\n=== Processing Summary ===")
    print(f"Successfully processed: {len(all_results)} files")
    print(f"Failed to process: {len(failed_files)} files")
    print(f"Total residues analyzed: {len(combined_results)}")
    print(f"Unique PDB structures: {combined_results['pdb_id'].nunique()}")
    
    if include_metadata and 'v_gene' in combined_results.columns:
        v_gene_count = combined_results['v_gene'].notna().sum()
        j_gene_count = combined_results['j_gene'].notna().sum()
        print(f"Residues with V gene data: {v_gene_count}")
        print(f"Residues with J gene data: {j_gene_count}")
        print(f"Unique organisms: {combined_results['organism'].nunique()}")
    
    # Save to file if requested
    if output_file:
        combined_results.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    
    if verbose and failed_files:
        print(f"\nFailed files: {failed_files[:10]}{'...' if len(failed_files) > 10 else ''}")
    
    return combined_results


def print_sample_results(results, num_samples=5):
    """Print sample results with metadata."""
    if results.empty:
        print("No results to display")
        return
    
    print(f"\n=== Sample Results ({num_samples} rows) ===")
    
    # Select relevant columns for display
    display_cols = ['pdb_id', 'organism', 'v_gene', 'j_gene', 'chain_id', 
                   'residue_number', 'amino_acid', 'sasa_difference', 'rsa_difference']
    available_cols = [col for col in display_cols if col in results.columns]
    
    sample_data = results[available_cols].head(num_samples)
    print(sample_data.to_string(index=False))


def main():
    """Main function for command line execution."""
    parser = argparse.ArgumentParser(
        description="Analyze solvent accessibility (SASA/RSA) for PDB structures with V/D/J gene metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_sasa_analysis.py --scheme opig-imgt --output results.csv
  python run_sasa_analysis.py --scheme rcsb --max-files 10 --verbose
  python run_sasa_analysis.py --pdb-dir /path/to/pdbs --no-metadata --output basic_results.csv
        """
    )
    
    parser.add_argument('--scheme', default='opig-imgt', 
                       choices=['opig-imgt', 'opig-chothia', 'rcsb'],
                       help='PDB numbering scheme to use (default: opig-imgt)')
    
    parser.add_argument('--pdb-dir', type=str, 
                       help='Custom PDB directory path (overrides --scheme)')
    
    parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output CSV file path')
    
    parser.add_argument('--antibody-chains', nargs='+', default=['H', 'L'],
                       help='Antibody chain IDs (default: H L)')
    
    parser.add_argument('--site-range', nargs=2, type=int, metavar=('MIN', 'MAX'),
                       help='Residue number range to analyze (e.g., --site-range 1 150)')
    
    parser.add_argument('--max-files', type=int,
                       help='Maximum number of PDB files to process (for testing)')
    
    parser.add_argument('--no-metadata', action='store_true',
                       help='Disable V/D/J gene metadata extraction')
    
    parser.add_argument('--sabdab-summary', type=str,
                       default='/fh/fast/matsen_e/shared/bcr-mut-sel/sabdab/sabdab_summary_for_dnsm.tsv',
                       help='Path to SAbDab summary file')
    
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Print verbose output')
    
    args = parser.parse_args()
    
    # Set up environment
    os.environ.setdefault('LIBCIFPP_DATA_DIR', '/home/nharel/miniforge3/envs/netam_env/share/libcifpp')
    
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
        metadata_extractor = ProteinMetadataExtractor(
            sabdab_summary_path=args.sabdab_summary,
            verbose=args.verbose
        )
    
    # Print configuration
    print(f"=== SASA Analysis Configuration ===")
    print(f"PDB directory: {pdb_directory}")
    print(f"Antibody chains: {args.antibody_chains}")
    print(f"Site range: {site_range if site_range else 'All residues'}")
    print(f"Max files: {args.max_files if args.max_files else 'All files'}")
    print(f"Include metadata: {include_metadata}")
    print(f"Output file: {args.output}")
    print()
    
    # Run analysis
    try:
        results = process_pdb_directory(
            pdb_directory=pdb_directory,
            antibody_chains=args.antibody_chains,
            site_range=site_range,
            output_file=args.output,
            max_files=args.max_files,
            include_metadata=include_metadata,
            metadata_extractor=metadata_extractor,
            verbose=args.verbose
        )
        
        if not results.empty:
            print_sample_results(results, num_samples=10)
            
            # Basic statistics
            print(f"\n=== Basic Statistics ===")
            print(f"Mean SASA difference: {results['sasa_difference'].mean():.3f} ± {results['sasa_difference'].std():.3f} A^2")
            print(f"Mean RSA difference: {results['rsa_difference'].mean():.4f} ± {results['rsa_difference'].std():.4f}")
            
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