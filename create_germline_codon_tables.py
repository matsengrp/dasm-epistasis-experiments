"""
Create germline codon tables for IMGT and Chothia numbering schemes.

This script downloads OGRDB germline sequences and creates codon-level dataframes
with site numbering for both IMGT and Chothia schemes.

Purpose:
--------
We need germline annotation files for:
- In current analysis repo:
    - Add germline annotation per site (is_germline_codon/aa)
    - Annotate sites for out of frame data - ANARCI cannot be used to number these
      sequences, so for the V gene area (until CDR3) for sequences that have no
      frameshifts (filtered for earlier), we can use the germline annotation files
      to annotate the sequences according to PARTIS v gene identification
- Later this can be used for pcp generation files (pcp-pipeline) with Chothia-based
  CDRs start and end sites.

Approach:
---------
IMGT gapped files can be downloaded from the OGRDB website.
There is no such resource for Chothia.
We create Chothia numbering with ANARCI:
- Download the ungapped version from OGRDB website
- Filter for IGHV
- Run in ANARCI with human and heavy chain parameters
- Note: Chothia uses insertion codes (e.g., 52A, 52B) while IMGT is just int-based
  for IGHV gene area

Output files:
-------------
- germline/germline_codons_imgt.csv: IMGT numbering scheme
- germline/germline_codons_chothia.csv: Chothia numbering scheme
"""

import os
import sys
import pandas as pd
import requests
from Bio import SeqIO
import subprocess
from collections import Counter

# Add parent directory to path to import from dnsmex
sys.path.insert(0, '/home/nharel/re/dnsm-experiments-1')
from netam.sequences import translate_sequence


def download_ogrdb_germline_fasta(output_dir="germline"):
    """Download OGRDB germline sequences (gapped and ungapped versions)."""
    os.makedirs(output_dir, exist_ok=True)

    # Download gapped version (for IMGT)
    gapped_url = "https://ogrdb.airr-community.org/download_germline_set/Human/IGH_VDJ/published/gapped_ex"
    gapped_file = os.path.join(output_dir, "ogrdb_human_IGH_gapped.fasta")

    if os.path.exists(gapped_file):
        print(f"Gapped FASTA file {gapped_file} already exists. Skipping download.")
    else:
        print("Downloading gapped FASTA file from OGRDB...")
        response = requests.get(gapped_url)
        response.raise_for_status()
        with open(gapped_file, "wb") as f:
            f.write(response.content)
        print(f"Gapped FASTA downloaded successfully as {gapped_file}")

    # Download ungapped version (for Chothia with ANARCI)
    ungapped_url = "https://ogrdb.airr-community.org/download_germline_set/Human/IGH_VDJ/published/ungapped_ex"
    ungapped_file = os.path.join(output_dir, "ogrdb_human_IGH_ungapped.fasta")

    if os.path.exists(ungapped_file):
        print(f"Ungapped FASTA file {ungapped_file} already exists. Skipping download.")
    else:
        print("Downloading ungapped FASTA file from OGRDB...")
        response = requests.get(ungapped_url)
        response.raise_for_status()
        with open(ungapped_file, "wb") as f:
            f.write(response.content)
        print(f"Ungapped FASTA downloaded successfully as {ungapped_file}")

    return gapped_file, ungapped_file


def filter_ighv_sequences(fasta_file, output_file):
    """Filter FASTA file to only include IGHV sequences."""
    print(f"Filtering IGHV sequences from {fasta_file}...")
    SeqIO.write(
        (record for record in SeqIO.parse(fasta_file, 'fasta') if record.id.startswith('IGHV')),
        output_file,
        'fasta'
    )
    count = sum(1 for _ in SeqIO.parse(output_file, 'fasta'))
    print(f"Filtered {count} IGHV sequences saved to {output_file}")
    return output_file


def translate_to_amino_acids(input_fasta, output_fasta):
    """Translate nucleotide sequences to amino acids."""
    print(f"Translating sequences to amino acids...")
    with open(output_fasta, 'w') as f:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            # Trim to multiple of 3
            seq_len = len(record.seq)
            trimmed_len = seq_len - (seq_len % 3)
            trimmed_seq = record.seq[:trimmed_len]

            # Translate
            aa_seq = trimmed_seq.translate()
            f.write(f'>{record.id}\n{str(aa_seq)}\n')

    count = sum(1 for _ in SeqIO.parse(output_fasta, 'fasta'))
    print(f"Translated {count} sequences saved to {output_fasta}")
    return output_fasta


def run_anarci(input_fasta, output_prefix, numbering_scheme, custom_load_bash_command='source ~/.bashrc && conda activate netam_env'):
    """Run ANARCI to number sequences."""
    print(f"\nRunning ANARCI with {numbering_scheme} numbering...")

    # Construct the full command with custom environment loading
    anarci_command = (
        f'ANARCI -i {input_fasta} -o {output_prefix} '
        f'-s {numbering_scheme} -r heavy --use_species human --csv --assign_germline'
    )

    if custom_load_bash_command != None:
        anarci_command = f'{custom_load_bash_command} && {anarci_command}'

    result = subprocess.run(
        anarci_command,
        shell=True,
        executable='/bin/bash',
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"ANARCI failed with exit code {result.returncode}")
        print(f"Stderr: {result.stderr}")
        return None

    # ANARCI creates files with _H.csv suffix for heavy chain
    output_file = f"{output_prefix}_H.csv"

    if os.path.exists(output_file):
        print(f"ANARCI completed successfully! Output: {output_file}")
        return output_file
    else:
        print(f"ANARCI output file not found: {output_file}")
        return None


def verify_codon_gaps(fasta_file):
    """Verify that gaps in the gapped FASTA are codon-aligned (3 nucleotides per gap)."""
    print("Verifying that gaps are codon-aligned...")
    all_valid = True
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id.startswith('IGHV'):
            seq = str(record.seq)
            counts = Counter([i // 3 for i in range(len(seq)) if seq[i] == '.'])
            # Each codon position should have exactly 3 gaps (one per nucleotide)
            if not all(count == 3 for count in counts.values()):
                print(f"WARNING: {record.id} has non-codon-aligned gaps!")
                all_valid = False

    if all_valid:
        print("✓ All gaps are codon-aligned (3 nucleotides per codon)")
    return all_valid


def create_imgt_codon_table(gapped_fasta, output_csv):
    """Create IMGT codon table from gapped FASTA file."""
    print("\n=== Creating IMGT codon table ===")

    # Verify gaps are codon-aligned
    verify_codon_gaps(gapped_fasta)

    # Read FASTA and convert to DataFrame
    records = []
    for record in SeqIO.parse(gapped_fasta, "fasta"):
        if record.id.startswith('IGHV'):
            records.append({
                'v_gene': record.id,
                'sequence': str(record.seq),
                'v_family': record.id.split('-')[0]
            })

    df = pd.DataFrame(records)
    print(f"Loaded {len(df)} IGHV sequences")

    # Create codon dataframe
    codons_df = []
    for row in df.itertuples():
        sequence = row.sequence
        seq_df = pd.DataFrame(
            list(enumerate([sequence[i:i+3] for i in range(0, len(sequence), 3)], start=1)),
            columns=['site', 'codon']
        )
        seq_df['v_gene'] = row.v_gene
        seq_df['v_family'] = row.v_family
        codons_df.append(seq_df)

    codons_df = pd.concat(codons_df, ignore_index=True)

    # Filter to ensure codons are of length 3 and remove gaps
    codons_df = codons_df[(codons_df['codon'] != '...') & (codons_df['codon'].str.len() == 3)]

    # Remove v genes that have stop codons in them - likely pseudogenes
    stop_codons = ['TAA', 'TAG', 'TGA']
    genes_with_stops = codons_df[codons_df['codon'].isin(stop_codons)]['v_gene'].unique()
    print(f"Removing {len(genes_with_stops)} genes with stop codons")
    codons_df = codons_df[~codons_df['v_gene'].isin(genes_with_stops)]

    # Translate codons to amino acids
    codons_df['amino_acid'] = codons_df['codon'].apply(translate_sequence)

    # Save to CSV
    codons_df.to_csv(output_csv, index=False)
    print(f"IMGT codon table saved: {output_csv}")
    print(f"  Rows: {len(codons_df)}")
    print(f"  Unique genes: {codons_df['v_gene'].nunique()}")

    return codons_df


def create_chothia_codon_table(anarci_csv, ungapped_fasta, output_csv):
    """Create Chothia codon table from ANARCI output and ungapped FASTA."""
    print("\n=== Creating Chothia codon table ===")

    # Read ANARCI results
    anarci_df = pd.read_csv(anarci_csv)
    print(f"Loaded ANARCI results: {len(anarci_df)} sequences")

    # Read ungapped FASTA to get nucleotide sequences
    seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(ungapped_fasta, 'fasta')}

    # Get site column names (numeric positions, skipping metadata columns)
    metadata_cols = ['Id', 'domain_no', 'hmm_species', 'chain_type', 'e-value', 'score',
                     'seqstart_index', 'seqend_index', 'identity_species', 'v_gene',
                     'v_identity', 'j_gene', 'j_identity']
    site_cols = [col for col in anarci_df.columns if col not in metadata_cols]

    # Build long-format dataframe by zipping ANARCI positions with nucleotides
    rows = []
    initial_gene_count = len(anarci_df)

    for idx, row in anarci_df.iterrows():
        v_gene = row['Id']  # Use Id from OGRDB instead of v_gene ANARCI assignment column
        v_family = v_gene.split('-')[0] if pd.notna(v_gene) else None

        # Get nucleotide sequence
        if v_gene in seq_dict:
            nuc_seq = seq_dict[v_gene]
            # Trim to multiple of 3
            nuc_seq = nuc_seq[:len(nuc_seq) - (len(nuc_seq) % 3)]

            # Track nucleotide position as we iterate through ANARCI sites
            nuc_idx = 0

            for site in site_cols:
                amino_acid = row[site]

                # If there's an amino acid at this site (not a gap)
                if amino_acid != '-' and pd.notna(amino_acid):
                    # Extract codon from nucleotide sequence
                    codon_start = nuc_idx * 3
                    codon_end = codon_start + 3

                    if codon_end <= len(nuc_seq):
                        codon = nuc_seq[codon_start:codon_end]
                    else:
                        codon = None

                    rows.append({
                        'site': site,
                        'codon': codon,
                        'v_gene': v_gene,
                        'v_family': v_family,
                        'amino_acid': amino_acid
                    })

                    # Move to next codon in nucleotide sequence
                    nuc_idx += 1
                # If it's a gap, skip this site (don't increment nuc_idx)

    # Create dataframe
    result_df = pd.DataFrame(rows)
    print(f"Created dataframe with {len(result_df)} rows from {initial_gene_count} sequences")

    # Remove v genes that have stop codons in them - likely pseudogenes
    stop_codons = ['TAA', 'TAG', 'TGA']
    genes_with_stops = result_df[result_df['codon'].isin(stop_codons)]['v_gene'].unique()
    print(f"Removing {len(genes_with_stops)} genes with stop codons")
    result_df = result_df[~result_df['v_gene'].isin(genes_with_stops)]

    # Remove v genes that are not aligned until the end of the v gene
    # Remove anything where the last site is not '94' (Chothia numbering)
    last_sites = result_df.groupby('v_gene')['site'].max()
    genes_not_to_94 = last_sites[last_sites != '94'].index.tolist()
    print(f"Removing {len(genes_not_to_94)} genes not extending to site 94")
    result_df = result_df[~result_df['v_gene'].isin(genes_not_to_94)]

    # Save to CSV
    result_df.to_csv(output_csv, index=False)
    print(f"Chothia codon table saved: {output_csv}")
    print(f"  Rows: {len(result_df)}")
    print(f"  Unique genes: {result_df['v_gene'].nunique()}")

    return result_df


def main():
    """Main function to create both IMGT and Chothia codon tables."""
    output_dir = "germline"

    print("="*80)
    print("Creating Germline Codon Tables")
    print("="*80)

    # Step 1: Download OGRDB files
    print("\n### Step 1: Downloading OGRDB files ###")
    gapped_fasta, ungapped_fasta = download_ogrdb_germline_fasta(output_dir)

    # Step 2: Create IMGT codon table
    print("\n### Step 2: Creating IMGT codon table ###")
    imgt_output = os.path.join(output_dir, "germline_codons_imgt.csv")
    imgt_df = create_imgt_codon_table(gapped_fasta, imgt_output)

    # Step 3: Prepare files for Chothia (ANARCI)
    print("\n### Step 3: Preparing files for Chothia numbering ###")

    # Filter ungapped FASTA to IGHV only
    ungapped_ighv = os.path.join(output_dir, "ogrdb_human_IGH_ungapped_ighv.fasta")
    filter_ighv_sequences(ungapped_fasta, ungapped_ighv)

    # Translate to amino acids for ANARCI
    ighv_aa = os.path.join(output_dir, "ogrdb_human_IGH_ungapped_ighv_aa.fasta")
    translate_to_amino_acids(ungapped_ighv, ighv_aa)

    # Step 4: Run ANARCI with Chothia numbering
    print("\n### Step 4: Running ANARCI with Chothia numbering ###")
    anarci_output_prefix = os.path.join(output_dir, "ogrdb_human_IGH_ungapped_ighv_chothia")
    anarci_csv = run_anarci(ighv_aa, anarci_output_prefix, "chothia")

    if anarci_csv is None:
        print("ERROR: ANARCI failed. Cannot create Chothia codon table.")
        return

    # Step 5: Create Chothia codon table
    print("\n### Step 5: Creating Chothia codon table ###")
    chothia_output = os.path.join(output_dir, "germline_codons_chothia.csv")
    chothia_df = create_chothia_codon_table(anarci_csv, ungapped_ighv, chothia_output)

    # Final summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"IMGT codon table:")
    print(f"  File: {imgt_output}")
    print(f"  Genes: {imgt_df['v_gene'].nunique()}")
    print(f"  Rows: {len(imgt_df)}")
    print(f"\nChothia codon table:")
    print(f"  File: {chothia_output}")
    print(f"  Genes: {chothia_df['v_gene'].nunique()}")
    print(f"  Rows: {len(chothia_df)}")
    print("\nDone!")


if __name__ == "__main__":
    main()
