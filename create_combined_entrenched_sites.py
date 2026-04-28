"""
Combine all between-family and within-family entrenched site CSVs into a single table.
"""

import glob
import os
import pandas as pd

NUMBERING_SCHEME = 'chothia'
DATA_DIR = f'_output/entrenchment_analysis/{NUMBERING_SCHEME}'
OUTPUT_DIR = '_output'


def load_and_combine():
    within_files = sorted(glob.glob(f'{DATA_DIR}/entrenched_aa_sites_within_*.csv'))
    vs_files = sorted(glob.glob(f'{DATA_DIR}/entrenched_aa_sites_*_vs_*.csv'))

    dfs = []
    for f in within_files:
        df = pd.read_csv(f, dtype={'site': str})
        name = os.path.basename(f).replace('entrenched_aa_sites_', '').replace('.csv', '')
        df['source'] = name
        df['type'] = 'within'
        dfs.append(df)

    for f in vs_files:
        df = pd.read_csv(f, dtype={'site': str})
        name = os.path.basename(f).replace('entrenched_aa_sites_', '').replace('.csv', '')
        df['source'] = name
        df['type'] = 'between'
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(within_files)} within-family and {len(vs_files)} between-family files")
    print(f"  Within rows: {(combined['type'] == 'within').sum()}")
    print(f"  Between rows: {(combined['type'] == 'between').sum()}")
    print(f"  Total rows: {len(combined)}")
    return combined


if __name__ == '__main__':
    combined = load_and_combine()
    output_path = f'{OUTPUT_DIR}/combined_entrenched_sites.csv'
    combined.to_csv(output_path, index=False)
    print(f"\nSaved to: {output_path}")
