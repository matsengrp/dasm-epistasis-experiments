# dasm-epistasis-experiments

Analysis code for reproducing the results presented in:

**"Entrenchment of germline amino-acid differences in antibody affinity maturation"** by Noam Harel, Kevin Sung, Will Dumm, Mackenzie M. Johnson, David Rich, Julia Fukuyama, Hugh K. Haddox, Frederick A. Matsen IV

bioRxiv preprint: [10.64898/2026.04.21.720000](https://doi.org/10.64898/2026.04.21.720000)

## Installation

### 1. Install netam

Follow the instructions to install [netam](https://github.com/matsengrp/netam) into a virtual environment.

### 2. Install additional dependencies

```bash
pip install aaindex Levenshtein logomaker statsmodels
```

### 3. Install shmex

The rate validation analysis requires [shmex](https://github.com/matsengrp/thrifty-experiments-1) from the thrifty-experiments-1 repository:

```bash
git clone https://github.com/matsengrp/thrifty-experiments-1.git
cd thrifty-experiments-1
pip install -e .
```

## Data

All data required to reproduce the analyses are available as a single archive on Zenodo.

### 1. Download and extract data

Download `dasm-epistasis-data.tar` from [Zenodo](https://doi.org/10.5281/zenodo.20171634) and extract it:

```bash
tar xf dasm-epistasis-data.tar
```

This creates a `dasm-epistasis-data/` directory containing:

| Directory | Contents |
|-----------|----------|
| `v3/` | PCP files (Jaffe, Tang, Rodriguez) and Chothia-numbered ANARCI outputs |
| `v1/` | Out-of-frame PCP dataset (Tang-SHM, used by shmex) |
| `trained_models/` | DASM (4M params) and DNSM (1M params) model weights and branch lengths |
| `dasm_test_output/` | Cached DASM evaluation results on Rodriguez (see note below) |
| `dnsm_test_output/` | Cached observed mutation datasets on Jaffe, Tang, and Rodriguez, computed via the DNSM framework for convenience (see note below) |
| `neutral_mutability_cache/` | Pre-computed neutral mutability DataFrames for Jaffe, Tang, and Rodriguez (see note below) |

**Note on cached outputs:** The `dasm_test_output/`, `dnsm_test_output/`, and `neutral_mutability_cache/` directories contain pre-computed intermediate results that are included for convenience. The analysis code will regenerate these automatically if they are missing, but recomputation is expensive. Only the PCP files, ANARCI outputs, and trained models are strictly required.

### 2. Configure local paths

```bash
cp dnsmex/local_config.py.template dnsmex/local_config.py
```

Edit `dnsmex/local_config.py` if needed. The default assumes the archive was extracted to `~/dasm-epistasis-data`; update `ZENODO_DATA_DIR` if you extracted it elsewhere.

### 3. Set up shmex data path

The [`rates_analysis_productive_non_productive.ipynb`](rates_analysis_productive_non_productive.ipynb) notebook uses the [shmex](https://github.com/matsengrp/thrifty-experiments-1) package, which expects data at `~/data/v1/`. Create a symlink from the extracted data:

```bash
mkdir -p ~/data
ln -s /path/to/dasm-epistasis-data/v1 ~/data/v1
```

### Solvent Accessible Surface Area (SASA) data (optional)

Pre-computed SASA results are included in the repository (`_output/sasa_human_chothia_anarci.csv.gz`) and are loaded automatically by [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb). Most users do not need to regenerate this data.

If you would like to regenerate this data from scratch, there are two steps:

#### Step 1: Generate V/J gene annotations with ANARCI

[`annotate_sabdab_anarci.py`](annotate_sabdab_anarci.py) assigns V/J germline genes to SAbDab structures using [ANARCI](https://github.com/oxpig/ANARCI). It reads the SAbDab bulk download TSV, extracts heavy and light chain sequences from PDB files using BioPython, and runs ANARCI with species-constrained germline assignment.

```bash
pip install anarci
python annotate_sabdab_anarci.py --output data/sabdab/sabdab_anarci_annotations.tsv
```

Pre-computed annotations are included in [`data/sabdab/sabdab_anarci_annotations.tsv`](data/sabdab/sabdab_anarci_annotations.tsv), so this step can be skipped unless you want to regenerate them.

#### Step 2: Run the SASA analysis pipeline

[`run_sasa_analysis.py`](run_sasa_analysis.py) computes per-residue solvent accessibility across six structural scenarios (full complex, antibody alone, heavy chain alone, etc.) and derives burial effects from each binding partner. It applies a 6-stage filtering pipeline: file discovery, organism filtering (exact match), chain completeness, heavy-chain-only antibody exclusion, V/J gene consistency, and optional file limits.

You will need:

- **SAbDab summary tables and ANARCI annotations**: included in [`data/sabdab/`](data/sabdab/) (2024-01-26 snapshot)
- **PDB structure files**: download from [RCSB PDB](https://www.rcsb.org/) or via SAbDab's bulk download. Update `PDB_BASE_DIR` in [`run_sasa_analysis.py`](run_sasa_analysis.py) to point to your local PDB directory.
- **DSSP**: install via `conda install -c conda-forge dssp`

```bash
python run_sasa_analysis.py --organism "homo sapiens"
```

## Reproducing the analysis

### Step 1: Run main analysis notebooks

Run [`v_families_entrenchment_dasm.ipynb`](v_families_entrenchment_dasm.ipynb) first — it produces Figs 2, 5, S6, S12, and writes entrenchment results to `_output/entrenchment_analysis/` that are required by the other notebooks.

The remaining notebooks can then be run in any order:

- [`shannon_entropy_entrenchment.ipynb`](shannon_entropy_entrenchment.ipynb) — Shannon entropy at entrenched vs. non-entrenched sites (Figs 2C, S11)
- [`grantham_distance_analysis.ipynb`](grantham_distance_analysis.ipynb) — Physicochemical distance (Grantham) of entrenched substitutions (Figs 4, S9)
- [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) — Relative solvent accessibility (RSA) and partner contact analysis at entrenched sites (Figs 3, 5C, S7, S8, S10, S13)
- [`germline.ipynb`](germline.ipynb) — V-gene pairwise amino acid similarity (Fig S1)
- [`within_family_validation.ipynb`](within_family_validation.ipynb) — Validates that pooling V gene alleles within a family does not create false entrenchment calls (Fig S5)
- [`rates_analysis_productive_non_productive.ipynb`](rates_analysis_productive_non_productive.ipynb) — Mutation rate validation using out-of-frame sequences as neutral baseline (Fig 7A)
- [`rates_analysis_productive_w_thrifty_multi.ipynb`](rates_analysis_productive_w_thrifty_multi.ipynb) — Mutation rate validation using Thrifty-predicted neutral rates as baseline (Figs 7B-C, S14). **Note:** This notebook loads a large dataframe so run with enough memory.

### Step 2: Aggregate results and create supplementary data files

```bash
python create_combined_entrenched_sites.py    # combined_entrenched_sites.csv
python entrenchment_threshold_sensitivity.py  # Figs S2-S4
python create_combined_validation_table.py    # combined_validation_table_entrenched_for_paper.csv
```

## Figure-to-source mapping

### Main figures

| Figure | Panel | Description | Source |
|--------|-------|-------------|--------|
| Fig 1 | | Entrenchment flow diagram | Manually created (`entrenchment_flow_figure.pdf`, SVG in tex repo) |
| Fig 2 | A,B | Reciprocal selection + entrenchment within IGHV1/3 | [`v_families_entrenchment_dasm.ipynb`](v_families_entrenchment_dasm.ipynb) |
| Fig 2 | C | Shannon entropy at within-family entrenched sites | [`shannon_entropy_entrenchment.ipynb`](shannon_entropy_entrenchment.ipynb) |
| Fig 3 | A–C | Structural images | Manually created (ChimeraX) |
| Fig 3 | D | RSA at within-family entrenched sites | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig 4 | | Grantham distance at entrenched sites | [`grantham_distance_analysis.ipynb`](grantham_distance_analysis.ipynb) |
| Fig 5 | A,B | Between-family entrenchment with within-family overlay | [`v_families_entrenchment_dasm.ipynb`](v_families_entrenchment_dasm.ipynb) |
| Fig 5 | C | RSA at between-family entrenched sites | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig 6 | | Structural entrenchment at sites 9 and 73 | Manually created (ChimeraX) |
| Fig 7 | A | Mutation rate validation using out-of-frame baseline | [`rates_analysis_productive_non_productive.ipynb`](rates_analysis_productive_non_productive.ipynb) |
| Fig 7 | B,C | Mutation rate validation using Thrifty baseline | [`rates_analysis_productive_w_thrifty_multi.ipynb`](rates_analysis_productive_w_thrifty_multi.ipynb) |

### Supplementary figures and tables

| Figure | Description | Source |
|--------|-------------|--------|
| Fig S1 | V-gene pairwise similarity | [`germline.ipynb`](germline.ipynb) |
| Fig S2 | Threshold sensitivity within IGHV1 | [`entrenchment_threshold_sensitivity.py`](entrenchment_threshold_sensitivity.py) |
| Fig S3 | Threshold sensitivity within IGHV3 | [`entrenchment_threshold_sensitivity.py`](entrenchment_threshold_sensitivity.py) |
| Fig S4 | Threshold sensitivity IGHV1 vs IGHV3 | [`entrenchment_threshold_sensitivity.py`](entrenchment_threshold_sensitivity.py) |
| Fig S5 | Per-V-gene-allele median selection factors at entrenched sites | [`within_family_validation.ipynb`](within_family_validation.ipynb) |
| Fig S6 | Within-family entrenchment for IGHV4 | [`v_families_entrenchment_dasm.ipynb`](v_families_entrenchment_dasm.ipynb) |
| Fig S7 | RSA at IGHV1 within-family entrenched sites | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S8 | RSA grid by site and amino acid | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S9 | Reciprocal selection factors at within-family entrenched sites | [`grantham_distance_analysis.ipynb`](grantham_distance_analysis.ipynb) |
| Fig S10 | RSA at between-family entrenched sites | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S11 | Shannon entropy (3-category entrenchment) | [`shannon_entropy_entrenchment.ipynb`](shannon_entropy_entrenchment.ipynb) |
| Fig S12 | Between-family entrenchment for additional V-family pairs | [`v_families_entrenchment_dasm.ipynb`](v_families_entrenchment_dasm.ipynb) |
| Fig S13 | Backbone angles at site 9 | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S14 | Pairwise validation for IGHV4 comparisons | [`rates_analysis_productive_w_thrifty_multi.ipynb`](rates_analysis_productive_w_thrifty_multi.ipynb) |
| Table S1 | Hydrogen bonding at sites 73–75 | Manually created (ChimeraX); see [`site_73_structural/`](site_73_structural/) for commands and notes |

### Supplementary data files

| File | Description |
|------|-------------|
| [`_output/combined_entrenched_sites.csv`](_output/combined_entrenched_sites.csv) | Combined list of all entrenched amino acid substitutions identified from both within-family and between-family comparisons. Each row specifies the V family, site, parent and target amino acid, the source comparison, and whether it is a within- or between-family result. Generated by [`create_combined_entrenched_sites.py`](create_combined_entrenched_sites.py). |
| [`_output/combined_validation_table_entrenched_for_paper.csv`](_output/combined_validation_table_entrenched_for_paper.csv) | Summary of all entrenched sites with DASM selection factors and mutation rate ratio validation. Each row represents one entrenched substitution (parent amino acid → child amino acid) at a given site and V family. Columns include: DASM median log selection factor; out-of-frame validation columns (prefix `oof_`) with pseudocount-adjusted mutation counts, total branch lengths, and log rate ratios; Thrifty validation columns (prefix `thrifty_`) with pseudocount-adjusted counts and log ratios. Generated by [`create_combined_validation_table.py`](create_combined_validation_table.py). See the paper methods for details. |

## Supporting code

| File | Role |
|------|------|
| [`utils.py`](utils.py) | Core utilities imported by most notebooks |
| [`rates_analysis_util.py`](rates_analysis_util.py) | Rate analysis utilities for validation notebooks |
| [`sasa_plotting.py`](sasa_plotting.py) | RSA/SASA plotting functions |
| [`annotate_sabdab_anarci.py`](annotate_sabdab_anarci.py) | V/J gene annotation of SAbDab structures using ANARCI |
| [`run_sasa_analysis.py`](run_sasa_analysis.py) | SASA data generation from PDB structures |
| [`create_germline_codon_tables.py`](create_germline_codon_tables.py) | Generates germline reference data |
| [`site9_discrepancy_analysis.ipynb`](site9_discrepancy_analysis.ipynb) | Investigation of site 9 validation discrepancy |
| [`dnsmex/`](dnsmex/) | Core library for DNSM/DASM analysis |
| [`germline/`](germline/) | Germline reference data (OGRDB-derived); to regenerate, run [`create_germline_codon_tables.py`](create_germline_codon_tables.py) (requires [ANARCI](https://github.com/oxpig/ANARCI)) |
| [`branch-length-regression/`](branch-length-regression/) | Linear relationship between DASM branch lengths and mutation frequency, used in the Thrifty validation |

## Exploratory notebooks (not used in paper)

| Notebook | Description |
|----------|-------------|
| [`grab_motifs_with_dnsm.ipynb`](grab_motifs_with_dnsm.ipynb) | Motif extraction, may generate intermediate data |
| [`dasm_model_comparison.ipynb`](dasm_model_comparison.ipynb) | Model comparison exploration |
| [`light_chain_selection_at_entrenched_sites.ipynb`](light_chain_selection_at_entrenched_sites.ipynb) | Light chain selection analysis |
| [`light_chain_pairing_bias_with_entrenched_sites.ipynb`](light_chain_pairing_bias_with_entrenched_sites.ipynb) | Light chain pairing bias |
| [`neutral_rates_for_interesting_sites.ipynb`](neutral_rates_for_interesting_sites.ipynb) | Neutral rates at specific sites |

The [`old/`](old/) directory contains superseded analyses that are not used in the paper.
