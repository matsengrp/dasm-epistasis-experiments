# dasm-epistasis-experiments

Analysis code for reproducing the results presented in:

**"Entrenchment of germline amino-acid differences in antibody affinity maturation"** by Noam Harel, Kevin Sung, Will Dumm, Mackenzie M. Johnson, David Rich, Julia Fukuyama, Hugh K. Haddox, Frederick A. Matsen IV

bioRxiv preprint: [10.64898/2026.04.21.720000](https://doi.org/10.64898/2026.04.21.720000)

## Installation

First follow the instructions to install [netam](https://github.com/matsengrp/netam) into a virtual environment.

Then install additional dependencies:

```bash
pip install aaindex Levenshtein logomaker statsmodels
pip install -e /path/to/thrifty-experiments-1  # for shmex; see https://github.com/matsengrp/thrifty-experiments-1
```

## Figure-to-source mapping

### Main figures

| Figure | Description | Source |
|--------|-------------|--------|
| Fig 1 | Entrenchment flow diagram | Manually created (`entrenchment_flow_figure.pdf`, SVG in tex repo) |
| Fig 2 | Reciprocal selection + Entrenchment within IGHV1/3 + Shannon entropy at within-family entrenched sites | [`v_families_dasm.ipynb`](v_families_dasm.ipynb) (reciprocal selection), [`shannon_entropy_entrenchment.ipynb`](shannon_entropy_entrenchment.ipynb) (entropy) |
| Fig 3 | RSA at within-family entrenched sites + structural images | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) (RSA), structural images manually created |
| Fig 4 | Grantham distance at entrenched sites | [`grantham_distance_analysis.ipynb`](grantham_distance_analysis.ipynb) |
| Fig 5 | Between-family entrenchment with within-family overlay | [`v_families_dasm.ipynb`](v_families_dasm.ipynb) |
| Fig 6 | Structural entrenchment at sites 9 and 73 | Structural images manually created (ChimeraX) |
| Fig 7 | Mutation rate validation | [`rates_analysis_productive_non_productive.ipynb`](rates_analysis_productive_non_productive.ipynb) (Rodriguez), [`rates_analysis_productive_w_thrifty_multi.ipynb`](rates_analysis_productive_w_thrifty_multi.ipynb) (Jaffe+Tang thrifty) |

### Supplementary figures

| Figure | Description | Source |
|--------|-------------|--------|
| Fig S1 | V-gene pairwise similarity | [`germline.ipynb`](germline.ipynb) |
| Fig S2 | Threshold sensitivity within IGHV1 | [`entrenchment_threshold_sensitivity.py`](entrenchment_threshold_sensitivity.py) |
| Fig S3 | Within-family validation strip plot | [`within_family_validation.ipynb`](within_family_validation.ipynb) |
| Fig S4 | Within-family entrenchment for IGHV4 | [`v_families_dasm.ipynb`](v_families_dasm.ipynb) |
| Fig S5 | RSA at IGHV1 within-family entrenched sites | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S6 | RSA grid by site and amino acid | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S7 | Selection factors at within-family entrenched sites | [`grantham_distance_analysis.ipynb`](grantham_distance_analysis.ipynb) |
| Fig S8 | RSA at between-family entrenched sites | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S9 | Shannon entropy (3-category entrenchment) | [`shannon_entropy_entrenchment.ipynb`](shannon_entropy_entrenchment.ipynb) |
| Fig S10 | Between-family entrenchment for additional V-family pairs | [`v_families_dasm.ipynb`](v_families_dasm.ipynb) |
| Fig S11 | Backbone angles at site 9 | [`solvent_accessibility_analysis.ipynb`](solvent_accessibility_analysis.ipynb) |
| Fig S12 | Pairwise validation for IGHV4 comparisons | [`rates_analysis_productive_w_thrifty_multi.ipynb`](rates_analysis_productive_w_thrifty_multi.ipynb) |

### Tables

| Table | Description | Source |
|-------|-------------|--------|
| Table S1 | Combined entrenched sites | [`create_combined_entrenched_sites.py`](create_combined_entrenched_sites.py) |
| Table S2 | Hydrogen bonding at sites 73-75 | Manually created (ChimeraX) |
| Table S3 | Combined validation table | [`create_combined_validation_table.py`](create_combined_validation_table.py) |

## Supporting code

| File | Role |
|------|------|
| [`utils.py`](utils.py) | Core utilities imported by most notebooks |
| [`rates_analysis_util.py`](rates_analysis_util.py) | Rate analysis utilities for validation notebooks |
| [`sasa_plotting.py`](sasa_plotting.py) | RSA/SASA plotting functions |
| [`run_sasa_analysis.py`](run_sasa_analysis.py) | SASA data generation from PDB structures |
| [`create_germline_codon_tables.py`](create_germline_codon_tables.py) | Generates germline reference data |
| [`run_thrifty_neutral.py`](run_thrifty_neutral.py) | Generates cached neutral rate data for thrifty validation |
| [`site9_discrepancy_analysis.ipynb`](site9_discrepancy_analysis.ipynb) | Investigation of site 9 validation discrepancy |
| [`dnsmex/`](dnsmex/) | Core library for DNSM/DASM analysis |
| [`germline/`](germline/) | Germline reference data (OGRDB-derived) |
| `_output/entrenchment_analysis/` | Intermediate entrenchment results consumed by notebooks |

## Exploratory notebooks (not in paper)

| Notebook | Description |
|----------|-------------|
| [`grab_motifs_with_dnsm.ipynb`](grab_motifs_with_dnsm.ipynb) | Motif extraction, may generate intermediate data |
| [`dasm_model_comparison.ipynb`](dasm_model_comparison.ipynb) | Model comparison exploration |
| [`light_chain_selection_at_entrenched_sites.ipynb`](light_chain_selection_at_entrenched_sites.ipynb) | Light chain selection analysis |
| [`light_chain_pairing_bias_with_entrenched_sites.ipynb`](light_chain_pairing_bias_with_entrenched_sites.ipynb) | Light chain pairing bias |
| [`neutral_rates_for_interesting_sites.ipynb`](neutral_rates_for_interesting_sites.ipynb) | Neutral rates at specific sites |
