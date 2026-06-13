"""
Test: does DASM predict different selection factors when the parent amino acid
is germline-encoded vs somatically mutated?

At entrenched sites, compare:
  Group A (germline): parent has amino acid X and it matches the germline codon
  Group B (mutated-to): parent has amino acid X but arrived there via somatic mutation
                        (germline at this site was a different amino acid)

If selection factors are similar between groups, then germline misassignment
doesn't affect our entrenchment signal — DASM predicts the same regardless of
whether the parent amino acid is "truly germline" or acquired via mutation.

Follows the same filtering as v_families_entrenchment_dasm.ipynb:
  - depth == 2 (MRCA sequences)
  - one_mutation_away == True (single-nucleotide substitutions)
  - ≥10 observations per group

Usage:
    cd dasm-epistasis-experiments
    python germline_vs_mutated_selection_test.py
"""

import pandas as pd
import numpy as np
from scipy import stats

from utils import load_and_process_dasm_data

MODEL_NAME = "dasm_4m-v1jaffeCC+v1tangCC-joint"
DATASET_NAME = "v1rodriguez"
NUMBERING_SCHEME = "chothia"
ENTRENCHED_DIR = "_output/entrenchment_analysis/chothia"
WITHIN_FAMILIES = ["IGHV1", "IGHV3", "IGHV4"]
MIN_OBS = 10


def load_entrenched_pairs():
    """Load entrenched (site, amino_acid, target_amino_acid) per family."""
    pairs = {}
    for fam in WITHIN_FAMILIES:
        path = f"{ENTRENCHED_DIR}/entrenched_aa_sites_within_{fam}.csv"
        df = pd.read_csv(path, dtype={"site": str})
        pairs[fam] = df
    return pairs


def main():
    print("Loading data (this takes a minute)...")
    site_df, pcp_df, aa_df = load_and_process_dasm_data(
        model_name=MODEL_NAME,
        dataset_name=DATASET_NAME,
        numbering_scheme=NUMBERING_SCHEME,
    )
    print(f"  Loaded {len(aa_df)} rows")

    entrenched_pairs = load_entrenched_pairs()

    # Filter to depth=2 and one_mutation_away, matching the notebook
    df = aa_df[aa_df["depth"] == 2].copy()
    if "one_mutation_away" in df.columns:
        df = df[df["one_mutation_away"] == True]
    df = df.dropna(subset=["log_selection_factor"])

    # Split into germline and non-germline parent
    df_germline = df[df["is_germline_codon"] == True].copy()
    df_mutated = df[df["is_germline_codon"] == False].copy()

    print(f"  depth=2 after filtering: {len(df)} rows")
    print(f"    germline codon: {len(df_germline)}")
    print(f"    non-germline codon: {len(df_mutated)}")

    # For each entrenched (v_family, site, amino_acid, target_amino_acid):
    # - Group A: germline codon, parent_aa == amino_acid, target = target_amino_acid
    # - Group B: non-germline codon, parent_aa == amino_acid, target = target_amino_acid
    #   (this means the parent mutated TO this amino acid from a different germline)

    print(f"\n{'=' * 90}")
    print("GERMLINE vs SOMATICALLY-MUTATED PARENT: SELECTION FACTOR COMPARISON")
    print(f"{'=' * 90}")
    print(
        "\nFor each entrenched substitution X→Y, compare DASM selection factors when:\n"
        "  Group A: X is the germline amino acid (current analysis)\n"
        "  Group B: X was acquired via somatic mutation (germline was something else)\n"
    )

    results = []

    for v_family in WITHIN_FAMILIES:
        ent_df = entrenched_pairs[v_family]

        print(f"\n{'─' * 80}")
        print(f"{v_family}")
        print(f"{'─' * 80}")

        for _, row in ent_df.iterrows():
            site = str(row["site"])
            parent_aa = row["amino_acid"]
            target_aa = row["target_amino_acid"]

            # Group A: germline parent
            mask_a = (
                (df_germline["v_family"] == v_family)
                & (df_germline["site"] == site)
                & (df_germline["parent_aa"] == parent_aa)
                & (df_germline["selection_factor_target_aa"] == target_aa)
            )
            group_a = df_germline.loc[mask_a, "log_selection_factor"]

            # Group B: somatically mutated parent (same amino acid, not germline)
            mask_b = (
                (df_mutated["v_family"] == v_family)
                & (df_mutated["site"] == site)
                & (df_mutated["parent_aa"] == parent_aa)
                & (df_mutated["selection_factor_target_aa"] == target_aa)
            )
            group_b = df_mutated.loc[mask_b, "log_selection_factor"]

            n_a = len(group_a)
            n_b = len(group_b)
            median_a = group_a.median() if n_a > 0 else np.nan
            median_b = group_b.median() if n_b > 0 else np.nan

            # Mann-Whitney U test if both groups have enough observations
            if n_a >= MIN_OBS and n_b >= MIN_OBS:
                stat, pval = stats.mannwhitneyu(group_a, group_b, alternative="two-sided")
                diff = median_b - median_a
                status = "TESTABLE"
            else:
                pval = np.nan
                diff = np.nan
                status = "too few obs"

            results.append({
                "v_family": v_family,
                "site": site,
                "parent_aa": parent_aa,
                "target_aa": target_aa,
                "substitution": f"{parent_aa}→{target_aa}",
                "n_germline": n_a,
                "n_mutated": n_b,
                "median_germline": median_a,
                "median_mutated": median_b,
                "diff_mutated_minus_germline": diff,
                "mannwhitney_p": pval,
                "status": status,
            })

            if n_b > 0:
                p_str = f"p={pval:.3g}" if not np.isnan(pval) else "insufficient n"
                print(
                    f"  site {site:>3} {parent_aa}→{target_aa}: "
                    f"germline median={median_a:+.2f} (n={n_a}), "
                    f"mutated median={median_b:+.2f} (n={n_b}), "
                    f"diff={diff:+.3f}, {p_str}"
                )
            else:
                print(
                    f"  site {site:>3} {parent_aa}→{target_aa}: "
                    f"germline median={median_a:+.2f} (n={n_a}), "
                    f"mutated n=0 (no sequences mutated TO {parent_aa} at this site)"
                )

    results_df = pd.DataFrame(results)
    testable = results_df[results_df["status"] == "TESTABLE"]

    print(f"\n\n{'=' * 90}")
    print("SUMMARY")
    print(f"{'=' * 90}")
    print(f"\nTotal entrenched substitutions: {len(results_df)}")
    print(f"Testable (≥{MIN_OBS} in both groups): {len(testable)}")

    if not testable.empty:
        print(f"\nAmong testable substitutions:")
        print(f"  Median absolute difference: {testable['diff_mutated_minus_germline'].abs().median():.3f}")
        print(f"  Max absolute difference:    {testable['diff_mutated_minus_germline'].abs().max():.3f}")
        print(f"  Significant at p<0.05:      {(testable['mannwhitney_p'] < 0.05).sum()}/{len(testable)}")
        print(f"  Significant at p<0.01:      {(testable['mannwhitney_p'] < 0.01).sum()}/{len(testable)}")

        print(f"\n  Direction of differences (mutated - germline):")
        print(f"    More negative (stronger purifying in mutated): "
              f"{(testable['diff_mutated_minus_germline'] < 0).sum()}")
        print(f"    More positive (weaker purifying in mutated): "
              f"{(testable['diff_mutated_minus_germline'] > 0).sum()}")

        print(f"\n  Germline median log selection factors for testable substitutions:")
        print(f"    All below -1 (entrenched): "
              f"{(testable['median_germline'] < -1).sum()}/{len(testable)}")
        print(f"    Mutated-to also below -1: "
              f"{(testable['median_mutated'] < -1).sum()}/{len(testable)}")

    output_path = "_output/germline_vs_mutated_selection_test.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
