"""
Robustness test: re-run entrenchment analysis after excluding alleles
that are confusable (within a nucleotide distance threshold) with another
allele at an entrenched site.

If entrenchment persists after excluding at-risk alleles, then germline
misassignment cannot be driving the signal.

Replicates the within-family entrenchment logic from
v_families_entrenchment_dasm.ipynb and produces comparison plots.

Usage:
    cd dasm-epistasis-experiments

    # Confusable-pair mode (original):
    python entrenchment_exclusion_robustness_test.py

    # Exclude specific V genes (all alleles):
    python entrenchment_exclusion_robustness_test.py --exclude-genes IGHV1-69
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import combinations

from utils import load_and_process_dasm_data, sort_antibody_sites, add_cdr_shading

MODEL_NAME = "dasm_4m-v1jaffeCC+v1tangCC-joint"
DATASET_NAME = "v1rodriguez"
NUMBERING_SCHEME = "chothia"
WITHIN_FAMILIES = ["IGHV1", "IGHV3", "IGHV4"]
GERMLINE_CODONS_PATH = "germline/germline_codons_chothia.csv"
ENTRENCHED_DIR = "_output/entrenchment_analysis/chothia"
OUTPUT_DIR = "_output/exclusion_robustness"
CONFUSABLE_THRESHOLD = 10
MIN_OBS = 10
ENTRENCHMENT_THRESHOLD = -1


def compute_at_risk_alleles(threshold):
    """Find alleles at risk of misassignment at entrenched sites."""
    codons = pd.read_csv(GERMLINE_CODONS_PATH, dtype={"site": str})
    codon_wide = codons.pivot(index="v_gene", columns="site", values="codon")
    gene_to_family = (
        codons[["v_gene", "v_family"]]
        .drop_duplicates()
        .set_index("v_gene")["v_family"]
    )

    entrenched_sites = {}
    for fam in WITHIN_FAMILIES:
        path = f"{ENTRENCHED_DIR}/entrenched_aa_sites_within_{fam}.csv"
        df = pd.read_csv(path, dtype={"site": str})
        entrenched_sites[fam] = set(df["site"].unique())

    all_genes = sorted(codon_wide.index)
    at_risk = {}  # (site, family) -> set of genes

    for g1, g2 in combinations(all_genes, 2):
        r1, r2 = codon_wide.loc[g1], codon_wide.loc[g2]
        total_nt = 0
        diff_sites = []
        for site in codon_wide.columns:
            c1, c2 = r1[site], r2[site]
            if pd.isna(c1) or pd.isna(c2):
                continue
            if c1 != c2:
                d = sum(a != b for a, b in zip(c1, c2))
                total_nt += d
                diff_sites.append(site)
                if total_nt > threshold:
                    break

        if total_nt > threshold:
            continue

        f1 = gene_to_family.get(g1, "")
        f2 = gene_to_family.get(g2, "")

        for site in diff_sites:
            for fam in [f1, f2]:
                if fam in entrenched_sites and site in entrenched_sites[fam]:
                    key = (site, fam)
                    at_risk.setdefault(key, set())
                    if gene_to_family.get(g1) == fam:
                        at_risk[key].add(g1)
                    if gene_to_family.get(g2) == fam:
                        at_risk[key].add(g2)

    # Build per-family set of all at-risk alleles (union across all sites)
    at_risk_per_family = {}
    for (site, fam), genes in at_risk.items():
        at_risk_per_family.setdefault(fam, set()).update(genes)

    return at_risk, at_risk_per_family


def run_entrenchment_analysis(aa_df, v_family):
    """Replicate the within-family entrenchment logic from the notebook."""
    df = aa_df.copy()
    df = df[df["depth"] == 2]
    if "one_mutation_away" in df.columns:
        df = df[df["one_mutation_away"] == True]
    df = df[df["is_germline_codon"] == True]
    df = df.dropna(subset=["log_selection_factor"])

    # Count filter: >=10 unique PCPs per (v_family, site, parent_aa)
    counts = (
        df[["v_family", "site", "parent_aa", "pcp_index"]]
        .drop_duplicates()
        .groupby(["v_family", "site", "parent_aa"])
        .size()
        .reset_index(name="count")
    )
    counts = counts[counts["count"] >= MIN_OBS]
    df = df.merge(
        counts[["v_family", "site", "parent_aa"]],
        on=["v_family", "site", "parent_aa"],
        how="inner",
    )

    # Compute medians
    medians = (
        df.groupby(["v_family", "site", "parent_aa", "selection_factor_target_aa"])
        ["log_selection_factor"]
        .median()
        .reset_index()
    )
    medians = medians[medians["v_family"] == v_family].copy()

    # Self-merge for reciprocal pairs
    m1 = medians.rename(columns={
        "parent_aa": "aa_A",
        "selection_factor_target_aa": "aa_B",
        "log_selection_factor": "lsf_A_to_B",
    })
    m2 = medians.rename(columns={
        "parent_aa": "aa_B",
        "selection_factor_target_aa": "aa_A",
        "log_selection_factor": "lsf_B_to_A",
    })

    compare = m1.merge(m2, on=["v_family", "site", "aa_A", "aa_B"], how="inner")
    compare["sum"] = compare["lsf_A_to_B"] + compare["lsf_B_to_A"]
    compare["is_entrenched"] = (
        (compare["lsf_A_to_B"] < ENTRENCHMENT_THRESHOLD)
        & (compare["lsf_B_to_A"] < ENTRENCHMENT_THRESHOLD)
    )

    return compare


def plot_entrenchment_comparison(full_compare, excl_compare, v_family,
                                 germline_codons_df, output_dir=OUTPUT_DIR):
    """Plot entrenchment scatter: full vs excluded, side by side."""
    family_data = germline_codons_df[germline_codons_df.v_family == v_family]
    all_sites = sorted(set(family_data["site"].unique()))
    sorted_sites = sort_antibody_sites(all_sites)
    site_to_position = {site: i for i, site in enumerate(sorted_sites)}

    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True,
                             gridspec_kw={"hspace": 0.15})

    for ax, compare, title_suffix in [
        (axes[0], full_compare, "all alleles"),
        (axes[1], excl_compare, f"safe alleles only (excluded ≤{CONFUSABLE_THRESHOLD}nt confusable)"),
    ]:
        add_cdr_shading(ax, sorted_sites, numbering_scheme=NUMBERING_SCHEME)

        # Deduplicate reciprocal pairs for plotting (A→B and B→A are the same pair)
        plot_df = compare[compare["aa_A"] < compare["aa_B"]].copy()

        for is_ent in [False, True]:
            subset = plot_df[plot_df["is_entrenched"] == is_ent]
            if subset.empty:
                continue
            x = [site_to_position.get(s, -1) for s in subset["site"]]
            color = "#262626" if is_ent else "#D4D2D2"
            label = "entrenched" if is_ent else "not entrenched"
            ax.scatter(x, subset["sum"], alpha=0.85, label=label,
                       color=color, edgecolor="black", linewidth=0.5, s=50)

        n_ent = plot_df["is_entrenched"].sum()
        ax.set_title(f"{v_family} — {title_suffix} ({n_ent} entrenched pairs)",
                     fontsize=16)
        ax.set_ylabel("Selection A→B\n+\nSelection B→A", fontsize=13)
        ax.grid(True, axis="y", alpha=0.4, linewidth=0.8)
        for i in range(0, len(sorted_sites), 2):
            ax.axvline(i, color="lightgray", alpha=0.4, linewidth=0.8)
        ax.legend(bbox_to_anchor=(1.02, 0.8), loc="upper left", fontsize=13)
        ax.set_xlim(-0.5, len(sorted_sites) - 0.5)

    axes[1].set_xticks(range(len(sorted_sites)))
    axes[1].set_xticklabels(sorted_sites, rotation=90, fontsize=10)
    axes[1].set_xlabel("Site position", fontsize=13)

    plt.tight_layout()
    fig.savefig(f"{output_dir}/entrenchment_scatter_{v_family}.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_reciprocal_comparison(full_compare, excl_compare, v_family, output_dir=OUTPUT_DIR):
    """Plot reciprocal selection factors (A→B vs B→A): full vs excluded."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, compare, title_suffix in [
        (axes[0], full_compare, "all alleles"),
        (axes[1], excl_compare, "safe alleles only"),
    ]:
        plot_df = compare[compare["aa_A"] < compare["aa_B"]].copy()

        for is_ent in [False, True]:
            subset = plot_df[plot_df["is_entrenched"] == is_ent]
            if subset.empty:
                continue
            color = "#262626" if is_ent else "#D4D2D2"
            label = "entrenched" if is_ent else "not entrenched"
            ax.scatter(subset["lsf_A_to_B"], subset["lsf_B_to_A"],
                       alpha=0.85, label=label, color=color,
                       edgecolor="black", linewidth=0.5, s=50)

        ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
        ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
        ax.axhline(ENTRENCHMENT_THRESHOLD, color="red", linestyle=":",
                   alpha=0.3, linewidth=1)
        ax.axvline(ENTRENCHMENT_THRESHOLD, color="red", linestyle=":",
                   alpha=0.3, linewidth=1)
        n_ent = plot_df["is_entrenched"].sum()
        ax.set_title(f"{v_family} — {title_suffix}\n({n_ent} entrenched)",
                     fontsize=13)
        ax.set_xlabel("Selection A→B", fontsize=12)
        ax.set_ylabel("Selection B→A", fontsize=12)
        ax.set_xlim(-4.2, 2)
        ax.set_ylim(-4.2, 2)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=10)

    plt.tight_layout()
    fig.savefig(f"{output_dir}/reciprocal_scatter_{v_family}.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--exclude-genes", nargs="+", default=None,
        help="V gene names to exclude (all alleles). E.g. --exclude-genes IGHV1-69. "
             "When set, skips confusable-pair computation and excludes these genes directly."
    )
    return parser.parse_args()


def main():
    import os

    args = parse_args()
    exclude_genes = args.exclude_genes

    print("\nLoading data (this takes a minute)...")
    _, _, aa_df = load_and_process_dasm_data(
        model_name=MODEL_NAME,
        dataset_name=DATASET_NAME,
        numbering_scheme=NUMBERING_SCHEME,
    )

    germline_codons_df = pd.read_csv(GERMLINE_CODONS_PATH, dtype={"site": str})

    if exclude_genes is not None:
        # Exclude specific genes mode
        label = "genes_" + "+".join(exclude_genes)
        output_dir = f"{OUTPUT_DIR}/{label}"
        os.makedirs(output_dir, exist_ok=True)

        # Find all alleles matching the gene names
        all_alleles = aa_df["v_gene"].unique()
        excluded_alleles = {a for a in all_alleles
                            if any(a.startswith(g + "*") or a == g for g in exclude_genes)}

        print(f"\n{'=' * 90}")
        print(f"EXCLUSION ROBUSTNESS TEST (excluding genes: {', '.join(exclude_genes)})")
        print(f"  Matched alleles to exclude: {len(excluded_alleles)}")
        print(f"{'=' * 90}")

        # Build per-family exclusion sets
        at_risk_per_family = {}
        for allele in excluded_alleles:
            fam = allele.split("-")[0]  # e.g. IGHV1-69*01 -> IGHV1
            # Actually need the v_family from aa_df
            rows = aa_df[aa_df["v_gene"] == allele]
            if len(rows) > 0:
                fam = rows.iloc[0]["v_family"]
                at_risk_per_family.setdefault(fam, set()).add(allele)
    else:
        # Original confusable-pair mode
        output_dir = OUTPUT_DIR
        os.makedirs(output_dir, exist_ok=True)

        label = f"{CONFUSABLE_THRESHOLD}nt"
        print(f"Computing at-risk alleles (threshold <= {CONFUSABLE_THRESHOLD} nt)...")
        _, at_risk_per_family = compute_at_risk_alleles(CONFUSABLE_THRESHOLD)

        print(f"\n{'=' * 90}")
        print(f"EXCLUSION ROBUSTNESS TEST (threshold <= {CONFUSABLE_THRESHOLD} nt)")
        print(f"{'=' * 90}")

    all_results = []

    for v_family in WITHIN_FAMILIES:
        risk_genes = at_risk_per_family.get(v_family, set())

        # Full analysis (original)
        full_compare = run_entrenchment_analysis(aa_df, v_family)

        # Excluded analysis (remove at-risk alleles by v_gene)
        aa_df_safe = aa_df[~aa_df["v_gene"].isin(risk_genes)]
        excl_compare = run_entrenchment_analysis(aa_df_safe, v_family)

        full_entrenched = full_compare[full_compare["is_entrenched"]]
        excl_entrenched = excl_compare[excl_compare["is_entrenched"]]

        print(f"\n{'─' * 80}")
        print(f"{v_family} ({len(risk_genes)} alleles excluded)")
        print(f"{'─' * 80}")

        # Compare: which entrenched pairs survived?
        full_pairs = set(
            zip(full_entrenched["site"], full_entrenched["aa_A"], full_entrenched["aa_B"])
        )
        excl_pairs = set(
            zip(excl_entrenched["site"], excl_entrenched["aa_A"], excl_entrenched["aa_B"])
        )
        lost_pairs = full_pairs - excl_pairs
        survived_pairs = full_pairs & excl_pairs

        print(f"  Original entrenched pairs: {len(full_pairs)}")
        print(f"  After exclusion:           {len(excl_pairs)}")
        print(f"  Survived:                  {len(survived_pairs)}")
        print(f"  Lost:                      {len(lost_pairs)}")

        # Detail for each original entrenched pair
        print(f"\n  {'Site':>5} {'Pair':>8} {'Full A→B':>10} {'Full B→A':>10} {'Excl A→B':>10} {'Excl B→A':>10} {'Status':>12}")
        print(f"  {'─'*5} {'─'*8} {'─'*10} {'─'*10} {'─'*10} {'─'*10} {'─'*12}")

        for site, aa_a, aa_b in sorted(full_pairs):
            full_row = full_compare[
                (full_compare["site"] == site)
                & (full_compare["aa_A"] == aa_a)
                & (full_compare["aa_B"] == aa_b)
            ].iloc[0]

            excl_rows = excl_compare[
                (excl_compare["site"] == site)
                & (excl_compare["aa_A"] == aa_a)
                & (excl_compare["aa_B"] == aa_b)
            ]

            if len(excl_rows) > 0:
                excl_row = excl_rows.iloc[0]
                excl_ab = f"{excl_row['lsf_A_to_B']:+.2f}"
                excl_ba = f"{excl_row['lsf_B_to_A']:+.2f}"
                status = "SURVIVED" if excl_row["is_entrenched"] else "WEAKENED"
            else:
                excl_row = None
                excl_ab = "n/a"
                excl_ba = "n/a"
                status = "DROPPED"

            print(
                f"  {site:>5} {aa_a}→{aa_b:>5} "
                f"{full_row['lsf_A_to_B']:>+10.2f} {full_row['lsf_B_to_A']:>+10.2f} "
                f"{excl_ab:>10} {excl_ba:>10} {status:>12}"
            )

            all_results.append({
                "v_family": v_family,
                "site": site,
                "aa_A": aa_a,
                "aa_B": aa_b,
                "full_lsf_A_to_B": full_row["lsf_A_to_B"],
                "full_lsf_B_to_A": full_row["lsf_B_to_A"],
                "excl_lsf_A_to_B": excl_row["lsf_A_to_B"] if excl_row is not None else np.nan,
                "excl_lsf_B_to_A": excl_row["lsf_B_to_A"] if excl_row is not None else np.nan,
                "status": status,
            })

        # Generate plots
        plot_entrenchment_comparison(full_compare, excl_compare, v_family,
                                     germline_codons_df, output_dir=output_dir)
        plot_reciprocal_comparison(full_compare, excl_compare, v_family,
                                   output_dir=output_dir)
        print(f"  Plots saved to {output_dir}/")

    # Summary
    results_df = pd.DataFrame(all_results)
    print(f"\n\n{'=' * 90}")
    print("OVERALL SUMMARY")
    print(f"{'=' * 90}")
    total = len(results_df)
    survived = (results_df["status"] == "SURVIVED").sum()
    weakened = (results_df["status"] == "WEAKENED").sum()
    dropped = (results_df["status"] == "DROPPED").sum()
    print(f"Total originally entrenched pairs: {total}")
    print(f"  SURVIVED (still entrenched):     {survived} ({100*survived/total:.0f}%)")
    print(f"  WEAKENED (no longer entrenched):  {weakened} ({100*weakened/total:.0f}%)")
    print(f"  DROPPED (insufficient data):      {dropped} ({100*dropped/total:.0f}%)")

    output_path = f"{output_dir}/exclusion_robustness_test_{label}.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
