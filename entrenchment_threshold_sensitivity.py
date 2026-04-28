"""Threshold sensitivity analysis for reciprocal-entrenchment calls.

Reads the precomputed comparison_*.csv files produced by
`compare_group_2_sites_within_vfamily_with_distribution` and
`compare_group_2_sites_between_vfamilies_with_distribution`, then re-classifies
each reciprocal pair as entrenched under a grid of log-selection-factor
thresholds. For every comparison CSV we emit a multi-panel figure that mirrors
the published plot style (same x axis, same CDR shading, same germline amino
acid distribution panel at the bottom) but stacks one scatter row per threshold
so a reviewer can see how the called set shifts around the published -1 cutoff.

No model re-runs are required: all scatter inputs live in
`_output/entrenchment_analysis/chothia/comparison_*.csv`, and the distribution
panel is rebuilt directly from `germline/germline_codons_chothia.csv`.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from netam.sequences import AA_STR_SORTED

from utils import add_cdr_shading, sort_antibody_sites

REPO_ROOT = Path(__file__).resolve().parent
INPUT_DIR = REPO_ROOT / "_output" / "entrenchment_analysis" / "chothia"
GERMLINE_CSV = REPO_ROOT / "germline" / "germline_codons_chothia.csv"
OUTPUT_DIR = INPUT_DIR  # save alongside existing outputs

# Grid requested by the reviewer response: brackets the published -1 cutoff
# symmetrically in 0.2 steps.
THRESHOLDS = [-0.6, -0.8, -1.0, -1.2, -1.4]
PUBLISHED_THRESHOLD = -1.0

# Fine-grained grid used for the count-vs-threshold line plot.
SWEEP_THRESHOLDS = [round(-2.0 + i * 0.1, 2) for i in range(21)]

NUMBERING_SCHEME = "chothia"
CHAIN = "heavy"

WITHIN_FAMILIES = ["IGHV1", "IGHV3", "IGHV4"]
BETWEEN_PAIRS = [("IGHV1", "IGHV3"), ("IGHV1", "IGHV4"), ("IGHV3", "IGHV4")]

# Match the notebook's palette definition exactly.
PALETTE_AA = dict(zip(AA_STR_SORTED, sns.color_palette("tab20", len(AA_STR_SORTED))))


def dedupe_within_family(df: pd.DataFrame) -> pd.DataFrame:
    """Within a single V family each reciprocal pair appears twice (A->B and
    B->A). Keep a single canonical copy, matching the filter used in
    ``compare_group_2_sites_within_vfamily_with_distribution``."""
    mask = df["parent_aa_1_and_target_aa_2"] < df["parent_aa_2_and_target_aa_1"]
    return df[mask].copy()


def draw_scatter_panel(
    ax,
    compare_df: pd.DataFrame,
    sorted_sites: list,
    site_to_position: dict,
    threshold: float,
    total_pairs: int,
) -> None:
    """Render a single sum-vs-site scatter row for the given threshold."""
    add_cdr_shading(
        ax, sorted_sites, numbering_scheme=NUMBERING_SCHEME, chain=CHAIN
    )

    is_entrenched = (
        (compare_df["log_selection_factor_1"] < threshold)
        & (compare_df["log_selection_factor_2"] < threshold)
    )

    for flag, color, label in [
        (False, "#D4D2D2", "not entrenched"),
        (True, "#262626", "entrenched"),
    ]:
        subset = compare_df[is_entrenched == flag]
        if subset.empty:
            continue
        xs = [site_to_position[s] for s in subset["site"]]
        ax.scatter(
            xs,
            subset["sum"],
            alpha=0.85,
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.5,
        )

    ax.set_ylabel("Sum of log\nselection factors", fontsize=10)
    ax.grid(True, alpha=0.4, linewidth=0.8)
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)

    ax.set_title(
        f"threshold = {threshold:g}",
        fontsize=10,
        loc="left",
    )

    if threshold == PUBLISHED_THRESHOLD:
        for spine in ax.spines.values():
            spine.set_edgecolor("#b22222")
            spine.set_linewidth(1.6)


def draw_distribution_panel_within(
    ax, germline_codons_df: pd.DataFrame, v_family: str,
    sorted_sites: list, site_to_position: dict,
) -> None:
    """Stacked bar panel of germline AA composition for a single V family.

    Mirrors the bottom panel of
    ``compare_group_2_sites_within_vfamily_with_distribution``."""
    add_cdr_shading(
        ax, sorted_sites, numbering_scheme=NUMBERING_SCHEME, chain=CHAIN
    )

    family_data = germline_codons_df[germline_codons_df["v_family"] == v_family]
    site_total = family_data["v_gene"].nunique()

    family_site_counts = {
        aa: family_data.loc[family_data["amino_acid"] == aa, "site"].value_counts()
        for aa in AA_STR_SORTED
    }

    for site in sorted_sites:
        site_pos = site_to_position[site]
        bottom = 0.0
        for aa in AA_STR_SORTED:
            count = family_site_counts[aa].get(site, 0)
            if count > 0 and site_total > 0:
                percentage = (count / site_total) * 100
                ax.bar(
                    site_pos,
                    percentage,
                    bottom=bottom,
                    color=PALETTE_AA[aa],
                    width=0.8,
                    edgecolor="black",
                    linewidth=0.5,
                    align="center",
                )
                bottom += percentage

    ax.set_ylabel("Germline AA\ndistribution", fontsize=10)
    ax.axhline(0, color="black", linestyle="-", alpha=0.8, linewidth=1.0, zorder=0)
    ax.set_yticks([])


def draw_distribution_panel_between(
    ax, germline_codons_df: pd.DataFrame, v_family1: str, v_family2: str,
    sorted_sites: list, site_to_position: dict,
) -> None:
    """Mirrored stacked bar panel for a pair of V families.

    Mirrors the bottom panel of
    ``compare_group_2_sites_between_vfamilies_with_distribution``."""
    add_cdr_shading(
        ax, sorted_sites, numbering_scheme=NUMBERING_SCHEME, chain=CHAIN
    )

    family1_data = germline_codons_df[germline_codons_df["v_family"] == v_family1]
    family2_data = germline_codons_df[germline_codons_df["v_family"] == v_family2]
    site_total1 = family1_data["v_gene"].nunique()
    site_total2 = family2_data["v_gene"].nunique()

    family1_site_counts = {
        aa: family1_data.loc[family1_data["amino_acid"] == aa, "site"].value_counts()
        for aa in AA_STR_SORTED
    }
    family2_site_counts = {
        aa: family2_data.loc[family2_data["amino_acid"] == aa, "site"].value_counts()
        for aa in AA_STR_SORTED
    }

    for site in sorted_sites:
        site_pos = site_to_position[site]
        bottom1 = 0.0
        for aa in AA_STR_SORTED:
            count = family1_site_counts[aa].get(site, 0)
            if count > 0 and site_total1 > 0:
                percentage = (count / site_total1) * 100
                ax.bar(
                    site_pos,
                    percentage,
                    bottom=bottom1,
                    color=PALETTE_AA[aa],
                    width=0.8,
                    edgecolor="black",
                    linewidth=0.5,
                    align="center",
                )
                bottom1 += percentage

        bottom2 = 0.0
        for aa in AA_STR_SORTED:
            count = family2_site_counts[aa].get(site, 0)
            if count > 0 and site_total2 > 0:
                percentage = -(count / site_total2) * 100
                ax.bar(
                    site_pos,
                    percentage,
                    bottom=bottom2,
                    color=PALETTE_AA[aa],
                    width=0.8,
                    edgecolor="black",
                    linewidth=0.5,
                    align="center",
                )
                bottom2 += percentage

    ax.set_ylabel("Germline AA\ndistribution", fontsize=10)
    ax.axhline(0, color="black", linestyle="-", alpha=0.8, linewidth=1.5, zorder=0)
    ax.text(
        -0.01, 0.75, v_family1, transform=ax.transAxes, fontsize=10,
        ha="right", va="center", rotation=90,
    )
    ax.text(
        -0.01, 0.25, v_family2, transform=ax.transAxes, fontsize=10,
        ha="right", va="center", rotation=90,
    )
    ax.set_yticks([])


def build_count_sweep_figure(
    compare_df: pd.DataFrame,
    title: str,
    output_stem: Path,
) -> None:
    """Single-panel line plot: number of entrenched pairs vs threshold.

    Sweeps ``SWEEP_THRESHOLDS`` and counts how many reciprocal pairs satisfy
    the entrenchment criterion at each value, with the published ``-1`` cutoff
    marked by a red dashed vertical line. Input ``compare_df`` should already
    be deduplicated to one row per reciprocal pair so the counts agree with
    the subtitles in the multi-panel ``threshold_sensitivity_*`` figures.
    """
    counts = []
    for threshold in SWEEP_THRESHOLDS:
        is_entrenched = (
            (compare_df["log_selection_factor_1"] < threshold)
            & (compare_df["log_selection_factor_2"] < threshold)
        )
        counts.append(int(is_entrenched.sum()))

    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.0))
    ax.plot(
        SWEEP_THRESHOLDS,
        counts,
        marker="o",
        linewidth=1.5,
        markersize=4,
        color="#1f77b4",
    )
    ax.axvline(
        x=PUBLISHED_THRESHOLD,
        color="#b22222",
        linestyle="--",
        linewidth=1.2,
        alpha=0.8,
    )
    ax.set_xlabel("Threshold (log selection factor)", fontsize=9)
    ax.set_ylabel("Number of entrenched pairs", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    pdf_path = output_stem.with_suffix(".pdf")
    png_path = output_stem.with_suffix(".png")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {pdf_path.relative_to(REPO_ROOT)}")


def build_sensitivity_figure(
    compare_df: pd.DataFrame,
    germline_codons_df: pd.DataFrame,
    sorted_sites: list,
    title: str,
    output_stem: Path,
    distribution_kind: str,
    family1: str,
    family2: str | None = None,
) -> None:
    """Assemble the stacked scatter + distribution figure for one comparison."""
    site_to_position = {site: i for i, site in enumerate(sorted_sites)}
    total_pairs = len(compare_df)
    n_thresholds = len(THRESHOLDS)

    # One scatter row per threshold plus one distribution row at the bottom.
    height_ratios = [1.0] * n_thresholds + [1.3]
    fig_height = 2.6 * n_thresholds + 3.0
    fig, axes = plt.subplots(
        n_thresholds + 1,
        1,
        figsize=(20, fig_height),
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.30},
    )

    scatter_axes = axes[:-1]
    dist_ax = axes[-1]

    for ax, threshold in zip(scatter_axes, THRESHOLDS):
        draw_scatter_panel(
            ax, compare_df, sorted_sites, site_to_position, threshold, total_pairs
        )

    if distribution_kind == "within":
        draw_distribution_panel_within(
            dist_ax, germline_codons_df, family1, sorted_sites, site_to_position
        )
    elif distribution_kind == "between":
        assert family2 is not None
        draw_distribution_panel_between(
            dist_ax, germline_codons_df, family1, family2,
            sorted_sites, site_to_position,
        )
    else:
        raise ValueError(distribution_kind)

    # X-axis settings mirror the notebook: tick at every site, labels on the
    # bottom panel only.
    for ax in scatter_axes:
        ax.set_xlim(-0.5, len(sorted_sites) - 0.5)
        ax.set_xticks(range(0, len(sorted_sites), 2))
        ax.set_xticklabels([])
    dist_ax.set_xlim(-0.5, len(sorted_sites) - 0.5)
    dist_ax.set_xticks(range(len(sorted_sites)))
    dist_ax.set_xticklabels(sorted_sites, rotation=90, fontsize=8)
    dist_ax.set_xlabel("Site position (Chothia numbering)", fontsize=11)

    # Shared entrenched/not-entrenched legend on the top scatter row.
    handles, labels = scatter_axes[0].get_legend_handles_labels()
    order = {lbl: h for h, lbl in zip(handles, labels)}
    ordered_labels = [lbl for lbl in ["entrenched", "not entrenched"] if lbl in order]
    if ordered_labels:
        scatter_axes[0].legend(
            [order[lbl] for lbl in ordered_labels],
            ordered_labels,
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
            frameon=False,
        )

    # Amino acid legend on the distribution row, matching the notebook figure.
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=PALETTE_AA[aa],
                      edgecolor="black", linewidth=0.5)
        for aa in AA_STR_SORTED
    ]
    dist_ax.legend(
        legend_handles,
        list(AA_STR_SORTED),
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        title="Amino Acids",
        fontsize=9,
        ncol=2,
    )

    fig.suptitle(title, fontsize=14, y=0.995)
    # The distribution panel uses axes-transform text, which is incompatible
    # with tight_layout, so we set padding manually.
    fig.subplots_adjust(
        left=0.06, right=0.90, top=0.965, bottom=0.07, hspace=0.30,
    )

    pdf_path = output_stem.with_suffix(".pdf")
    png_path = output_stem.with_suffix(".png")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {pdf_path.relative_to(REPO_ROOT)}")


def run_within(v_family: str, germline_codons_df: pd.DataFrame) -> None:
    csv_path = INPUT_DIR / f"comparison_within_{v_family}.csv"
    compare_df = pd.read_csv(csv_path, dtype={"site": str})
    compare_df = dedupe_within_family(compare_df)

    family_sites = (
        germline_codons_df.loc[
            germline_codons_df["v_family"] == v_family, "site"
        ]
        .unique()
        .tolist()
    )
    sorted_sites = sort_antibody_sites(family_sites)

    title = f"Threshold sensitivity of entrenchment calls within {v_family}"
    output_stem = OUTPUT_DIR / f"threshold_sensitivity_within_{v_family}"
    build_sensitivity_figure(
        compare_df,
        germline_codons_df,
        sorted_sites,
        title,
        output_stem,
        distribution_kind="within",
        family1=v_family,
    )

    sweep_title = f"Entrenched pairs vs threshold within {v_family}"
    sweep_stem = OUTPUT_DIR / f"threshold_sensitivity_sweep_within_{v_family}"
    build_count_sweep_figure(compare_df, sweep_title, sweep_stem)


def run_between(
    v_family1: str, v_family2: str, germline_codons_df: pd.DataFrame
) -> None:
    csv_path = INPUT_DIR / f"comparison_{v_family1}_vs_{v_family2}.csv"
    compare_df = pd.read_csv(csv_path, dtype={"site": str})

    family_sites = set(
        germline_codons_df.loc[
            germline_codons_df["v_family"] == v_family1, "site"
        ].unique()
    ) | set(
        germline_codons_df.loc[
            germline_codons_df["v_family"] == v_family2, "site"
        ].unique()
    )
    sorted_sites = sort_antibody_sites(family_sites)

    title = (
        f"Threshold sensitivity of entrenchment calls: {v_family1} vs {v_family2}"
    )
    output_stem = (
        OUTPUT_DIR / f"threshold_sensitivity_{v_family1}_vs_{v_family2}"
    )
    build_sensitivity_figure(
        compare_df,
        germline_codons_df,
        sorted_sites,
        title,
        output_stem,
        distribution_kind="between",
        family1=v_family1,
        family2=v_family2,
    )

    sweep_title = (
        f"Entrenched pairs vs threshold: {v_family1} vs {v_family2}"
    )
    sweep_stem = (
        OUTPUT_DIR / f"threshold_sensitivity_sweep_{v_family1}_vs_{v_family2}"
    )
    build_count_sweep_figure(compare_df, sweep_title, sweep_stem)


def main() -> None:
    germline_codons_df = pd.read_csv(GERMLINE_CSV, dtype={"site": str})

    for family in WITHIN_FAMILIES:
        run_within(family, germline_codons_df)
    for f1, f2 in BETWEEN_PAIRS:
        run_between(f1, f2, germline_codons_df)


if __name__ == "__main__":
    main()
