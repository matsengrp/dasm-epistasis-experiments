"""
SASA (Solvent Accessibility Surface Area) Plotting Utilities

Reusable plotting functions for solvent accessibility analysis of antibody structures.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

from utils import sort_antibody_sites, add_cdr_shading


def plot_rsa_by_vfamily(
    df,
    y_col,
    v_families=None,
    title=None,
    figsize=(20, 10),
    hue_col="is_cdr",
    save_path=None,
):
    """
    Plot RSA boxplots by site for multiple V families.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns: site, v_family_heavy, and the y_col
    y_col : str
        Column name to plot on y-axis (e.g., 'rsa_heavy_light', 'rsa_heavy_light_antigen')
    v_families : list, optional
        List of V families to plot. Defaults to ['IGHV1', 'IGHV3', 'IGHV4']
    title : str, optional
        Plot title. Auto-generated if not provided.
    figsize : tuple
        Figure size (width, height)
    hue_col : str
        Column to use for hue coloring
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig, axes : matplotlib figure and axes
    """
    if v_families is None:
        v_families = ["IGHV1", "IGHV3", "IGHV4"]

    fig, axes = plt.subplots(len(v_families), 1, figsize=figsize, sharex=True)
    if len(v_families) == 1:
        axes = [axes]

    sorted_sites = sort_antibody_sites(df["site"].unique())

    for i, gene_type in enumerate(v_families):
        filtered_data = df[df.v_family_heavy == gene_type]

        sns.boxplot(
            data=filtered_data,
            x="site",
            y=y_col,
            showfliers=False,
            ax=axes[i],
            hue=hue_col,
            order=sorted_sites,
        )

        axes[i].set_title(f"{gene_type}")
        axes[i].tick_params(axis="x", rotation=90)
        axes[i].grid()
        axes[i].legend(loc="upper left", bbox_to_anchor=(1, 1))

    plt.tight_layout()

    if title is None:
        title = f"{y_col} by Residue Number and V Family"
    plt.suptitle(title, fontsize=16, y=1.05)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=800)

    plt.show()
    return fig, axes


def plot_entrenched_rsa_comparison(
    background_data,
    highlighted_data,
    y_col,
    v_families=None,
    palette_aa=None,
    numbering_scheme="imgt",
    title=None,
    figsize=(20, 12),
    save_path=None,
):
    """
    Plot RSA comparison showing entrenched sites highlighted against background.

    Parameters
    ----------
    background_data : pd.DataFrame
        DataFrame with non-highlighted (background) data points
    highlighted_data : pd.DataFrame
        DataFrame with highlighted (entrenched) data points
    y_col : str
        Column name to plot on y-axis
    v_families : list, optional
        List of V families to plot. Defaults to ['IGHV1', 'IGHV3', 'IGHV4']
    palette_aa : dict
        Dictionary mapping amino acids to colors
    numbering_scheme : str
        Numbering scheme for CDR shading ('imgt' or 'chothia')
    title : str, optional
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig, axes : matplotlib figure and axes
    """
    if v_families is None:
        v_families = ["IGHV1", "IGHV3", "IGHV4"]

    fig, axes = plt.subplots(len(v_families), 1, figsize=figsize, sharex=True)
    if len(v_families) == 1:
        axes = [axes]

    # Compute global sorted sites for consistent x-axis
    all_sites_global = set()
    for vf in v_families:
        bg = background_data[background_data.v_family_heavy == vf]
        hl = highlighted_data[highlighted_data.v_family_heavy == vf]
        all_sites_global |= set(bg["site"].unique()) | set(hl["site"].unique())
    global_sorted_sites = sort_antibody_sites(list(all_sites_global))

    for i, gene_type in enumerate(v_families):
        filtered_bg = background_data[
            background_data.v_family_heavy == gene_type
        ].copy()
        filtered_hl = highlighted_data[
            highlighted_data.v_family_heavy == gene_type
        ].copy()

        # Apply global categories for consistent x-axis
        filtered_bg["site"] = pd.Categorical(
            filtered_bg["site"], categories=global_sorted_sites, ordered=True
        )
        filtered_hl["site"] = pd.Categorical(
            filtered_hl["site"], categories=global_sorted_sites, ordered=True
        )

        # Plot background points
        sns.scatterplot(
            data=filtered_bg,
            x="site",
            y=y_col,
            color="black",
            s=70,
            alpha=0.9,
            ax=axes[i],
            legend=False,
        )

        # Plot highlighted points on top
        sns.scatterplot(
            data=filtered_hl,
            x="site",
            y=y_col,
            hue="amino_acid",
            s=70,
            alpha=0.9,
            ax=axes[i],
            palette=palette_aa,
            legend=False,
        )

        axes[i].set_title(f"{gene_type}")
        axes[i].tick_params(axis="x", rotation=90)
        axes[i].grid()
        add_cdr_shading(axes[i], global_sorted_sites, numbering_scheme=numbering_scheme)

    # Create single shared legend
    if palette_aa:
        handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=palette_aa[aa],
                markersize=8,
            )
            for aa in sorted(palette_aa.keys())
        ]
        labels = sorted(palette_aa.keys())
        fig.legend(
            handles,
            labels,
            bbox_to_anchor=(1.02, 0.5),
            loc="center left",
            borderaxespad=0.0,
            title="Amino Acid",
            fontsize=12,
            title_fontsize=14,
        )

    plt.tight_layout()

    if title:
        plt.suptitle(title, fontsize=24, y=1.05)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=800)

    plt.show()
    return fig, axes


def plot_single_vfamily_rsa_comparison(
    v_family,
    background_data,
    highlighted_data,
    numbering_scheme,
    palette_aa,
    figsize=(20, 12),
    save_path=None,
    title=None,
):
    """
    Plot RSA comparison for a single V family with three subplots showing different RSA metrics.

    Parameters
    ----------
    v_family : str
        V family to plot (e.g., 'IGHV1', 'IGHV3', 'IGHV4')
    background_data : pd.DataFrame
        DataFrame containing background (non-highlighted) data
    highlighted_data : pd.DataFrame
        DataFrame containing highlighted (entrenched) data
    numbering_scheme : str
        Numbering scheme for antibody sites ('imgt' or 'chothia')
    palette_aa : dict
        Dictionary mapping amino acids to colors
    figsize : tuple, optional
        Figure size (width, height)
    save_path : str, optional
        Path to save the figure
    title : str, optional
        Plot title

    Returns
    -------
    fig, axes : matplotlib figure and axes objects
    """
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    fig.subplots_adjust(hspace=0.4)

    # Filter data for the specified V family
    filtered_bg = background_data[background_data.v_family_heavy == v_family].copy()
    filtered_hl = highlighted_data[highlighted_data.v_family_heavy == v_family].copy()

    # Compute unified sorted sites
    all_sites = set(filtered_bg["site"].unique()) | set(filtered_hl["site"].unique())
    sorted_sites = sort_antibody_sites(list(all_sites))
    filtered_bg["site"] = pd.Categorical(
        filtered_bg["site"], categories=sorted_sites, ordered=True
    )
    filtered_hl["site"] = pd.Categorical(
        filtered_hl["site"], categories=sorted_sites, ordered=True
    )

    # Define the three y-variables to plot
    y_vars = [
        ("rsa_heavy_light_antigen", "RSA (Heavy + Light + Antigen)"),
        ("burial_by_antigen", "Δ RSA Antigen"),
        ("burial_by_light_chain", "Δ RSA Light Chain"),
    ]

    for i, (y_var, y_label) in enumerate(y_vars):
        # Put grid behind points
        axes[i].set_axisbelow(True)

        # Plot background points (unentrenched) as hollow circles
        # Convert categorical to numeric positions for scatter
        site_positions = {site: idx for idx, site in enumerate(sorted_sites)}
        bg_x = filtered_bg["site"].map(site_positions)
        axes[i].scatter(
            bg_x,
            filtered_bg[y_var],
            facecolors="none",
            edgecolors="black",
            s=70,
            alpha=0.9,
            linewidths=1,
        )

        # Plot highlighted points (entrenched) as filled circles (no alpha)
        for aa in filtered_hl["amino_acid"].unique():
            aa_data = filtered_hl[filtered_hl["amino_acid"] == aa]
            hl_x = aa_data["site"].map(site_positions)
            axes[i].scatter(
                hl_x,
                aa_data[y_var],
                facecolors=palette_aa.get(aa, "gray"),
                edgecolors="black",
                s=180,
                alpha=1.0,
                linewidths=0.5,
            )

        # Set x-axis ticks to match site labels
        axes[i].set_xticks(range(len(sorted_sites)))
        axes[i].set_xticklabels(sorted_sites)
        axes[i].set_ylabel(y_label, fontsize=16)
        axes[i].set_title(f"{y_label}", fontsize=20)
        axes[i].tick_params(axis="x", rotation=90, labelsize=14)
        axes[i].grid(True)
        add_cdr_shading(axes[i], sorted_sites, numbering_scheme=numbering_scheme)

    axes[-1].set_xlabel("Site", fontsize=16)

    # Create shared legend with larger markers and font
    # Title "Entrenched" with amino acid entries first
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=palette_aa[aa],
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=12,
        )
        for aa in sorted(palette_aa.keys())
    ]
    labels = list(sorted(palette_aa.keys()))

    # Add blank spacer line
    handles.append(plt.Line2D([0], [0], marker="", color="w", linestyle=""))
    labels.append("")

    # Add "not entrenched" hollow circle
    handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1,
            markersize=12,
        )
    )
    labels.append("Not entrenched")

    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        borderaxespad=0.0,
        title="Entrenched",
        fontsize=12,
        title_fontsize=14,
    )

    plt.tight_layout()

    if title is None:
        title = f"RSA Analysis for {v_family}\nGermline Encoded Sites with Entrenched Amino Acids"
    fig.suptitle(title, fontsize=22, y=1.1)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=800)

    plt.show()
    return fig, axes


def plot_single_vfamily_rsa_complex_only(
    v_family,
    background_data,
    highlighted_data,
    numbering_scheme,
    palette_aa,
    figsize=(20, 5),
    save_path=None,
    title=None,
):
    """
    Plot RSA in complex only for a single V family (single subplot version).

    Parameters
    ----------
    v_family : str
        V family to plot (e.g., 'IGHV1', 'IGHV3', 'IGHV4')
    background_data : pd.DataFrame
        DataFrame containing background (non-highlighted) data
    highlighted_data : pd.DataFrame
        DataFrame containing highlighted (entrenched) data
    numbering_scheme : str
        Numbering scheme for antibody sites ('imgt' or 'chothia')
    palette_aa : dict
        Dictionary mapping amino acids to colors
    figsize : tuple, optional
        Figure size (width, height)
    save_path : str, optional
        Path to save the figure
    title : str, optional
        Plot title

    Returns
    -------
    fig, ax : matplotlib figure and axes objects
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Filter data for the specified V family
    filtered_bg = background_data[background_data.v_family_heavy == v_family].copy()
    filtered_hl = highlighted_data[highlighted_data.v_family_heavy == v_family].copy()

    # Compute unified sorted sites
    all_sites = set(filtered_bg["site"].unique()) | set(filtered_hl["site"].unique())
    sorted_sites = sort_antibody_sites(list(all_sites))
    filtered_bg["site"] = pd.Categorical(
        filtered_bg["site"], categories=sorted_sites, ordered=True
    )
    filtered_hl["site"] = pd.Categorical(
        filtered_hl["site"], categories=sorted_sites, ordered=True
    )

    # Put grid behind points
    ax.set_axisbelow(True)

    # Plot background points (unentrenched) as hollow circles
    site_positions = {site: idx for idx, site in enumerate(sorted_sites)}
    bg_x = filtered_bg["site"].map(site_positions)
    ax.scatter(
        bg_x,
        filtered_bg["rsa_heavy_light_antigen"],
        facecolors="none",
        edgecolors="black",
        s=70,
        alpha=0.9,
        linewidths=1,
    )

    # Plot highlighted points (entrenched) as filled circles (no alpha)
    for aa in filtered_hl["amino_acid"].unique():
        aa_data = filtered_hl[filtered_hl["amino_acid"] == aa]
        hl_x = aa_data["site"].map(site_positions)
        ax.scatter(
            hl_x,
            aa_data["rsa_heavy_light_antigen"],
            facecolors=palette_aa.get(aa, "gray"),
            edgecolors="black",
            s=180,
            alpha=1.0,
            linewidths=0.5,
        )

    # Set x-axis ticks to match site labels
    ax.set_xticks(range(len(sorted_sites)))
    ax.set_xticklabels(sorted_sites)
    ax.set_ylabel("RSA (Heavy + Light + Antigen)", fontsize=16)
    ax.set_xlabel("Site", fontsize=16)
    ax.tick_params(axis="x", rotation=90, labelsize=14)
    ax.grid(True)
    add_cdr_shading(ax, sorted_sites, numbering_scheme=numbering_scheme)

    # Create legend with larger markers and font
    # Title "Entrenched" with amino acid entries first
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=palette_aa[aa],
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=12,
        )
        for aa in sorted(palette_aa.keys())
    ]
    labels = list(sorted(palette_aa.keys()))

    # Add blank spacer line
    handles.append(plt.Line2D([0], [0], marker="", color="w", linestyle=""))
    labels.append("")

    # Add "not entrenched" hollow circle
    handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1,
            markersize=12,
        )
    )
    labels.append("Not entrenched")

    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        borderaxespad=0.0,
        title="Entrenched",
        fontsize=12,
        title_fontsize=14,
    )

    plt.tight_layout()

    if title is None:
        title = f"RSA in Complex for {v_family}\nGermline Encoded Sites with Entrenched Amino Acids"
    fig.suptitle(title, fontsize=22, y=1.05)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=800)

    plt.show()
    return fig, ax


def plot_rsa_summary_2x2(
    background_data,
    highlighted_data,
    v_families=None,
    x_col="phi",
    y_col="psi",
    palette_aa=None,
    title=None,
    figsize=(8, 8),
    save_path=None,
):
    """
    Plot 2x2 grid showing scatter plots (e.g., Ramachandran) for multiple V families.

    Parameters
    ----------
    background_data : pd.DataFrame
        Background data points
    highlighted_data : pd.DataFrame
        Highlighted (entrenched) data points
    v_families : list, optional
        V families to plot. Defaults to ['IGHV1', 'IGHV3', 'IGHV4']
    x_col, y_col : str
        Column names for x and y axes
    palette_aa : dict
        Amino acid color palette
    title : str, optional
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    fig, axes : matplotlib figure and axes
    """
    if v_families is None:
        v_families = ["IGHV1", "IGHV3", "IGHV4"]

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    for i, gene_type in enumerate(v_families):
        if i >= len(axes):
            break

        filtered_bg = background_data[
            background_data.v_family_heavy == gene_type
        ].copy()
        filtered_hl = highlighted_data[
            highlighted_data.v_family_heavy == gene_type
        ].copy()

        # Compute sorted sites
        all_sites = set(filtered_bg["site"].unique()) | set(
            filtered_hl["site"].unique()
        )
        vfamily_sorted_sites = sort_antibody_sites(list(all_sites))
        filtered_bg["site"] = pd.Categorical(
            filtered_bg["site"], categories=vfamily_sorted_sites, ordered=True
        )
        filtered_hl["site"] = pd.Categorical(
            filtered_hl["site"], categories=vfamily_sorted_sites, ordered=True
        )

        # Plot background
        sns.scatterplot(
            data=filtered_bg,
            x=x_col,
            y=y_col,
            color="black",
            s=70,
            alpha=0.9,
            ax=axes[i],
        )

        # Plot highlighted
        sns.scatterplot(
            data=filtered_hl,
            x=x_col,
            y=y_col,
            hue="amino_acid",
            s=70,
            alpha=0.9,
            ax=axes[i],
            palette=palette_aa,
            legend=False,
        )

        axes[i].set_title(f"{gene_type}")
        if x_col in ("phi", "psi"):
            axes[i].set_xlim(-180, 180)
            axes[i].set_ylim(-180, 180)
        axes[i].grid()

    # Hide unused subplot
    if len(v_families) < 4:
        axes[3].axis("off")

    # Create shared legend
    if palette_aa:
        handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=palette_aa[aa],
                markersize=8,
            )
            for aa in sorted(palette_aa.keys())
        ]
        labels = sorted(palette_aa.keys())
        fig.legend(
            handles,
            labels,
            bbox_to_anchor=(1.02, 0.5),
            loc="center left",
            borderaxespad=0.0,
            title="Amino Acid",
        )

    plt.tight_layout()

    if title:
        plt.suptitle(title, fontsize=20, y=1.1)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=800)

    plt.show()
    return fig, axes


def examine_site(
    df,
    vfamily,
    site,
    numbering_scheme,
    palette_aa,
    comparison_type="within",
    save_fig=False,
    output_dir="figures/site_analysis",
):
    """
    Create detailed examination plots for a specific site in a V family.

    Parameters
    ----------
    df : pd.DataFrame
        Main data frame with RSA and structural data
    vfamily : str
        V family to examine (e.g., 'IGHV3')
    site : str
        Site identifier
    numbering_scheme : str
        Numbering scheme ('imgt' or 'chothia')
    palette_aa : dict
        Amino acid color palette
    comparison_type : str
        Type of comparison ('within' for within-family)
    save_fig : bool
        Whether to save the figure
    output_dir : str
        Directory for saving figures

    Returns
    -------
    fig : matplotlib figure
    """
    cur_df = df[
        (df.is_germline == True) & (df.site == site) & (df.v_family_heavy == vfamily)
    ]

    # Use GridSpec for layout control
    fig = plt.figure(figsize=(7.5, 7))
    gs = fig.add_gridspec(2, 14, wspace=0.25, hspace=0.4)

    # Top row: 2 plots
    ax1 = fig.add_subplot(gs[0, 0:6])  # Degree of entrenchment (left)
    ax0 = fig.add_subplot(gs[0, 8:14])  # Germline AA identities (right)

    # Bottom row: 3 plots
    ax3 = fig.add_subplot(gs[1, 0:4])  # RSA in complex
    ax4 = fig.add_subplot(gs[1, 5:9])  # Burial by antigen
    ax5 = fig.add_subplot(gs[1, 10:14])  # Burial by light chain

    ax = [ax0, ax1, None, ax3, ax4, ax5]

    # Font sizes
    TITLE_SIZE = 11
    LABEL_SIZE = 10
    TICK_SIZE = 9
    ANNOT_SIZE = 8

    # Top-right: Germline amino acid identities
    sns.countplot(
        cur_df, x="amino_acid", hue="amino_acid", palette=palette_aa, ax=ax[0]
    )
    ax[0].set_title("Germline amino acid identities", fontsize=TITLE_SIZE)
    ax[0].set_xlabel("Amino acid", fontsize=LABEL_SIZE)
    ax[0].set_ylabel("Count", fontsize=LABEL_SIZE)
    ax[0].tick_params(labelsize=TICK_SIZE)

    # Bottom-left: RSA in complex
    sns.boxplot(
        cur_df,
        x="amino_acid",
        y="rsa_heavy_light_antigen",
        showfliers=False,
        hue="amino_acid",
        palette=palette_aa,
        ax=ax[3],
        whis=[5, 95],
    )
    ax[3].set_title("RSA in complex", fontsize=TITLE_SIZE)
    ax[3].set_xlabel("Amino acid", fontsize=LABEL_SIZE)
    ax[3].set_ylabel("")
    ax[3].tick_params(labelsize=TICK_SIZE)
    ax[3].grid(True, axis="y", alpha=0.3, linestyle="-")

    # Bottom-middle: Δ RSA Antigen (now positive = more buried)
    sns.boxplot(
        cur_df,
        x="amino_acid",
        y="burial_by_antigen",
        showfliers=False,
        hue="amino_acid",
        palette=palette_aa,
        ax=ax[4],
        whis=[5, 95],
    )
    ax[4].set_title("Δ RSA Antigen", fontsize=TITLE_SIZE)
    ax[4].set_xlabel("Amino acid", fontsize=LABEL_SIZE)
    ax[4].set_ylabel("")
    ax[4].tick_params(labelsize=TICK_SIZE)
    ax[4].tick_params(labelleft=False)
    ax[4].grid(True, axis="y", alpha=0.3, linestyle="-")

    # Bottom-right: Δ RSA Light Chain (now positive = more buried)
    sns.boxplot(
        cur_df,
        x="amino_acid",
        y="burial_by_light_chain",
        showfliers=False,
        hue="amino_acid",
        palette=palette_aa,
        ax=ax[5],
        whis=[5, 95],
    )
    ax[5].set_title("Δ RSA Light Chain", fontsize=TITLE_SIZE)
    ax[5].set_xlabel("Amino acid", fontsize=LABEL_SIZE)
    ax[5].set_ylabel("")
    ax[5].tick_params(labelsize=TICK_SIZE)
    ax[5].tick_params(labelleft=False)
    ax[5].grid(True, axis="y", alpha=0.3, linestyle="-")

    # Set shared y-axis limits for all three bottom panels
    all_ylims = [ax[3].get_ylim(), ax[4].get_ylim(), ax[5].get_ylim()]
    ymin = min(lim[0] for lim in all_ylims)
    ymax = max(lim[1] for lim in all_ylims)
    ax[3].set_ylim(ymin, ymax)
    ax[4].set_ylim(ymin, ymax)
    ax[5].set_ylim(ymin, ymax)

    if comparison_type == "within":
        # Load entrenchment results
        entrenchment_results_df = pd.read_csv(
            f"_output/entrenchment_analysis/{numbering_scheme}/comparison_within_{vfamily}.csv"
        )
        entrenchment_results_df["substitution"] = (
            entrenchment_results_df["parent_aa_1_and_target_aa_2"]
            + "→"
            + entrenchment_results_df["parent_aa_2_and_target_aa_1"]
        )

        site_data = entrenchment_results_df[entrenchment_results_df.site == site]
        # Remove duplicated pairs
        site_data = site_data[
            site_data["parent_aa_1_and_target_aa_2"]
            <= site_data["parent_aa_2_and_target_aa_1"]
        ]

        # Top-left: Degree of entrenchment
        # Color scheme matches compare_group_2_sites_between_vfamilies_with_within_overlay
        sns.scatterplot(
            data=site_data,
            x="log_selection_factor_1",
            y="log_selection_factor_2",
            hue="is_entrenched",
            palette={True: "#262626", False: "#D4D2D2"},
            ax=ax[1],
            legend=False,
            s=40,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.5,
        )
        ax[1].set_xlim(-4, 1)
        ax[1].set_ylim(-4, 1)
        ax[1].axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.3)
        ax[1].axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.3)
        ax[1].axline((0, 0), slope=1, linestyle="--", color="red", alpha=0.3)
        ax[1].set_title("Degree of entrenchment", fontsize=TITLE_SIZE)
        ax[1].set_xlabel("Forward selection (A→B)", fontsize=LABEL_SIZE)
        ax[1].set_ylabel("Backward selection (A←B)", fontsize=LABEL_SIZE)
        ax[1].tick_params(labelsize=TICK_SIZE)

        # Annotate points
        x_threshold = 0.5
        y_threshold = 0.3
        x_offset = 0.18

        points = list(
            zip(
                site_data["log_selection_factor_1"],
                site_data["log_selection_factor_2"],
                site_data["substitution"],
            )
        )

        # Ad-hoc overrides: specific labels to pull into a vertical list with lines
        callout_overrides = {
            ("33", "IGHV3"): {
                "labels": {"A→G", "D→E", "A→T", "G→S"},
                "label_x": -3.7,
                "label_y_start": 0.5,
                "label_spacing": 0.3,
            },
            ("50", "IGHV3"): {
                "labels": {"F→Y", "L→V", "A→V"},
                "label_x": -1.5,
                "label_y_start": 0.7,
                "label_spacing": 0.3,
            },
        }

        override = callout_overrides.get((site, vfamily))
        callout_labels = override["labels"] if override else set()

        for x, y, label in points:
            if label in callout_labels:
                continue  # handled below

            has_point_to_right = any(
                (ox > x) and (ox - x < x_threshold) and (abs(oy - y) < y_threshold)
                for ox, oy, _ in points
                if (ox, oy) != (x, y)
            )
            if has_point_to_right:
                ax[1].text(
                    x - x_offset,
                    y,
                    label,
                    fontsize=ANNOT_SIZE,
                    ha="right",
                    va="center",
                    zorder=5,
                )
            else:
                ax[1].text(
                    x + x_offset,
                    y,
                    label,
                    fontsize=ANNOT_SIZE,
                    ha="left",
                    va="center",
                    zorder=5,
                )

        # Draw callout list for overridden labels
        if override:
            callout_pts = [(x, y, l) for x, y, l in points if l in callout_labels]
            callout_pts.sort(key=lambda p: -p[1])  # top-to-bottom
            lx = override["label_x"]
            ly_start = override["label_y_start"]
            spacing = override["label_spacing"]

            for i, (x, y, label) in enumerate(callout_pts):
                ly = ly_start - i * spacing
                ax[1].annotate(
                    label,
                    xy=(x, y),
                    xytext=(lx, ly),
                    fontsize=ANNOT_SIZE,
                    ha="left",
                    va="center",
                    arrowprops=dict(
                        arrowstyle="-", color="black", alpha=0.4, linewidth=0.8
                    ),
                    zorder=5,
                )

    fig.suptitle(f"Site {site} in {vfamily}", fontsize=14, y=0.98)

    if save_fig:
        import os

        os.makedirs(output_dir, exist_ok=True)
        fig.savefig(
            f"{output_dir}/site_{site}_{vfamily}_{numbering_scheme}_rsa_analysis.png",
            bbox_inches="tight",
            dpi=800,
        )

    return fig


def plot_ramachandran_by_site_range(
    df,
    center_site,
    palette_aa=None,
    v_families=["IGHV1", "IGHV3"],
    site_range=5,
    save_fig=False,
    output_dir="figures",
    germline_filter=None,
):
    """
    Create Ramachandran plots for a range of sites centered around center_site.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing angle data with columns: v_family_heavy, site, amino_acid, phi, psi, is_germline, pdb_id
    center_site : int or str
        Center site number
    palette_aa : dict, optional
        Color palette for amino acids
    v_families : list, optional
        V families to plot. Defaults to ['IGHV1', 'IGHV3']
    site_range : int
        Number of sites to plot centered around center_site
    save_fig : bool
        Whether to save the figure
    output_dir : str
        Directory for saving figures
    germline_filter : list of tuples, optional
        Filter sequences by germline amino acid at a specific site.
        Each tuple is (v_family, site, amino_acid), e.g.:
        [('IGHV1', '9', 'A'), ('IGHV3', '9', 'G')]
        Only pdb_ids that match at least one filter condition (having the
        specified germline amino acid at the given site for that v_family)
        are included. Plots still only show germline residues.

    Returns
    -------
    fig, axes : matplotlib figure and axes
    """

    # Apply germline_filter: keep only pdb_ids matching at least one condition
    if germline_filter is not None:
        matching_pdb_ids = set()
        for v_family, filter_site, filter_aa in germline_filter:
            mask = (
                (df.v_family_heavy == v_family)
                & (df.site == str(filter_site))
                & (df.amino_acid == filter_aa)
                & (df.is_germline == True)
            )
            matching_pdb_ids.update(df.loc[mask, "pdb_id"].unique())
        df = df[df.pdb_id.isin(matching_pdb_ids)]

    # Prepare site range
    center_site_int = int(center_site)
    sites_to_plot = [
        str(i)
        for i in range(
            center_site_int - site_range // 2, center_site_int + site_range // 2 + 1
        )
    ]
    sites_to_plot = sort_antibody_sites(sites_to_plot)

    # Create figure
    n_rows = len(v_families)
    n_cols = len(sites_to_plot)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 3 * n_rows))

    # Handle array shape
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    # Track legend items
    all_amino_acids = []
    legend_handles = []
    legend_labels = []

    for row_idx, v_family in enumerate(v_families):
        for col_idx, site in enumerate(sites_to_plot):
            ax = axes[row_idx, col_idx]

            data = df[
                (df.v_family_heavy == v_family)
                & (df.site == site)
                & (df.is_germline == True)
            ].copy()

            if len(data) > 0:
                if palette_aa is not None:
                    for aa in sorted(data.amino_acid.unique()):
                        aa_data = data[data.amino_acid == aa]
                        color = palette_aa.get(aa, "gray")
                        scatter = ax.scatter(
                            aa_data.phi,
                            aa_data.psi,
                            alpha=0.4,
                            s=50,
                            color=color,
                            edgecolor="black",
                            linewidth=0.5,
                        )

                        if aa not in all_amino_acids:
                            all_amino_acids.append(aa)
                            legend_handles.append(scatter)
                            legend_labels.append(aa)
                else:
                    sns.scatterplot(
                        data=data,
                        x="phi",
                        y="psi",
                        hue="amino_acid",
                        ax=ax,
                        alpha=0.4,
                        s=50,
                        edgecolor="black",
                        linewidth=0.5,
                        legend=False,
                    )

                ax.axhline(
                    y=0, color="lightgray", linestyle="--", linewidth=0.5, alpha=0.3
                )
                ax.axvline(
                    x=0, color="lightgray", linestyle="--", linewidth=0.5, alpha=0.3
                )
                ax.set_xlim(-180, 180)
                ax.set_ylim(-180, 180)
                ax.set_aspect("equal")
                ax.set_xlabel("Phi (φ)", fontsize=10)
                ax.set_ylabel("Psi (ψ)", fontsize=10)
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_xlim(-180, 180)
                ax.set_ylim(-180, 180)

            if row_idx == 0:
                ax.set_title(f"Site {site}", fontsize=11, fontweight="bold")
            if col_idx == 0:
                ax.set_ylabel(f"{v_family}\nPsi (ψ)", fontsize=10, fontweight="bold")

    plt.tight_layout()

    if legend_handles:
        leg = fig.legend(
            legend_handles,
            legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.05),
            ncol=min(len(legend_labels), 10),
            fontsize=9,
            frameon=True,
            title="Amino Acid",
        )
        for handle in leg.legend_handles:
            handle.set_alpha(1.0)
        plt.subplots_adjust(bottom=0.1)

    if save_fig:
        import os

        os.makedirs(output_dir, exist_ok=True)
        fig.savefig(
            f"{output_dir}/site_{center_site}_angle_analysis.png",
            bbox_inches="tight",
            dpi=300,
        )

    return fig, axes
