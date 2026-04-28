#!/usr/bin/env python
"""Regression analysis: fit model to translate mutation frequency to DASM branch length."""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from dnsmex.dxsm_data import train_val_df_of_multiname
from dnsmex.local import localify


def mutation_frequency(parent_seq, child_seq):
    """Compute fraction of nucleotide sites that differ between parent and child."""
    assert len(parent_seq) == len(child_seq)
    mismatches = sum(1 for p, c in zip(parent_seq, child_seq) if p != c)
    return mismatches / len(parent_seq)


def load_data():
    """Load and prepare the data."""
    print("Loading PCP data...")
    pcp_df = train_val_df_of_multiname("v1jaffeCC+v1tangCC")

    pcp_df["mut_freq"] = [
        mutation_frequency(p, c)
        for p, c in zip(pcp_df["parent_heavy"], pcp_df["child_heavy"])
    ]

    model_dir = localify("DASM_TRAINED_MODELS_DIR")
    model_name = "dasm_4m-v1jaffeCC+v1tangCC-joint"

    train_bls = pd.read_csv(f"{model_dir}/{model_name}.train_branch_lengths.csv")
    val_bls = pd.read_csv(f"{model_dir}/{model_name}.val_branch_lengths.csv")

    train_df = pcp_df[pcp_df["in_train"]].reset_index(drop=True).copy()
    val_df = pcp_df[~pcp_df["in_train"]].reset_index(drop=True).copy()

    train_df["dasm_branch_length"] = train_bls["branch_length"].values
    val_df["dasm_branch_length"] = val_bls["branch_length"].values

    return train_df, val_df


def fit_linear(X, y, name):
    """Fit linear regression and report stats."""
    result = stats.linregress(X, y)
    y_pred = result.slope * X + result.intercept

    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    rmse = np.sqrt(np.mean((y - y_pred) ** 2))

    print(f"\n{name}")
    print(f"  y = {result.slope:.4f} * x + {result.intercept:.4f}")
    print(f"  R² = {r2:.4f}")
    print(f"  RMSE = {rmse:.6f}")

    return result, y_pred


def main():
    train_df, val_df = load_data()

    # Use training data for fitting
    df = train_df

    print("\n" + "=" * 60)
    print("LINEAR REGRESSION: Mutation Frequency → DASM Branch Length")
    print("=" * 60)

    model_mut, pred_mut = fit_linear(
        df["mut_freq"], df["dasm_branch_length"],
        "Mutation Frequency → DASM Branch Length"
    )

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    sort_idx = np.argsort(df["mut_freq"].values)

    ss_res = np.sum((df["dasm_branch_length"] - pred_mut) ** 2)
    ss_tot = np.sum((df["dasm_branch_length"] - df["dasm_branch_length"].mean()) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    # Trim to 99th percentile for clearer view
    xlim = df["mut_freq"].quantile(0.99) * 1.1
    ylim = df["dasm_branch_length"].quantile(0.99) * 1.1

    ax.hexbin(
        df["mut_freq"], df["dasm_branch_length"],
        gridsize=50, cmap='Blues', mincnt=1, alpha=0.7,
        extent=[0, xlim, 0, ylim]
    )
    # Plot regression line within view
    x_line = np.array([0, xlim])
    y_line = model_mut.slope * x_line + model_mut.intercept
    ax.plot(x_line, y_line, 'r-', linewidth=2,
            label=f'y = {model_mut.slope:.2f}x {model_mut.intercept:+.4f}\nR² = {r2:.4f}')
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    ax.set_xlabel("Mutation Frequency")
    ax.set_ylabel("DASM Branch Length")
    ax.set_title("Mutation Frequency → DASM Branch Length (99th percentile view)")
    ax.legend(loc='upper left')
    plt.colorbar(ax.collections[0], ax=ax, label='Count')

    plt.tight_layout()
    output_path = os.path.join(os.path.dirname(__file__), "branch_length_regression.png")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved to {output_path}")

    # Check residuals for nonlinearity
    print("\n" + "=" * 60)
    print("RESIDUAL ANALYSIS")
    print("=" * 60)

    resid = df["dasm_branch_length"] - pred_mut

    print("\nResiduals by decile:")
    df["decile"] = pd.qcut(df["mut_freq"], 10, labels=False)
    for decile in range(10):
        mask = df["decile"] == decile
        mean_resid = resid[mask].mean()
        print(f"  Decile {decile}: mean residual = {mean_resid:+.6f}")

    return model_mut, train_df, val_df


if __name__ == "__main__":
    main()
