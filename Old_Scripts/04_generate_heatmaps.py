#!/usr/bin/env python3

"""
Generate one heatmap per protein.

Each figure contains three heatmaps:
    Pocket 1
    Pocket 2
    Pocket 3

Rows:
    Mutations

Columns:
    Drugs

Cell color:
    Δ Docking Score (Mutation - WT)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ==========================================================
# Configuration
# ==========================================================

INPUT = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Analysis/02_compute_delta_scores/delta_scores_long.csv")

OUTPUT = Path("output/heatmaps")
OUTPUT.mkdir(exist_ok=True)

sns.set_theme(style="white")

# ==========================================================
# Load data
# ==========================================================

df = pd.read_csv(INPUT)

############################################################
# Remove WT rows
############################################################

df = df[df["Mutation"] != "WILD-TYE"].copy()

############################################################
# Global color scale
############################################################

VMIN = df["DeltaScore"].quantile(0.02)
VMAX = df["DeltaScore"].quantile(0.98)

print(f"Color scale: {VMIN:.2f} -> {VMAX:.2f}")

# ==========================================================
# Generate one figure per protein
# ==========================================================

for protein in sorted(df["Protein"].unique()):

    print(f"Generating {protein}...")

    protein_df = df[df["Protein"] == protein].copy()

    ########################################################
    # Order drugs
    #
    # Uses the maximum absolute ΔScore across ALL pockets.
    ########################################################

    drug_order = (
        protein_df
        .groupby("Drug")["DeltaScore"]
        .apply(lambda x: x.abs().max())
        .sort_values(ascending=False)
        .index
    )

    ########################################################
    # Order mutations
    #
    # Uses the maximum absolute ΔScore across ALL pockets.
    ########################################################

    mutation_order = (
        protein_df
        .groupby("Mutation")["DeltaScore"]
        .apply(lambda x: x.abs().max())
        .sort_values(ascending=False)
        .index
    )

    ########################################################
    # Pockets
    ########################################################

    pockets = sorted(protein_df["PocketRank"].unique())

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(pockets),
        figsize=(9 * len(pockets), 11),
        sharey=True
    )

    if len(pockets) == 1:
        axes = [axes]

    ########################################################
    # Draw one heatmap per pocket
    ########################################################

    for ax, pocket in zip(axes, pockets):

        pocket_df = protein_df[
            protein_df["PocketRank"] == pocket
        ]

        matrix = pocket_df.pivot_table(
            index="Mutation",
            columns="Drug",
            values="DeltaScore",
            aggfunc="first"
        )

        ####################################################
        # Apply SAME ordering to every pocket
        ####################################################

        matrix = matrix.reindex(
            index=mutation_order,
            columns=drug_order
        )

        ####################################################
        # Plot
        ####################################################

        hm = sns.heatmap(
            matrix,
            ax=ax,
            cmap="coolwarm",
            center=0,
            vmin=VMIN,
            vmax=VMAX,
            linewidths=0,
            linecolor="white",
            cbar=(pocket == pockets[-1]),
            cbar_kws={
                "label": "Δ Docking Score\n(Mutant − WT)"
            }
        )

        ax.set_title(
            f"Pocket {pocket}",
            fontsize=15,
            weight="bold"
        )

        ax.set_xlabel("Drug", fontsize=12)

        if pocket == pockets[0]:
            ax.set_ylabel("Mutation", fontsize=12)
        else:
            ax.set_ylabel("")

        ####################################################
        # Tick formatting
        ####################################################

        ax.tick_params(
            axis="x",
            rotation=90,
            labelsize=6
        )

        ax.tick_params(
            axis="y",
            labelsize=8
        )

    ########################################################
    # Overall title
    ########################################################

    fig.suptitle(
        protein,
        fontsize=22,
        weight="bold"
    )

    plt.tight_layout()

    plt.savefig(
        OUTPUT / f"{protein}_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

print("\nDone!")
print(f"Heatmaps saved to {OUTPUT}")