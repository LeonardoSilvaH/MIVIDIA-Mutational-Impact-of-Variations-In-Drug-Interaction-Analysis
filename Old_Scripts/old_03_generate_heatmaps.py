#!/usr/bin/env python3

"""
Generate one heatmap figure per protein.

Each figure contains one heatmap for each predicted pocket.

Input
-----
output/delta_scores_long.csv

Output
------
output/heatmaps/
    KRAS_heatmap.png
    MYC_heatmap.png
    TP53_heatmap.png
    ...
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ----------------------------------------------------------
# Configuration
# ----------------------------------------------------------

INPUT = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Analysis/02_compute_delta_scores/delta_scores_long.csv")

OUTPUT = Path("output/heatmaps")
OUTPUT.mkdir(exist_ok=True)

sns.set_theme(style="white")

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------

df = pd.read_csv(INPUT)

# Remove WT rows (their Δ is always zero)
df = df[df["Mutation"] != "WILD-TYE"].copy()

# Fixed color scale so all proteins are comparable
VMIN = df["DeltaScore"].quantile(0.02)
VMAX = df["DeltaScore"].quantile(0.98)

# ----------------------------------------------------------
# Generate one figure per protein
# ----------------------------------------------------------

for protein in sorted(df["Protein"].unique()):

    protein_df = df[df["Protein"] == protein].copy()

    ####################################################
    # Compute one drug ordering for the whole protein
    ####################################################

    drug_order = (
        protein_df
        .groupby("Drug")["DeltaScore"]
        .apply(lambda x: x.abs().mean())
        .sort_values(ascending=False)
        .index
    )

    pockets = sorted(protein_df["PocketRank"].unique())

    fig, axes = plt.subplots(
        1,
        len(pockets),
        figsize=(8 * len(pockets), 10),
        sharey=True
    )

    if len(pockets) == 1:
        axes = [axes]

    for ax, pocket in zip(axes, pockets):

        pocket_df = protein_df[
            protein_df["PocketRank"] == pocket
        ]

        ####################################################
        # Mutation x Drug matrix
        ####################################################

        matrix = pocket_df.pivot_table(
            index="Mutation",
            columns="Drug",
            values="DeltaScore",
            aggfunc="first"
        )

        ####################################################
        # Order mutations
        #
        # WT would normally be first, but we removed it.
        # Remaining mutations are sorted alphabetically.
        ####################################################

        matrix = matrix.sort_index()

        ####################################################
        # Order drugs by mean absolute effect
        ####################################################

        drug_order = (
            matrix.abs()
                  .mean(axis=0)
                  .sort_values(ascending=False)
                  .index
        )

        matrix = matrix[drug_order]

        ####################################################
        # Plot
        ####################################################

        sns.heatmap(
            matrix,
            ax=ax,
            cmap="coolwarm",
            center=0,
            vmin=VMIN,
            vmax=VMAX,
            linewidths=0.0,
            cbar=(pocket == pockets[-1]),
            cbar_kws={
                "label": "Δ Docking Score\n(Mutation − WT)"
            }
        )

        ax.set_title(f"Pocket {pocket}", fontsize=14)

        ax.set_xlabel("Drug")

        if pocket == pockets[0]:
            ax.set_ylabel("Mutation")
        else:
            ax.set_ylabel("")

        ax.tick_params(axis="x", labelrotation=90, labelsize=7)
        ax.tick_params(axis="y", labelsize=8)

    ########################################################
    # Figure title
    ########################################################

    fig.suptitle(
        protein,
        fontsize=20,
        fontweight="bold"
    )

    plt.tight_layout()

    plt.savefig(
        OUTPUT / f"{protein}_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

print("=" * 60)
print(f"Saved {len(df['Protein'].unique())} heatmaps to:")
print(OUTPUT)
print("=" * 60)