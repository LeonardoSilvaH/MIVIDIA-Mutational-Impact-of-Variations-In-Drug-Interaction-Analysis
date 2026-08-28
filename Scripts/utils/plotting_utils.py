import os

import matplotlib.pyplot as plt
import pandas as pd


# ==========================================================
# Internal score helper
# ==========================================================

def get_delta_column(df):

    if "DeltaDockingScore" in df.columns:

        return (
            "DeltaDockingScore",
            "Δ Docking Score"
        )

    elif "DeltaGlobalVinaScore" in df.columns:

        return (
            "DeltaGlobalVinaScore",
            "Δ Global Vina Score"
        )

    else:

        raise KeyError(
            "No valid delta-score column found. "
            "Expected 'DeltaDockingScore' or "
            "'DeltaGlobalVinaScore'."
        )


# ==========================================================
# Distribution
# ==========================================================

def plot_delta_distribution(
    df,
    output_dir
):

    delta_column, delta_label = (
        get_delta_column(df)
    )

    plt.figure(
        figsize=(8, 5)
    )

    plt.hist(
        df[delta_column],
        bins=40,
        edgecolor="black",
    )

    plt.xlabel(
        delta_label
    )

    plt.ylabel(
        "Count"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "delta_distribution.png"
        ),
        dpi=300,
    )

    plt.close()


# ==========================================================
# Protein barplot
# ==========================================================

def plot_protein_barplot(
    stats,
    output_dir
):

    plt.figure(
        figsize=(12, 6)
    )

    plt.bar(
        stats["Protein"],
        stats["MeanAbsDelta"],
    )

    plt.xticks(
        rotation=90
    )

    plt.ylabel(
        "Mean |Δ Score|"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "protein_barplot.png"
        ),
        dpi=300,
    )

    plt.close()


# ==========================================================
# Drug barplot
# ==========================================================

def plot_drug_barplot(
    stats,
    output_dir
):

    plt.figure(
        figsize=(14, 6)
    )

    plt.bar(
        stats["Drug"],
        stats["MeanAbsDelta"],
    )

    plt.xticks(
        rotation=90,
        fontsize=7
    )

    plt.ylabel(
        "Mean |Δ Score|"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "drug_barplot.png"
        ),
        dpi=300,
    )

    plt.close()


# ==========================================================
# Mutation barplot
# ==========================================================

def plot_mutation_barplot(
    stats,
    output_dir
):

    labels = (
        stats["Protein"]
        + "\n"
        + stats["Mutation"]
    )

    plt.figure(
        figsize=(16, 6)
    )

    plt.bar(
        labels,
        stats["MeanAbsDelta"],
    )

    plt.xticks(
        rotation=90,
        fontsize=6
    )

    plt.ylabel(
        "Mean |Δ Score|"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "mutation_barplot.png"
        ),
        dpi=300,
    )

    plt.close()


# ==========================================================
# Pocket transition heatmap
# ==========================================================

def plot_transition_heatmap(
    matrix,
    output_dir
):

    plt.figure(
        figsize=(5, 5)
    )

    plt.imshow(
        matrix.values
    )

    plt.xticks(
        range(
            len(matrix.columns)
        ),
        matrix.columns,
    )

    plt.yticks(
        range(
            len(matrix.index)
        ),
        matrix.index,
    )

    plt.xlabel(
        "Mutant Pocket"
    )

    plt.ylabel(
        "WT Pocket"
    )

    for i in range(
        matrix.shape[0]
    ):

        for j in range(
            matrix.shape[1]
        ):

            plt.text(
                j,
                i,
                matrix.iloc[i, j],
                ha="center",
                va="center",
                fontsize=12,
            )

    plt.colorbar()

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "pocket_transition_heatmap.png"
        ),
        dpi=300,
    )

    plt.close()


# ==========================================================
# Pocket distance distribution
# ==========================================================

def plot_distance_distribution(
    df,
    output_dir
):

    plt.figure(
        figsize=(8, 5)
    )

    plt.hist(
        df["Distance"],
        bins=40,
        edgecolor="black",
    )

    plt.xlabel(
        "Pocket Center Distance (Å)"
    )

    plt.ylabel(
        "Count"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "distance_distribution.png"
        ),
        dpi=300,
    )

    plt.close()


# ==========================================================
# Distance vs delta
# ==========================================================

def plot_distance_vs_delta(
    df,
    output_dir
):

    delta_column, delta_label = (
        get_delta_column(df)
    )

    plt.figure(
        figsize=(6, 6)
    )

    plt.scatter(
        df["Distance"],
        df[delta_column].abs(),
        alpha=0.5,
        s=10,
    )

    plt.xlabel(
        "Pocket Distance (Å)"
    )

    plt.ylabel(
        f"|{delta_label}|"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "distance_vs_delta.png"
        ),
        dpi=300,
    )

    plt.close()