import pandas as pd
from scipy.stats import pearsonr, spearmanr


# ==========================================================
# Protein statistics
# ==========================================================

def protein_statistics(df):

    return (
        df.groupby("Protein")
        .agg(
            NDrugs=("Drug", "nunique"),
            NMutations=("Mutation", "nunique"),
            MeanDelta=("DeltaDockingScore", "mean"),
            MeanAbsDelta=("DeltaDockingScore", lambda x: x.abs().mean()),
            Median=("DeltaDockingScore", "median"),
            Std=("DeltaDockingScore", "std"),
            MaxAbsDelta=("DeltaDockingScore", lambda x: x.abs().max()),
        )
        .reset_index()
        .sort_values("MeanAbsDelta", ascending=False)
    )


# ==========================================================
# Drug statistics
# ==========================================================

def drug_statistics(df):

    return (
        df.groupby("Drug")
        .agg(
            NProteins=("Protein", "nunique"),
            NMutations=("Mutation", "nunique"),
            MeanDelta=("DeltaDockingScore", "mean"),
            MeanAbsDelta=("DeltaDockingScore", lambda x: x.abs().mean()),
            Std=("DeltaDockingScore", "std"),
            MaxAbsDelta=("DeltaDockingScore", lambda x: x.abs().max()),
        )
        .reset_index()
        .sort_values("MeanAbsDelta", ascending=False)
    )


# ==========================================================
# Mutation statistics
# ==========================================================

def mutation_statistics(df):

    return (
        df.groupby(["Protein", "Mutation"])
        .agg(
            NDrugs=("Drug", "nunique"),
            MeanDelta=("DeltaDockingScore", "mean"),
            MeanAbsDelta=("DeltaDockingScore", lambda x: x.abs().mean()),
            Std=("DeltaDockingScore", "std"),
            MaxAbsDelta=("DeltaDockingScore", lambda x: x.abs().max()),
        )
        .reset_index()
        .sort_values("MeanAbsDelta", ascending=False)
    )


# ==========================================================
# Affinity summary
# ==========================================================

def classify_affinity(delta):

    if delta <= -1:
        return "Stronger"

    elif delta >= 1:
        return "Weaker"

    return "Neutral"


def affinity_summary(df):

    tmp = df.copy()

    tmp["AffinityChange"] = tmp["DeltaDockingScore"].apply(classify_affinity)

    summary = (
        tmp["AffinityChange"]
        .value_counts()
        .rename_axis("Class")
        .reset_index(name="Count")
    )

    summary["Percentage"] = (
        summary["Count"] / summary["Count"].sum() * 100
    )

    return summary


# ==========================================================
# Pocket transitions
# ==========================================================

def pocket_transition_matrix(df):

    return pd.crosstab(
        df["WT_PocketRank"],
        df["Mut_PocketRank"],
    )


# ==========================================================
# Distance statistics
# ==========================================================

def distance_statistics(df):

    stats = (
        df["Distance"]
        .describe()
        .to_frame()
        .T
    )

    return stats


# ==========================================================
# Correlations
# ==========================================================

def correlation_statistics(df):

    abs_delta = df["DeltaDockingScore"].abs()

    pearson_r, pearson_p = pearsonr(
        df["Distance"],
        abs_delta,
    )

    spearman_rho, spearman_p = spearmanr(
        df["Distance"],
        abs_delta,
    )

    return pd.DataFrame(
        {
            "Metric": [
                "Pearson",
                "Spearman",
            ],
            "Correlation": [
                pearson_r,
                spearman_rho,
            ],
            "PValue": [
                pearson_p,
                spearman_p,
            ],
        }
    )