#!/usr/bin/env python3

"""
Compute ΔDocking Scores (Mutation - WT)

Input
-----
output/pocket_scores_wide.csv

Outputs
-------
output/delta_scores_wide.csv
output/delta_scores_long.csv
"""

from pathlib import Path
import pandas as pd
import re

INPUT = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Analysis/01_extract_best_pose_per_pocket/pocket_scores_wide.csv")

OUTPUT_WIDE = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Analysis/02_compute_delta_scores/delta_scores_wide.csv")
OUTPUT_LONG = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Analysis/02_compute_delta_scores/delta_scores_long.csv")

df = pd.read_csv(INPUT)

##############################################################
# Identify WT rows
##############################################################

wt = df[
    df["Mutation"].str.upper().str.contains("WILD")
].copy()

if wt.empty:
    raise RuntimeError("No WT rows found.")

##############################################################
# Rename WT columns
##############################################################

rename = {}

for col in df.columns:

    if col.startswith("P") and col.endswith("Score"):

        rename[col] = f"WT_{col}"

wt = wt.rename(columns=rename)

keep_cols = ["Protein", "Drug"] + list(rename.values())

wt = wt[keep_cols]

##############################################################
# Merge WT with every mutation
##############################################################

merged = df.merge(
    wt,
    on=["Protein", "Drug"],
    how="left",
    validate="many_to_one"
)

##############################################################
# Compute deltas
##############################################################

pockets = sorted(
    [
        c for c in df.columns
        if re.fullmatch(r"P\d+Score", c)
    ]
)

for pocket in pockets:

    wt_col = f"WT_{pocket}"

    delta = pocket.replace("Score", "Delta")

    merged[delta] = merged[pocket] - merged[wt_col]

##############################################################
# Save wide table
##############################################################

merged.to_csv(
    OUTPUT_WIDE,
    index=False
)

##############################################################
# Build long table
##############################################################

rows = []

for _, row in merged.iterrows():

    for pocket in pockets:

        rank = pocket.replace("Score", "").replace("P", "")

        rows.append({

            "Protein": row["Protein"],
            "Drug": row["Drug"],
            "Mutation": row["Mutation"],

            "PocketRank": int(rank),

            "WTScore": row[f"WT_{pocket}"],

            "MutantScore": row[pocket],

            "DeltaScore": row[pocket.replace("Score", "Delta")],

            "PocketScore":
                row.get(f"P{rank}PocketScore", None),

            "PocketProgram":
                row.get(f"P{rank}PocketProgram", None),

            "OriginalPocketID":
                row.get(f"P{rank}OriginalPocketID", None)

        })

long = pd.DataFrame(rows)

long.to_csv(
    OUTPUT_LONG,
    index=False
)

##############################################################
# Summary
##############################################################

print("=" * 60)
print("Delta scores computed")
print("=" * 60)
print(f"Interactions : {merged.shape[0]}")
print(f"Long rows    : {long.shape[0]}")
print("=" * 60)

print("\nSaved:")
print("  output/delta_scores_wide.csv")
print("  output/delta_scores_long.csv")