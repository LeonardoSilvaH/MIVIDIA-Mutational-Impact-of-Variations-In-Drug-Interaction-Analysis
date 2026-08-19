#!/usr/bin/env python3

"""
Extract the best docking pose for each predicted pocket from CoBDock.

Inputs
------
input/
    pose_data.json
    top_3_pockets.json

Outputs
-------
output/
    pocket_scores_long.csv
    pocket_scores_wide.csv
"""

from pathlib import Path
import json
import pandas as pd

# ==========================================================
# Paths
# ==========================================================

INPUT_DIR = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Data/results_semifinal/pocket_locations")
OUTPUT_DIR = Path("/home/leo/Documents/GitHub/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis/Analysis/01_extract_best_pose_per_pocket")

POSE_FILE = INPUT_DIR / "pose_data.json"
POCKET_FILE = INPUT_DIR / "top_3_pockets.json"

OUTPUT_DIR.mkdir(exist_ok=True)

# ==========================================================
# Load JSONs
# ==========================================================

with open(POSE_FILE) as f:
    pose_data = json.load(f)

with open(POCKET_FILE) as f:
    pocket_data = json.load(f)

rows = []

missing_pockets = 0

# ==========================================================
# Parse all interactions
# ==========================================================

for interaction, mutations in pose_data.items():

    protein, drug = interaction.split("_", 1)

    if interaction not in pocket_data:
        print(f"WARNING: {interaction} not found in top_3_pockets.json")
        continue

    pocket_predictions = pocket_data[interaction]

    for mutation, chains in mutations.items():

        if mutation not in pocket_predictions:
            print(f"WARNING: {interaction} -> {mutation} missing from top_3_pockets.json")
            continue

        for chain_name, pockets in chains.items():

            if chain_name not in pocket_predictions[mutation]:
                print(f"WARNING: {interaction} -> {mutation} -> {chain_name} missing")
                continue

            predicted = pocket_predictions[mutation][chain_name]

            # Sort by PocketRank just to be safe
            predicted = sorted(predicted,
                               key=lambda x: x["pocket_rank"])

            for pocket in predicted:

                pocket_rank = int(pocket["pocket_rank"])

                # IMPORTANT:
                # pose_data.json uses PocketRank as key
                docking_key = str(pocket_rank)

                if docking_key not in pockets:

                    missing_pockets += 1

                    print(
                        f"WARNING: Missing docking results\n"
                        f" Protein : {protein}\n"
                        f" Drug    : {drug}\n"
                        f" Mutation: {mutation}\n"
                        f" Chain   : {chain_name}\n"
                        f" PocketRank: {pocket_rank}\n"
                    )

                    continue

                docking_poses = pockets[docking_key]

                best_pose = None
                best_score = float("inf")

                for pose_id, pose in docking_poses.items():

                    score = pose["total_score"]

                    if score < best_score:
                        best_score = score
                        best_pose = int(pose_id)

                rows.append({

                    "Protein": protein,
                    "Drug": drug,
                    "Mutation": mutation,

                    "InteractionID":
                        f"{protein}|{drug}|{mutation}",

                    "PocketRank": pocket_rank,

                    # key used inside pose_data.json
                    "DockingPocket": docking_key,

                    # original fpocket/p2rank ID
                    "OriginalPocketID": pocket["selected_pocket_id"],

                    "PocketProgram":
                        pocket["selected_pocket_program"],

                    "PocketScore":
                        pocket["pocket_score"],

                    "BestPose":
                        best_pose,

                    "BestScore":
                        best_score,

                    "NumPoses":
                        len(docking_poses)

                })

# ==========================================================
# Long dataframe
# ==========================================================

long_df = pd.DataFrame(rows)

long_df = long_df.sort_values(
    [
        "Protein",
        "Drug",
        "Mutation",
        "PocketRank"
    ]
)

long_df.to_csv(
    OUTPUT_DIR / "pocket_scores_long.csv",
    index=False
)

# ==========================================================
# Wide dataframe
# ==========================================================

wide_rows = []

for (protein, drug, mutation), group in long_df.groupby(
    ["Protein", "Drug", "Mutation"]
):

    row = {

        "Protein": protein,
        "Drug": drug,
        "Mutation": mutation,
        "InteractionID": f"{protein}|{drug}|{mutation}"

    }

    for _, pocket in group.iterrows():

        rank = int(pocket["PocketRank"])

        row[f"P{rank}Score"] = pocket["BestScore"]
        row[f"P{rank}PocketScore"] = pocket["PocketScore"]
        row[f"P{rank}PocketProgram"] = pocket["PocketProgram"]
        row[f"P{rank}OriginalPocketID"] = pocket["OriginalPocketID"]
        row[f"P{rank}DockingPocket"] = pocket["DockingPocket"]
        row[f"P{rank}BestPose"] = pocket["BestPose"]

    wide_rows.append(row)

wide_df = pd.DataFrame(wide_rows)

wide_df = wide_df.sort_values(
    [
        "Protein",
        "Drug",
        "Mutation"
    ]
)

wide_df.to_csv(
    OUTPUT_DIR / "pocket_scores_wide.csv",
    index=False
)

# ==========================================================
# Summary
# ==========================================================

print("=" * 70)
print("Extraction complete")
print("=" * 70)
print(f"Interactions : {wide_df.shape[0]}")
print(f"Pocket rows  : {long_df.shape[0]}")
print(f"Missing pockets : {missing_pockets}")
print("=" * 70)

print("\nFiles written:")
print("  output/pocket_scores_long.csv")
print("  output/pocket_scores_wide.csv")