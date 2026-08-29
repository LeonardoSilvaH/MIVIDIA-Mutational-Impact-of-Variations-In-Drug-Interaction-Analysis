#!/usr/bin/env python3

import json
import os
import pandas as pd


# ==========================================================
# INPUT / OUTPUT
# ==========================================================

POSE_JSON = (
    "Data/kras_article/Adagrasib_Daraxonrasib_results/merged_vina_results.json"
)

OUTPUT_CSV = (
    "Analysis/1.5_Docking_csv/all_docking_poses.csv"
)


# ==========================================================
# LOAD JSON
# ==========================================================

print()
print("Loading merged Vina results...")
print(f"Input: {POSE_JSON}")

with open(POSE_JSON, "r") as f:
    data = json.load(f)


# ==========================================================
# VALIDATE JSON
# ==========================================================

if "docking_runs" not in data:

    raise ValueError(
        "Invalid Vina JSON: 'docking_runs' was not found."
    )


docking_runs = data["docking_runs"]

print(
    f"Found {len(docking_runs)} docking runs."
)


# ==========================================================
# PARSE
# ==========================================================

rows = []


for run_id, run_data in docking_runs.items():

    # ------------------------------------------------------
    # Metadata
    # ------------------------------------------------------

    metadata = run_data.get(
        "metadata",
        {}
    )

    protein = metadata.get(
        "protein"
    )

    drug = metadata.get(
        "drug"
    )

    mutation_name = metadata.get(
        "mutation"
    )

    chain = metadata.get(
        "chain"
    )


    # ------------------------------------------------------
    # Mutation without protein prefix
    #
    # KRAS-G12A -> G12A
    #
    # But preserve the original mutation in
    # MutationFolder.
    # ------------------------------------------------------

    mutation = mutation_name

    if mutation and "-" in mutation:

        mutation = mutation.split(
            "-"
        )[-1]


    # ------------------------------------------------------
    # Source information
    # ------------------------------------------------------

    source = run_data.get(
        "source",
        {}
    )

    source_file = source.get(
        "file"
    )

    source_absolute_path = source.get(
        "absolute_path"
    )


    # ------------------------------------------------------
    # Docking information
    # ------------------------------------------------------

    docking = run_data.get(
        "docking",
        {}
    )

    number_of_modes = docking.get(
        "number_of_modes"
    )

    best_mode_data = docking.get(
        "best_mode",
        {}
    )

    best_mode = best_mode_data.get(
        "mode"
    )

    best_energy = best_mode_data.get(
        "energy"
    )

    best_pKd = best_mode_data.get(
        "pKd"
    )


    # ------------------------------------------------------
    # Individual poses
    # ------------------------------------------------------

    modes = docking.get(
        "modes",
        {}
    )


    for pose_id, pose in modes.items():

        row = {

            # ==============================================
            # Provenance
            # ==============================================

            "RunID": run_id,

            "SourceFile": source_file,

            "SourceAbsolutePath": source_absolute_path,

            "DockingEngine": data.get(
                "docking_engine"
            ),


            # ==============================================
            # Biological identity
            # ==============================================

            "Protein": protein,

            "Drug": drug,

            "Mutation": mutation,

            "MutationFolder": mutation_name,

            "Chain": chain,


            # ==============================================
            # Docking information
            # ==============================================

            "NumberOfModes": number_of_modes,

            "PoseID": int(
                pose.get(
                    "mode",
                    pose_id
                )
            ),

            "DockingScore": pose.get(
                "energy"
            ),

            "pKd": pose.get(
                "pKd"
            ),

            "RMSD_LB": pose.get(
                "rmsd_lb"
            ),

            "RMSD_UB": pose.get(
                "rmsd_ub"
            ),


            # ==============================================
            # Best pose information
            #
            # Useful for identifying whether this row
            # corresponds to the best Vina mode.
            # ==============================================

            "BestMode": best_mode,

            "BestEnergy": best_energy,

            "BestpKd": best_pKd,

            "IsBestMode": (
                pose.get(
                    "mode",
                    pose_id
                ) == best_mode
            ),

        }

        rows.append(row)


# ==========================================================
# CREATE DATAFRAME
# ==========================================================

df = pd.DataFrame(
    rows
)


if df.empty:

    raise ValueError(
        "No docking poses were extracted."
    )


# ==========================================================
# NUMERIC CONVERSION
# ==========================================================

numeric_columns = [

    "PoseID",

    "NumberOfModes",

    "DockingScore",

    "pKd",

    "RMSD_LB",

    "RMSD_UB",

    "BestMode",

    "BestEnergy",

    "BestpKd",

]


for column in numeric_columns:

    if column in df.columns:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce"
        )


# ==========================================================
# SORT
# ==========================================================

df = df.sort_values(
    [
        "Protein",
        "Drug",
        "Mutation",
        "Chain",
        "PoseID",
    ]
)


# ==========================================================
# OUTPUT DIRECTORY
# ==========================================================

output_directory = os.path.dirname(
    OUTPUT_CSV
)

os.makedirs(
    output_directory,
    exist_ok=True
)


# ==========================================================
# SAVE
# ==========================================================

df.to_csv(
    OUTPUT_CSV,
    index=False
)


# ==========================================================
# SUMMARY
# ==========================================================

print()
print("==============================================")
print("Docking extraction complete")
print("==============================================")
print()

print(
    f"Docking runs : {df['RunID'].nunique():,}"
)

print(
    f"Proteins     : {df['Protein'].nunique():,}"
)

print(
    f"Drugs        : {df['Drug'].nunique():,}"
)

print(
    f"Mutations    : {df['Mutation'].nunique():,}"
)

print(
    f"Total poses  : {len(df):,}"
)

print()

print(
    f"Best poses   : "
    f"{df['IsBestMode'].sum():,}"
)

print()

print(
    f"Saved to: {OUTPUT_CSV}"
)

print()