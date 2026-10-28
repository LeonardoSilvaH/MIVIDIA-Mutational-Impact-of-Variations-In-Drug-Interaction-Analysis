#!/usr/bin/env python3

import json
import pandas as pd
from tqdm import tqdm
from pathlib import Path


# ==========================================================
# INPUT / OUTPUT
# ==========================================================

MATCHED_POCKETS_FILE = (
    "Analysis/02_matched_pockets/02_matched_pockets.csv"
)

# ----------------------------------------------------------
# OLD CoBDock input
# ----------------------------------------------------------

POCKET_LOCATIONS_DIR = Path(
    "/hd/leo/MIVIDIA/Data/resultados_mestrado/pocket_locations"
)

POSE_DATA_FILE = (
    POCKET_LOCATIONS_DIR /
    "pose_data.json"
)

# ----------------------------------------------------------
# NEW fpocket + Vina input
# ----------------------------------------------------------

VINA_DATA_FILE = (
    "Data/kras_article/Adagrasib_Daraxonrasib_results/merged_vina_results.json"
)

# ----------------------------------------------------------
# Outputs
# ----------------------------------------------------------

OUTPUT_DIR = Path(
    "Analysis/03_compute_delta_scores"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_FULL = (
    OUTPUT_DIR /
    "matched_pockets_scored.csv"
)

OUTPUT_SUMMARY = (
    OUTPUT_DIR /
    "docking_summary.csv"
)


# ==========================================================
# HELPERS
# ==========================================================

def load_json(filename):

    with open(
        filename,
        "r"
    ) as f:

        return json.load(f)


# ==========================================================
# DETECT INPUT MODE
# ==========================================================

def detect_input_mode():

    pose_exists = Path(
        POSE_DATA_FILE
    ).exists()

    vina_exists = Path(
        VINA_DATA_FILE
    ).exists()

    if pose_exists:

        print(
            "Input mode: STANDARD CoBDock / PLANTS"
        )

        return "plants"

    elif vina_exists:

        print(
            "Input mode: fpocket + global Vina"
        )

        return "vina"

    else:

        raise FileNotFoundError(
            "\nNo valid docking input found.\n\n"
            "Expected either:\n"
            f"  {POSE_DATA_FILE}\n"
            "or:\n"
            f"  {VINA_DATA_FILE}\n"
        )


# ==========================================================
# STANDARD PLANTS FUNCTIONS
# ==========================================================

def get_chain_name(mutation_dict):
    """
    Returns the only chain available.
    """

    return next(
        iter(
            mutation_dict.keys()
        )
    )


def get_best_pose(pocket_dict):
    """
    Find the pose with the lowest PLANTS total_score.

    Returns
    -------
    best_pose : str
    best_score : float
    """

    best_pose = None

    best_score = float(
        "inf"
    )

    for pose_id, pose_data in pocket_dict.items():

        # Skip poses without a numeric score (failed or empty docking)
        try:
            score = float(
                pose_data["total_score"]
            )
        except (
            KeyError,
            TypeError,
            ValueError
        ):
            continue

        if score < best_score:

            best_score = score

            best_pose = pose_id

    if best_pose is None:

        raise ValueError(
            "No pose with a valid PLANTS total_score in this pocket "
            f"({len(pocket_dict)} poses found)"
        )

    return (
        best_pose,
        best_score
    )


def get_pocket(
    pose_data,
    protein,
    drug,
    mutation,
    pocket_rank,
):
    """
    Returns the dictionary corresponding
    to one pocket.
    """

    complex_name = (
        f"{protein}_{drug}"
    )

    mutation_dict = (
        pose_data[
            complex_name
        ][
            mutation
        ]
    )

    chain = get_chain_name(
        mutation_dict
    )

    return (
        mutation_dict[
            chain
        ][
            str(
                int(
                    pocket_rank
                )
            )
        ]
    )


# ==========================================================
# VINA FUNCTIONS
# ==========================================================

def get_vina_run(
    vina_data,
    protein,
    drug,
    mutation,
    chain="A"
):
    """
    Retrieve one docking run from
    merged_vina_results.json.

    Expected run ID:

        KRAS_Adagrasib_KRAS-G12A_A
    """

    run_id = (
        f"{protein}_"
        f"{drug}_"
        f"{mutation}_"
        f"{chain}"
    )

    runs = vina_data.get(
        "docking_runs",
        {}
    )

    return runs.get(
        run_id
    )


def get_best_vina_score(
    vina_data,
    protein,
    drug,
    mutation,
    chain="A"
):
    """
    Retrieve the best global Vina energy
    for one docking run.

    This is NOT a pocket-specific score.
    """

    run = get_vina_run(
        vina_data,
        protein,
        drug,
        mutation,
        chain
    )

    if run is None:

        raise KeyError(
            "Vina docking run not found: "
            f"{protein}_{drug}_{mutation}_{chain}"
        )

    docking = run.get(
        "docking",
        {}
    )

    best_mode = docking.get(
        "best_mode"
    )

    if best_mode is None:

        raise KeyError(
            "best_mode not found for: "
            f"{protein}_{drug}_{mutation}_{chain}"
        )

    energy = best_mode.get(
        "energy"
    )

    if energy is None:

        raise KeyError(
            "Energy not found for: "
            f"{protein}_{drug}_{mutation}_{chain}"
        )

    return float(
        energy
    )


# ==========================================================
# PROCESS ONE MATCH — PLANTS
# ==========================================================

def process_match_plants(
    row,
    pose_data
):

    protein = row["Protein"]

    drug = row["Drug"]

    mutation = row["Mutation"]

    # WT key is "KRAS-WILD-TYE" in older runs and
    # "KRAS.WILD-TYE" in newer ones, so look it up.
    wt_mutation = next(
        key
        for key in pose_data[f"{protein}_{drug}"]
        if "WILD" in key.upper()
    )

    wt_rank = row[
        "WT_PocketRank"
    ]

    mut_rank = row[
        "Mut_PocketRank"
    ]

    # ------------------------------------------------------
    # WT pocket
    # ------------------------------------------------------

    wt_pocket = get_pocket(
        pose_data,
        protein,
        drug,
        wt_mutation,
        wt_rank,
    )

    # ------------------------------------------------------
    # Mutant pocket
    # ------------------------------------------------------

    mut_pocket = get_pocket(
        pose_data,
        protein,
        drug,
        mutation,
        mut_rank,
    )

    # ------------------------------------------------------
    # Best poses
    # ------------------------------------------------------

    wt_pose, wt_score = (
        get_best_pose(
            wt_pocket
        )
    )

    mut_pose, mut_score = (
        get_best_pose(
            mut_pocket
        )
    )

    # ------------------------------------------------------
    # Store scores
    # ------------------------------------------------------

    row["WT_BestPose"] = int(
        wt_pose
    )

    row["WT_BestDockingScore"] = (
        wt_score
    )

    row["Mut_BestPose"] = int(
        mut_pose
    )

    row["Mut_BestDockingScore"] = (
        mut_score
    )

    # ------------------------------------------------------
    # Delta
    # ------------------------------------------------------

    delta = (
        mut_score -
        wt_score
    )

    row["DeltaDockingScore"] = delta

    row["AbsDeltaDockingScore"] = abs(
        delta
    )

    # ------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------

    if delta < 0:

        row["AffinityChange"] = (
            "Stronger"
        )

    elif delta > 0:

        row["AffinityChange"] = (
            "Weaker"
        )

    else:

        row["AffinityChange"] = (
            "Unchanged"
        )

    row["DockingScoreType"] = (
        "PLANTS_PocketSpecific"
    )

    return row


# ==========================================================
# PROCESS ONE MATCH — VINA
# ==========================================================

def process_match_vina(
    row,
    vina_data
):

    protein = row["Protein"]

    drug = row["Drug"]

    mutation = row["Mutation"]

    # ------------------------------------------------------
    # WT mutation
    # ------------------------------------------------------

    wt_mutation = (
        f"{protein}-WILD-TYE"
    )

    # ------------------------------------------------------
    # Chain
    #
    # Your current dataset uses chain A.
    # ------------------------------------------------------

    chain = "A"

    # ------------------------------------------------------
    # Retrieve global Vina scores
    # ------------------------------------------------------

    wt_score = get_best_vina_score(
        vina_data,
        protein,
        drug,
        wt_mutation,
        chain
    )

    mut_score = get_best_vina_score(
        vina_data,
        protein,
        drug,
        mutation,
        chain
    )

    # ------------------------------------------------------
    # Store global scores
    # ------------------------------------------------------

    row["WT_GlobalVinaScore"] = (
        wt_score
    )

    row["Mut_GlobalVinaScore"] = (
        mut_score
    )

    # ------------------------------------------------------
    # Calculate delta
    # ------------------------------------------------------

    delta = (
        mut_score -
        wt_score
    )

    row["DeltaGlobalVinaScore"] = (
        delta
    )

    row["AbsDeltaGlobalVinaScore"] = (
        abs(delta)
    )

    # ------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------

    if delta < 0:

        row["GlobalAffinityChange"] = (
            "Stronger"
        )

    elif delta > 0:

        row["GlobalAffinityChange"] = (
            "Weaker"
        )

    else:

        row["GlobalAffinityChange"] = (
            "Unchanged"
        )

    # ------------------------------------------------------
    # Explicitly document score type
    # ------------------------------------------------------

    row["DockingScoreType"] = (
        "Vina_GlobalBestMode"
    )

    return row


# ==========================================================
# MAIN
# ==========================================================

def main():

    print(
        "Loading data..."
    )

    # ------------------------------------------------------
    # Matched pockets
    # ------------------------------------------------------

    matched_df = pd.read_csv(
        MATCHED_POCKETS_FILE
    )

    print(
        f"Matched pocket rows: "
        f"{len(matched_df)}"
    )

    # ------------------------------------------------------
    # Detect mode
    # ------------------------------------------------------

    input_mode = (
        detect_input_mode()
    )

    # ------------------------------------------------------
    # Load appropriate JSON
    # ------------------------------------------------------

    if input_mode == "plants":

        print(
            "Loading pose_data.json..."
        )

        docking_data = load_json(
            POSE_DATA_FILE
        )

    else:

        print(
            "Loading merged_vina_results.json..."
        )

        docking_data = load_json(
            VINA_DATA_FILE
        )

    # ------------------------------------------------------
    # Process
    # ------------------------------------------------------

    processed_rows = []

    errors = []

    print(
        f"\nProcessing "
        f"{len(matched_df)} matched pockets...\n"
    )

    for index, row in tqdm(
        matched_df.iterrows(),
        total=len(matched_df)
    ):

        try:

            if input_mode == "plants":

                processed_row = (
                    process_match_plants(
                        row.copy(),
                        docking_data
                    )
                )

            else:

                processed_row = (
                    process_match_vina(
                        row.copy(),
                        docking_data
                    )
                )

            processed_rows.append(
                processed_row
            )

        except Exception as e:

            error_message = (
                f"{row['Protein']} | "
                f"{row['Drug']} | "
                f"{row['Mutation']} | "
                f"{e}"
            )

            errors.append(
                error_message
            )

            print(
                f"\nERROR: {error_message}"
            )

    # ------------------------------------------------------
    # DataFrame
    # ------------------------------------------------------

    scored_df = pd.DataFrame(
        processed_rows
    )

    if scored_df.empty:

        raise RuntimeError(
            "No rows were successfully processed."
        )

    # ------------------------------------------------------
    # Sort
    # ------------------------------------------------------

    scored_df = scored_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank"
        ]
    )

    # ------------------------------------------------------
    # Save full table
    # ------------------------------------------------------

    print(
        "\nSaving complete table..."
    )

    scored_df.to_csv(
        OUTPUT_FULL,
        index=False
    )

    # ======================================================
    # SUMMARY
    # ======================================================

    print(
        "Generating summary table..."
    )

    # ------------------------------------------------------
    # Important:
    #
    # In Vina mode, the score is GLOBAL and therefore
    # identical across the three matched pockets of the
    # same mutation/drug pair.
    #
    # We therefore select the first matched pocket rather
    # than sorting by docking score.
    # ------------------------------------------------------

    if input_mode == "vina":

        summary_df = (
            scored_df
            .sort_values(
                [
                    "Protein",
                    "Drug",
                    "Mutation",
                    "WT_PocketRank"
                ]
            )
            .groupby(
                [
                    "Protein",
                    "Drug",
                    "Mutation"
                ],
                as_index=False
            )
            .first()
        )

    else:

        summary_df = (
            scored_df
            .sort_values(
                "Mut_BestDockingScore"
            )
            .groupby(
                [
                    "Protein",
                    "Drug",
                    "Mutation"
                ],
                as_index=False
            )
            .first()
        )

    summary_df.to_csv(
        OUTPUT_SUMMARY,
        index=False
    )

    # ======================================================
    # FINAL REPORT
    # ======================================================

    print()
    print(
        "=" * 60
    )

    print(
        "DONE"
    )

    print(
        "=" * 60
    )

    print(
        f"Input mode        : {input_mode}"
    )

    print(
        f"Input rows        : {len(matched_df)}"
    )

    print(
        f"Processed rows    : {len(scored_df)}"
    )

    print(
        f"Errors            : {len(errors)}"
    )

    print()

    print(
        f"Full output:"
    )

    print(
        f"  {OUTPUT_FULL}"
    )

    print()

    print(
        f"Summary output:"
    )

    print(
        f"  {OUTPUT_SUMMARY}"
    )

    print(
        "=" * 60
    )

    # ------------------------------------------------------
    # Print errors
    # ------------------------------------------------------

    if errors:

        print()
        print(
            "Errors encountered:"
        )

        for error in errors:

            print(
                f"  - {error}"
            )


# ==========================================================
# ENTRY POINT
# ==========================================================

if __name__ == "__main__":

    main()