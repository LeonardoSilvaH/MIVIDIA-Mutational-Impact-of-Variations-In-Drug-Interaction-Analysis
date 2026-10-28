#!/usr/bin/env python3

"""
02_match_pockets.py

Matches WT and mutant binding pockets after structural alignment.

Supports two input modes:

1. STANDARD CoBDock MODE
   ----------------------
   top_3_pockets.json
   pose_data.json

2. ALTERNATIVE fpocket + Vina MODE
   ---------------------------------
   fpocket_data.json
   merged_vina_results.json

The script automatically detects which mode is being used.

Outputs
-------
Analysis/02_matched_pockets/02_matched_pockets.csv
Analysis/02_matched_pockets/02_distance_matrix.csv
"""


from pathlib import Path
from dataclasses import dataclass
import json
import pickle

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


###############################################################################
# USER SETTINGS
###############################################################################

ALIGNMENT_FILE = Path(
    "Analysis/01_align_structures/alignment_matrices.pkl"
)

# ---------------------------------------------------------------------------
# Standard Drug Screening CoBDock files
# ---------------------------------------------------------------------------

POSE_JSON = Path(
    "Data/results/pocket_locations/pose_data.json"
)

POCKET_JSON = Path(
    "Data/results/pocket_locations/top_3_pockets.json"
)

# ---------------------------------------------------------------------------
# Alternative fpocket + Vina files
# ---------------------------------------------------------------------------

FPOCKET_JSON = Path(
    "Data/kras_article/Adagrasib_Daraxonrasib_results/fpocket_data.json"
)

VINA_JSON = Path(
    "Data/kras_article/Adagrasib_Daraxonrasib_results/merged_vina_results.json"
)

# ---------------------------------------------------------------------------
# Number of fpocket pockets to use
# ---------------------------------------------------------------------------

TOP_N_FPOCKETS = 3

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path(
    "Analysis/02_matched_pockets"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

MATCHED_OUTPUT = (
    OUTPUT_DIR /
    "02_matched_pockets.csv"
)

DISTANCE_OUTPUT = (
    OUTPUT_DIR /
    "02_distance_matrix.csv"
)

SKIPPED_OUTPUT = (
    OUTPUT_DIR /
    "02_skipped_entries.csv"
)


###############################################################################
# DATA STRUCTURE
###############################################################################

@dataclass
class Pocket:

    rank: int

    pocket_id: str

    pocket_score: float

    center: np.ndarray

    docking_score: float = np.nan

    docking_score_type: str = "NA"


###############################################################################
# LOAD ALIGNMENTS
###############################################################################

print(
    "Loading alignment matrices..."
)

with open(
    ALIGNMENT_FILE,
    "rb"
) as f:

    alignments = pickle.load(f)


###############################################################################
# DETECT INPUT MODE
###############################################################################

standard_mode = (
    POCKET_JSON.exists()
    and POSE_JSON.exists()
)

alternative_mode = (
    FPOCKET_JSON.exists()
    and VINA_JSON.exists()
)


if standard_mode:

    MODE = "standard"

    print(
        "Input mode: STANDARD CoBDock"
    )

elif alternative_mode:

    MODE = "alternative"

    print(
        "Input mode: ALTERNATIVE fpocket + Vina"
    )

else:

    raise FileNotFoundError(
        "\nCould not identify a valid input configuration.\n\n"
        "For standard mode, expected:\n"
        f"  {POCKET_JSON}\n"
        f"  {POSE_JSON}\n\n"
        "For alternative mode, expected:\n"
        f"  {FPOCKET_JSON}\n"
        f"  {VINA_JSON}\n"
    )


###############################################################################
# LOAD JSON FILES
###############################################################################

if MODE == "standard":

    print(
        "Loading top_3_pockets.json..."
    )

    with open(
        POCKET_JSON
    ) as f:

        pocket_json = json.load(f)


    print(
        "Loading pose_data.json..."
    )

    with open(
        POSE_JSON
    ) as f:

        pose_json = json.load(f)


else:

    print(
        "Loading fpocket_data.json..."
    )

    with open(
        FPOCKET_JSON
    ) as f:

        pocket_json = json.load(f)


    print(
        "Loading merged_vina_results.json..."
    )

    with open(
        VINA_JSON
    ) as f:

        vina_json = json.load(f)


###############################################################################
# HELPER FUNCTIONS
###############################################################################

def get_wildtype_key(entry):
    """
    Find the WT entry.

    Examples:

        WILD-TYE
        KRAS-WILD-TYE
        RAF1-WILD-TYE
    """

    for key in entry.keys():

        if "WILD" in key.upper():

            return key

    raise ValueError(
        "Wild-type not found."
    )


def split_ligand_name(name):
    """
    Split:

        MYC_Prednisone

    into:

        MYC
        Prednisone
    """

    protein, drug = name.split(
        "_",
        1
    )

    return protein, drug


###############################################################################
# STANDARD DOCKING PARSER
###############################################################################

def best_docking_scores(structure_dict):
    """
    Extract the best PLANTS score for each pocket.

    Used only in STANDARD CoBDock mode.
    """

    scores = {}

    for pocket_number, poses in structure_dict.items():

        best = None

        for pose in poses.values():

            score = pose.get(
                "total_score"
            )

            if score is None:
                continue

            try:
                score = float(score)

            except (
                TypeError,
                ValueError
            ):
                continue

            if (
                best is None
                or score < best
            ):

                best = score

        scores[int(pocket_number)] = best

    return scores


###############################################################################
# ALTERNATIVE VINA PARSER
###############################################################################

def get_vina_run(
    protein,
    drug,
    mutation,
    chain
):
    """
    Retrieve a docking run from merged_vina_results.json.

    Expected ID:

        KRAS_Adagrasib_KRAS-G12A_A
    """

    run_id = (
        f"{protein}_"
        f"{drug}_"
        f"{mutation}_"
        f"{chain}"
    )

    runs = vina_json.get(
        "docking_runs",
        {}
    )

    if run_id not in runs:

        return None

    return runs[run_id]


def get_best_vina_score(
    protein,
    drug,
    mutation,
    chain
):
    """
    Return the best (lowest) Vina energy for a docking run.
    """

    run = get_vina_run(
        protein,
        drug,
        mutation,
        chain
    )

    if run is None:

        return np.nan

    best_mode = (
        run
        .get("docking", {})
        .get("best_mode")
    )

    if best_mode is None:

        return np.nan

    energy = best_mode.get(
        "energy"
    )

    if energy is None:

        return np.nan

    return float(energy)


###############################################################################
# BUILD POCKET LIST — STANDARD
###############################################################################

def build_standard_pocket_list(
    pocket_list,
    docking_scores
):
    """
    Build Pocket objects from top_3_pockets.json.
    """

    pockets = []

    for pocket in pocket_list:

        rank = int(
            pocket["pocket_rank"]
        )

        pockets.append(

            Pocket(

                rank=rank,

                pocket_id=str(
                    pocket["selected_pocket_id"]
                ),

                pocket_score=float(
                    pocket["pocket_score"]
                ),

                center=np.array(
                    [
                        pocket["center_x"],
                        pocket["center_y"],
                        pocket["center_z"]
                    ],
                    dtype=float
                ),

                docking_score=docking_scores.get(
                    rank,
                    np.nan
                ),

                docking_score_type="PLANTS"
            )
        )

    pockets.sort(
        key=lambda x: x.rank
    )

    return pockets


###############################################################################
# BUILD POCKET LIST — FPOCKET
###############################################################################

def build_fpocket_list(
    structure_pockets,
    global_vina_score
):
    """
    Build Pocket objects from fpocket_data.json.

    fpocket_data.json contains all detected pockets.

    The TOP_N_FPOCKETS pockets are selected according to
    fpocket pocket_score, with higher scores ranked first.

    IMPORTANT:
    The Vina score is a GLOBAL docking score for the
    protein-drug run and is therefore NOT interpreted as
    a pocket-specific score.
    """

    candidates = []

    for pocket_id, pocket in structure_pockets.items():

        try:

            center = np.array(
                [
                    float(
                        pocket["center_x"]
                    ),
                    float(
                        pocket["center_y"]
                    ),
                    float(
                        pocket["center_z"]
                    )
                ],
                dtype=float
            )

            pocket_score = float(
                pocket.get(
                    "pocket_score",
                    np.nan
                )
            )

        except (
            KeyError,
            TypeError,
            ValueError
        ):

            continue

        candidates.append(
            {
                "pocket_id": str(
                    pocket_id
                ),

                "pocket_score": pocket_score,

                "center": center
            }
        )

    # -------------------------------------------------------------------------
    # Rank by fpocket score.
    #
    # Higher pocket_score = higher rank.
    # -------------------------------------------------------------------------

    candidates.sort(
        key=lambda x: x["pocket_score"],
        reverse=True
    )

    candidates = candidates[
        :TOP_N_FPOCKETS
    ]

    pockets = []

    for rank, pocket in enumerate(
        candidates,
        start=1
    ):

        pockets.append(

            Pocket(

                rank=rank,

                pocket_id=pocket["pocket_id"],

                pocket_score=pocket[
                    "pocket_score"
                ],

                center=pocket["center"],

                # We intentionally DO NOT assign the global
                # Vina score as a pocket docking score.
                docking_score=np.nan,

                docking_score_type="GLOBAL_VINA"
            )
        )

    return pockets


###############################################################################
# TRANSFORM CENTERS
###############################################################################

def transform_centers(
    pockets,
    rotation,
    translation
):
    """
    Apply structural alignment to mutant pocket centers.
    """

    transformed = []

    for pocket in pockets:

        xyz = (
            pocket.center @ rotation
        ) + translation

        transformed.append(

            Pocket(

                rank=pocket.rank,

                pocket_id=pocket.pocket_id,

                pocket_score=pocket.pocket_score,

                center=xyz,

                docking_score=pocket.docking_score,

                docking_score_type=(
                    pocket.docking_score_type
                )
            )
        )

    return transformed


###############################################################################
# DISTANCE MATRIX
###############################################################################

def build_distance_matrix(
    wt_pockets,
    mut_pockets
):

    wt_xyz = np.array(
        [
            p.center
            for p in wt_pockets
        ]
    )

    mut_xyz = np.array(
        [
            p.center
            for p in mut_pockets
        ]
    )

    return cdist(
        wt_xyz,
        mut_xyz,
        metric="euclidean"
    )


###############################################################################
# HUNGARIAN MATCHING
###############################################################################

def match_pockets(
    distance_matrix
):

    rows, cols = linear_sum_assignment(
        distance_matrix
    )

    matches = []

    for r, c in zip(
        rows,
        cols
    ):

        matches.append(
            (
                r,
                c,
                distance_matrix[r, c]
            )
        )

    return matches


###############################################################################
# MAIN ANALYSIS
###############################################################################

matched_rows = []

distance_rows = []

skipped_rows = []


def skip_entry(protein, drug, mutation, reason):
    """
    Record an entry that could not be processed instead of crashing.
    """

    print(
        f"    Skipping {protein} - {drug} - {mutation}: {reason}"
    )

    skipped_rows.append(
        {
            "Protein": protein,
            "Drug": drug,
            "Mutation": mutation,
            "Reason": reason,
        }
    )


###############################################################################
# STANDARD MODE
###############################################################################

if MODE == "standard":

    for ligand_name in sorted(
        pocket_json
    ):

        protein, drug = split_ligand_name(
            ligand_name
        )

        print(
            f"\nProcessing {protein} - {drug}"
        )

        pocket_entry = pocket_json[
            ligand_name
        ]

        if ligand_name not in pose_json:

            skip_entry(
                protein,
                drug,
                "ALL",
                "Missing in pose_data.json"
            )

            continue

        pose_entry = pose_json[
            ligand_name
        ]

        # ---------------------------------------------------------------------
        # WT
        # ---------------------------------------------------------------------

        try:

            wt_key = get_wildtype_key(
                pocket_entry
            )

        except ValueError:

            skip_entry(
                protein,
                drug,
                "ALL",
                "Wild-type missing in top_3_pockets.json"
            )

            continue

        wt_structure = next(
            iter(
                pocket_entry[
                    wt_key
                ]
            )
        )

        if wt_structure not in pose_entry.get(wt_key, {}):

            skip_entry(
                protein,
                drug,
                "ALL",
                f"Wild-type ({wt_key}/{wt_structure}) "
                f"missing in pose_data.json"
            )

            continue

        wt_docking = best_docking_scores(
            pose_entry[
                wt_key
            ][
                wt_structure
            ]
        )

        wt_pockets = (
            build_standard_pocket_list(
                pocket_entry[
                    wt_key
                ][
                    wt_structure
                ],

                wt_docking
            )
        )

        # ---------------------------------------------------------------------
        # Mutants
        # ---------------------------------------------------------------------

        for mutation in pocket_entry:

            if mutation == wt_key:

                continue

            alignment_key = (
                protein,
                mutation
            )

            if alignment_key not in alignments:

                print(
                    f"Alignment missing for "
                    f"{protein} {mutation}"
                )

                continue

            rotation = alignments[
                alignment_key
            ]["rotation"]

            translation = alignments[
                alignment_key
            ]["translation"]

            mut_structure = next(
                iter(
                    pocket_entry[
                        mutation
                    ]
                )
            )

            if mut_structure not in pose_entry.get(mutation, {}):

                skip_entry(
                    protein,
                    drug,
                    mutation,
                    f"Mutant ({mutation}/{mut_structure}) "
                    f"missing in pose_data.json"
                )

                continue

            mut_docking = best_docking_scores(
                pose_entry[
                    mutation
                ][
                    mut_structure
                ]
            )

            mut_pockets = (
                build_standard_pocket_list(
                    pocket_entry[
                        mutation
                    ][
                        mut_structure
                    ],

                    mut_docking
                )
            )

            mut_pockets = transform_centers(
                mut_pockets,
                rotation,
                translation
            )

            # -------------------------------------------------------------
            # Distance matrix
            # -------------------------------------------------------------

            D = build_distance_matrix(
                wt_pockets,
                mut_pockets
            )

            # -------------------------------------------------------------
            # Save complete distance matrix
            # -------------------------------------------------------------

            for i, wt in enumerate(
                wt_pockets
            ):

                for j, mut in enumerate(
                    mut_pockets
                ):

                    distance_rows.append(

                        {
                            "Protein": protein,
                            "Drug": drug,
                            "Mutation": mutation,

                            "WT_PocketRank": wt.rank,
                            "WT_PocketID": wt.pocket_id,

                            "Mut_PocketRank": mut.rank,
                            "Mut_PocketID": mut.pocket_id,

                            "Distance": float(
                                D[i, j]
                            )
                        }
                    )

            # -------------------------------------------------------------
            # Hungarian assignment
            # -------------------------------------------------------------

            matches = match_pockets(
                D
            )

            # -------------------------------------------------------------
            # Save matches
            # -------------------------------------------------------------

            for (
                wt_idx,
                mut_idx,
                dist
            ) in matches:

                wt = wt_pockets[
                    wt_idx
                ]

                mut = mut_pockets[
                    mut_idx
                ]

                matched_rows.append(

                    {
                        "Protein": protein,
                        "Drug": drug,
                        "Mutation": mutation,

                        "WT_PocketRank": wt.rank,
                        "WT_PocketID": wt.pocket_id,

                        "Mut_PocketRank": mut.rank,
                        "Mut_PocketID": mut.pocket_id,

                        "Distance": float(
                            dist
                        ),

                        "WT_PocketScore": (
                            wt.pocket_score
                        ),

                        "Mut_PocketScore": (
                            mut.pocket_score
                        ),

                        "WT_DockingScore": (
                            wt.docking_score
                        ),

                        "Mut_DockingScore": (
                            mut.docking_score
                        ),

                        "ScoreType": "PLANTS",

                        "WT_X": wt.center[0],
                        "WT_Y": wt.center[1],
                        "WT_Z": wt.center[2],

                        "Mut_X": mut.center[0],
                        "Mut_Y": mut.center[1],
                        "Mut_Z": mut.center[2],
                    }
                )


###############################################################################
# ALTERNATIVE FPOCKET + VINA MODE
###############################################################################

if MODE == "alternative":

    # ---------------------------------------------------------------------------
    # First pass:
    # Build a WT dictionary independently of the order of the Vina JSON.
    # ---------------------------------------------------------------------------

    wt_database = {}

    print()
    print("Building WT pocket database...")

    for run_id, run in sorted(
        vina_json.get("docking_runs", {}).items()
    ):

        metadata = run.get(
            "metadata",
            {}
        )

        protein = metadata.get("protein")
        drug = metadata.get("drug")
        mutation = metadata.get("mutation")
        chain = metadata.get("chain")

        if not all(
            [
                protein,
                drug,
                mutation,
                chain
            ]
        ):
            continue

        # Identify WT independently of ordering.

        if "WILD" not in mutation.upper():
            continue

        # -------------------------------------------------------------
        # Find corresponding fpocket structure
        # -------------------------------------------------------------

        if mutation not in pocket_json:

            print(
                f"WARNING: WT {mutation} not found in fpocket_data.json"
            )

            continue

        mutation_entry = pocket_json[
            mutation
        ]

        expected_structure = (
            f"{mutation}_{chain}"
        )

        if expected_structure in mutation_entry:

            structure_key = expected_structure

        elif mutation_entry:

            structure_key = next(
                iter(
                    mutation_entry
                )
            )

        else:

            print(
                f"WARNING: No structure found for WT "
                f"{protein} {drug} {mutation}"
            )

            continue

        structure_pockets = (
            mutation_entry[
                structure_key
            ]
        )

        # -------------------------------------------------------------
        # Build WT fpocket list
        # -------------------------------------------------------------

        wt_pockets = build_fpocket_list(
            structure_pockets,
            get_best_vina_score(
                protein,
                drug,
                mutation,
                chain
            )
        )

        # -------------------------------------------------------------
        # Store using protein + drug
        # -------------------------------------------------------------

        wt_key = (
            protein,
            drug
        )

        wt_database[wt_key] = {

            "protein": protein,

            "drug": drug,

            "mutation": mutation,

            "chain": chain,

            "structure": structure_key,

            "pockets": wt_pockets,

            "vina_score": get_best_vina_score(
                protein,
                drug,
                mutation,
                chain
            )
        }

        print(
            f"    WT loaded | "
            f"{protein} | "
            f"{drug} | "
            f"{mutation} | "
            f"{len(wt_pockets)} pockets | "
            f"Vina: "
            f"{wt_database[wt_key]['vina_score']}"
        )


    print(
        f"\nWT structures loaded: "
        f"{len(wt_database)}"
    )


    # ---------------------------------------------------------------------------
    # Second pass:
    # Process every MUTANT against its corresponding WT.
    # ---------------------------------------------------------------------------

    print()
    print("Matching mutant pockets against WT...")

    for run_id, run in sorted(
        vina_json.get("docking_runs", {}).items()
    ):

        metadata = run.get(
            "metadata",
            {}
        )

        protein = metadata.get("protein")
        drug = metadata.get("drug")
        mutation = metadata.get("mutation")
        chain = metadata.get("chain")

        if not all(
            [
                protein,
                drug,
                mutation,
                chain
            ]
        ):
            continue

        # -------------------------------------------------------------
        # Ignore WT entries
        # -------------------------------------------------------------

        if "WILD" in mutation.upper():
            continue

        # -------------------------------------------------------------
        # Find corresponding WT
        # -------------------------------------------------------------

        wt_key = (
            protein,
            drug
        )

        if wt_key not in wt_database:

            print(
                f"WARNING: WT not available for "
                f"{protein} {drug} {mutation}"
            )

            continue

        wt_data = wt_database[
            wt_key
        ]

        # -------------------------------------------------------------
        # Find mutant in fpocket data
        # -------------------------------------------------------------

        if mutation not in pocket_json:

            print(
                f"WARNING: fpocket data missing "
                f"{mutation}"
            )

            continue

        mutation_entry = pocket_json[
            mutation
        ]

        expected_structure = (
            f"{mutation}_{chain}"
        )

        if expected_structure in mutation_entry:

            structure_key = expected_structure

        elif mutation_entry:

            structure_key = next(
                iter(
                    mutation_entry
                )
            )

        else:

            print(
                f"WARNING: No structure found for "
                f"{protein} {drug} {mutation}"
            )

            continue

        structure_pockets = (
            mutation_entry[
                structure_key
            ]
        )

        # -------------------------------------------------------------
        # Build mutant pockets
        # -------------------------------------------------------------

        mut_vina_score = get_best_vina_score(
            protein,
            drug,
            mutation,
            chain
        )

        mut_pockets = build_fpocket_list(
            structure_pockets,
            mut_vina_score
        )

        # -------------------------------------------------------------
        # Alignment
        # -------------------------------------------------------------

        alignment_key = (
            protein,
            mutation
        )

        if alignment_key not in alignments:

            print(
                f"WARNING: Alignment missing for "
                f"{protein} {mutation}"
            )

            continue

        rotation = alignments[
            alignment_key
        ]["rotation"]

        translation = alignments[
            alignment_key
        ]["translation"]

        # -------------------------------------------------------------
        # Transform mutant pocket centers
        # -------------------------------------------------------------

        mut_pockets = transform_centers(
            mut_pockets,
            rotation,
            translation
        )

        wt_pockets = wt_data[
            "pockets"
        ]

        # -------------------------------------------------------------
        # Calculate distance matrix
        # -------------------------------------------------------------

        D = build_distance_matrix(
            wt_pockets,
            mut_pockets
        )

        # -------------------------------------------------------------
        # Save every WT-mutant pocket distance
        # -------------------------------------------------------------

        for i, wt in enumerate(
            wt_pockets
        ):

            for j, mut in enumerate(
                mut_pockets
            ):

                distance_rows.append(

                    {
                        "Protein": protein,

                        "Drug": drug,

                        "Mutation": mutation,

                        "WT_PocketRank": wt.rank,

                        "WT_PocketID": wt.pocket_id,

                        "Mut_PocketRank": mut.rank,

                        "Mut_PocketID": mut.pocket_id,

                        "Distance": float(
                            D[i, j]
                        )
                    }
                )

        # -------------------------------------------------------------
        # Hungarian matching
        # -------------------------------------------------------------

        matches = match_pockets(
            D
        )

        # -------------------------------------------------------------
        # WT and mutant global Vina scores
        # -------------------------------------------------------------

        wt_vina_score = wt_data[
            "vina_score"
        ]

        # -------------------------------------------------------------
        # Save matched pockets
        # -------------------------------------------------------------

        for (
            wt_idx,
            mut_idx,
            dist
        ) in matches:

            wt = wt_pockets[
                wt_idx
            ]

            mut = mut_pockets[
                mut_idx
            ]

            matched_rows.append(

                {
                    "Protein": protein,

                    "Drug": drug,

                    "Mutation": mutation,

                    "WT_PocketRank": wt.rank,

                    "WT_PocketID": wt.pocket_id,

                    "Mut_PocketRank": mut.rank,

                    "Mut_PocketID": mut.pocket_id,

                    "Distance": float(
                        dist
                    ),

                    "WT_PocketScore": (
                        wt.pocket_score
                    ),

                    "Mut_PocketScore": (
                        mut.pocket_score
                    ),

                    "WT_GlobalVinaScore": (
                        wt_vina_score
                    ),

                    "Mut_GlobalVinaScore": (
                        mut_vina_score
                    ),

                    "ScoreType": (
                        "GLOBAL_VINA"
                    ),

                    "WT_X": wt.center[0],

                    "WT_Y": wt.center[1],

                    "WT_Z": wt.center[2],

                    "Mut_X": mut.center[0],

                    "Mut_Y": mut.center[1],

                    "Mut_Z": mut.center[2],
                }
            )

        print(
            f"    OK | "
            f"{protein} | "
            f"{drug} | "
            f"{mutation} | "
            f"matched {len(matches)} pockets | "
            f"distance range: "
            f"{D.min():.2f}-{D.max():.2f} Å"
        )


###############################################################################
# CREATE DATAFRAMES
###############################################################################

matched_df = pd.DataFrame(
    matched_rows
)

distance_df = pd.DataFrame(
    distance_rows
)


###############################################################################
# HANDLE EMPTY OUTPUT
###############################################################################

if matched_df.empty:

    print(
        "\nWARNING: No matched pockets were generated."
    )

else:

    matched_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank"
        ],
        inplace=True
    )


if distance_df.empty:

    print(
        "\nWARNING: No distance matrix entries generated."
    )

else:

    distance_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank",
            "Mut_PocketRank"
        ],
        inplace=True
    )


###############################################################################
# SAVE
###############################################################################

matched_df.to_csv(
    MATCHED_OUTPUT,
    index=False
)

distance_df.to_csv(
    DISTANCE_OUTPUT,
    index=False
)

skipped_df = pd.DataFrame(
    skipped_rows,
    columns=["Protein", "Drug", "Mutation", "Reason"]
)

skipped_df.to_csv(
    SKIPPED_OUTPUT,
    index=False
)


###############################################################################
# SUMMARY
###############################################################################

print()
print("=" * 70)
print("Finished.")
print("=" * 70)

print(
    f"Input mode       : {MODE}"
)

print(
    f"Matched pockets  : {len(matched_df)}"
)

print(
    f"Distance entries  : {len(distance_df)}"
)

print(
    f"Skipped entries  : {len(skipped_df)}"
)

print()
print(
    f"Matched output:"
)

print(
    f"    {MATCHED_OUTPUT}"
)

print(
    f"Distance output:"
)

print(
    f"    {DISTANCE_OUTPUT}"
)

print(
    f"Skipped entries log:"
)

print(
    f"    {SKIPPED_OUTPUT}"
)

print("=" * 70)

dist = matched_df["Distance"]

print("\nPocket Matching Statistics")
print("---------------------------")
print(f"Number of matches: {len(dist)}")
print(f"Mean distance:     {dist.mean():.2f} Å")
print(f"Median distance:   {dist.median():.2f} Å")
print(f"Std distance:      {dist.std():.2f} Å")
print(f"Min distance:      {dist.min():.2f} Å")
print(f"Max distance:      {dist.max():.2f} Å")

print(f"≤1 Å : {(dist <= 1).mean()*100:.2f}%")
print(f"≤2 Å : {(dist <= 2).mean()*100:.2f}%")
print(f"≤3 Å : {(dist <= 3).mean()*100:.2f}%")
print(f"≤5 Å : {(dist <= 5).mean()*100:.2f}%")
print(f"≤5 Å : {(dist <= 6).mean()*100:.2f}%")
