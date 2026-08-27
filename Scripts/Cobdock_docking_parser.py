#!/usr/bin/env python3

"""
merge_vina_jsons.py

Recursively finds all AutoDock Vina *.log.json files and merges them
into ONE standardized JSON file while preserving provenance.

Example input:

ROOT/
├── Adagrasib/
│   └── docking/
│       └── vina/
│           └── KRAS_Adagrasib/
│               ├── KRAS-G12A/
│               │   └── KRAS-G12A_A/
│               │       └── KRAS_Adagrasib_KRAS-G12A_A.log.json
│               ├── KRAS-G12C/
│               │   └── ...
│               └── ...
│
└── Daraxonrasib/
    └── docking/
        └── vina/
            └── ...

Output:

ROOT/
└── merged_vina_results.json

The resulting JSON contains ALL docking runs and ALL docking modes.
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime


# -------------------------------------------------------------------------
# Metadata extraction
# -------------------------------------------------------------------------

def extract_metadata(json_path: Path):
    """
    Extract protein, drug, mutation, and chain from the directory structure.

    Expected structure:

        <protein>_<drug>/
            <mutation>/
                <mutation>_<chain>/
                    file.log.json

    Example:

        KRAS_Adagrasib/
            KRAS-G12A/
                KRAS-G12A_A/
                    KRAS_Adagrasib_KRAS-G12A_A.log.json

    Also handles:

        KRAS-WILD-TYE/
            KRAS-WILD-TYE_A/
                KRAS_Daraxonrasib_KRAS-WILD-TYE_A.log.json
    """

    filename = json_path.name

    if not filename.endswith(".log.json"):
        raise ValueError(
            f"File does not end with .log.json: {filename}"
        )

    # -------------------------------------------------------------
    # Mutation
    # -------------------------------------------------------------

    # Example:
    #
    # parent:
    #     KRAS-G12A_A
    #
    # parent.parent:
    #     KRAS-G12A

    mutation_dir = json_path.parent.parent

    mutation = mutation_dir.name

    if not mutation:
        raise ValueError(
            f"Could not determine mutation from: {json_path}"
        )

    # -------------------------------------------------------------
    # Protein + drug
    # -------------------------------------------------------------

    # Example:
    #
    # KRAS-G12A/
    #     ↑
    # mutation_dir
    #
    # KRAS_Adagrasib/
    # ↑
    # mutation_dir.parent

    protein_drug_dir = mutation_dir.parent.name

    if "_" not in protein_drug_dir:
        raise ValueError(
            f"Could not determine protein/drug from directory: "
            f"{protein_drug_dir}"
        )

    protein, drug = protein_drug_dir.split("_", 1)

    # -------------------------------------------------------------
    # Chain
    # -------------------------------------------------------------

    # Example:
    #
    # KRAS-G12A_A
    #
    # Split from the final underscore.

    pose_directory = json_path.parent.name

    if "_" not in pose_directory:
        raise ValueError(
            f"Could not determine chain from directory: "
            f"{pose_directory}"
        )

    chain = pose_directory.rsplit("_", 1)[1]

    return {
        "protein": protein,
        "drug": drug,
        "mutation": mutation,
        "chain": chain,
    }


# -------------------------------------------------------------------------
# Load Vina JSON
# -------------------------------------------------------------------------

def load_vina_json(json_path: Path):
    """
    Load the original Vina JSON and preserve all docking modes.

    Expected format:

        {
            "1": {
                "mode": 1,
                "energy": -11.8,
                "rmsd_lb": 0.0,
                "rmsd_ub": 0.0,
                "pKd": 8.65
            },
            ...
        }
    """

    with open(
        json_path,
        "r",
        encoding="utf-8"
    ) as f:

        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError(
            "JSON top-level object is not a dictionary."
        )

    valid_modes = {}

    for key, value in data.items():

        if not isinstance(value, dict):
            continue

        if "mode" not in value:
            continue

        if "energy" not in value:
            continue

        try:
            mode = int(value["mode"])
            energy = float(value["energy"])

        except (TypeError, ValueError):
            continue

        mode_data = {
            "mode": mode,
            "energy": energy
        }

        # Preserve optional values.

        if "rmsd_lb" in value:
            mode_data["rmsd_lb"] = value["rmsd_lb"]

        if "rmsd_ub" in value:
            mode_data["rmsd_ub"] = value["rmsd_ub"]

        if "pKd" in value:
            mode_data["pKd"] = value["pKd"]

        valid_modes[str(mode)] = mode_data

    if not valid_modes:
        raise ValueError(
            "No valid docking modes found."
        )

    return valid_modes


# -------------------------------------------------------------------------
# Find best docking mode
# -------------------------------------------------------------------------

def find_best_mode(modes):
    """
    Select the mode with the lowest (most negative) Vina energy.
    """

    return min(
        modes.values(),
        key=lambda mode: mode["energy"]
    )


# -------------------------------------------------------------------------
# Build docking-run record
# -------------------------------------------------------------------------

def build_record(
    json_path: Path,
    root_dir: Path
):

    metadata = extract_metadata(json_path)

    modes = load_vina_json(json_path)

    best_mode = find_best_mode(modes)

    try:
        relative_path = str(
            json_path.relative_to(root_dir)
        )

    except ValueError:
        relative_path = str(json_path)

    # Unique identifier for this docking run.

    run_id = (
        f"{metadata['protein']}_"
        f"{metadata['drug']}_"
        f"{metadata['mutation']}_"
        f"{metadata['chain']}"
    )

    record = {

        "metadata": {
            "protein": metadata["protein"],
            "drug": metadata["drug"],
            "mutation": metadata["mutation"],
            "chain": metadata["chain"]
        },

        "source": {
            "file": relative_path,
            "absolute_path": str(
                json_path.resolve()
            )
        },

        "docking": {

            "number_of_modes": len(modes),

            "best_mode": {
                "mode": best_mode["mode"],
                "energy": best_mode["energy"]
            },

            "modes": modes
        }
    }

    # Preserve pKd for best mode.

    if "pKd" in best_mode:

        record["docking"]["best_mode"]["pKd"] = (
            best_mode["pKd"]
        )

    return run_id, record


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Merge all AutoDock Vina *.log.json files "
            "into one standardized JSON."
        )
    )

    parser.add_argument(
        "root",
        type=Path,
        help="Root directory containing docking results."
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output JSON filename. "
            "Default: <root>/merged_vina_results.json"
        )
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output JSON."
    )

    args = parser.parse_args()

    # -------------------------------------------------------------
    # Validate root
    # -------------------------------------------------------------

    root_dir = args.root.resolve()

    if not root_dir.exists():

        print(
            f"ERROR: Directory does not exist:\n"
            f"{root_dir}"
        )

        sys.exit(1)

    if not root_dir.is_dir():

        print(
            f"ERROR: Path is not a directory:\n"
            f"{root_dir}"
        )

        sys.exit(1)

    # -------------------------------------------------------------
    # Output
    # -------------------------------------------------------------

    if args.output is None:

        output_path = (
            root_dir /
            "merged_vina_results.json"
        )

    else:

        output_path = args.output.resolve()

    # -------------------------------------------------------------
    # Find JSON files
    # -------------------------------------------------------------

    json_files = sorted(
        root_dir.rglob("*.log.json")
    )

    # Do not accidentally include JSONs inside an output directory
    # if the user places the output somewhere below root.

    json_files = [
        path
        for path in json_files
        if path.resolve() != output_path.resolve()
    ]

    print("=" * 70)
    print("AutoDock Vina JSON Merger")
    print("=" * 70)

    print(f"Root directory : {root_dir}")
    print(f"Output file    : {output_path}")
    print()

    print(
        f"Found {len(json_files)} .log.json files."
    )

    print()

    if not json_files:

        print(
            "ERROR: No *.log.json files were found."
        )

        sys.exit(1)

    # -------------------------------------------------------------
    # Check existing output
    # -------------------------------------------------------------

    if output_path.exists() and not args.overwrite:

        print(
            f"ERROR: Output file already exists:\n"
            f"{output_path}\n\n"
            f"Use --overwrite if you want to replace it."
        )

        sys.exit(1)

    # -------------------------------------------------------------
    # Initialize merged structure
    # -------------------------------------------------------------

    merged = {

        "schema_version": "1.0",

        "docking_engine": "AutoDock Vina",

        "created_at": datetime.now().isoformat(
            timespec="seconds"
        ),

        "source_directory": str(root_dir),

        "number_of_docking_runs": 0,

        "docking_runs": {}

    }

    # -------------------------------------------------------------
    # Statistics
    # -------------------------------------------------------------

    successful = 0
    errors = 0

    total_modes = 0

    # -------------------------------------------------------------
    # Process files
    # -------------------------------------------------------------

    for index, json_path in enumerate(
        json_files,
        start=1
    ):

        print(
            f"[{index}/{len(json_files)}] "
            f"{json_path.name}"
        )

        try:

            run_id, record = build_record(
                json_path=json_path,
                root_dir=root_dir
            )

            # -----------------------------------------------------
            # Detect duplicate run IDs
            # -----------------------------------------------------

            if run_id in merged["docking_runs"]:

                existing = merged["docking_runs"][run_id]

                print(
                    f"    ERROR: Duplicate docking run ID:\n"
                    f"    {run_id}\n"
                    f"    Existing:\n"
                    f"    {existing['source']['file']}\n"
                    f"    New:\n"
                    f"    {json_path}"
                )

                errors += 1

                continue

            # -----------------------------------------------------
            # Store record
            # -----------------------------------------------------

            merged["docking_runs"][run_id] = record

            number_modes = (
                record["docking"]["number_of_modes"]
            )

            best_energy = (
                record["docking"]["best_mode"]["energy"]
            )

            total_modes += number_modes

            successful += 1

            print(
                f"    OK"
                f" | {record['metadata']['protein']}"
                f" | {record['metadata']['drug']}"
                f" | {record['metadata']['mutation']}"
                f" | chain {record['metadata']['chain']}"
            )

            print(
                f"    Modes: {number_modes}"
                f" | Best mode: "
                f"{record['docking']['best_mode']['mode']}"
                f" | Energy: {best_energy}"
            )

        except Exception as e:

            errors += 1

            print(
                f"    ERROR: {e}"
            )

    # -------------------------------------------------------------
    # Update summary
    # -------------------------------------------------------------

    merged["number_of_docking_runs"] = successful

    merged["summary"] = {

        "files_found": len(json_files),

        "files_merged": successful,

        "files_with_errors": errors,

        "total_docking_modes": total_modes

    }

    # -------------------------------------------------------------
    # Write single merged JSON
    # -------------------------------------------------------------

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            merged,
            f,
            indent=4
        )

        f.write("\n")

    # -------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------

    print()
    print("=" * 70)
    print("MERGE COMPLETE")
    print("=" * 70)

    print(
        f"Files found       : {len(json_files)}"
    )

    print(
        f"Files merged      : {successful}"
    )

    print(
        f"Files with errors : {errors}"
    )

    print(
        f"Total Vina modes  : {total_modes}"
    )

    print()
    print(
        f"Merged JSON:"
    )

    print(
        f"    {output_path}"
    )

    print("=" * 70)

    if errors > 0:

        print()
        print(
            "WARNING: Some files could not be merged."
        )


if __name__ == "__main__":
    main()