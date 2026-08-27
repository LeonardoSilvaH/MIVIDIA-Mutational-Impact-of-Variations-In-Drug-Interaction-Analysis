#!/usr/bin/env python3

"""
merge_vina_jsons.py

Recursively finds AutoDock Vina *.log.json files and converts them into
standardized JSON files while preserving docking-mode information and
source provenance.

Expected input structure (approximately):

    ROOT/
    ├── Adagrasib/
    │   └── docking/
    │       └── vina/
    │           └── KRAS_Adagrasib/
    │               └── KRAS-G12A/
    │                   └── KRAS-G12A_A/
    │                       └── KRAS_Adagrasib_KRAS-G12A_A.log.json

The resulting output will be:

    OUTPUT/
    ├── KRAS_Adagrasib_KRAS-G12A_A.json
    ├── KRAS_Adagrasib_KRAS-G12C_A.json
    └── ...

    OUTPUT/
    └── vina_merge_summary.csv


Usage:

    python merge_vina_jsons.py /path/to/results

Optional:

    python merge_vina_jsons.py /path/to/results \
        --output merged_vina

    python merge_vina_jsons.py /path/to/results \
        --output merged_vina \
        --overwrite
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from datetime import datetime


# -------------------------------------------------------------------------
# Filename parsing
# -------------------------------------------------------------------------

# Expected filename:
#
# KRAS_Adagrasib_KRAS-G12A_A.log.json
#
# We deliberately use the final components rather than assuming that
# protein/drug names contain no underscores.
#
# Pattern:
#   <protein>_<drug>_<mutation>_<chain>.log.json
#
# Mutation examples:
#   KRAS-G12A
#   BRAF-V600E
#   RAF1-S257L
#   MYC-P74A
#
FILENAME_PATTERN = re.compile(
    r"^(?P<prefix>.+)_(?P<mutation>[A-Za-z0-9]+-[A-Za-z0-9]+)_(?P<chain>[A-Za-z0-9]+)\.log\.json$"
)


# -------------------------------------------------------------------------
# Metadata extraction
# -------------------------------------------------------------------------

def extract_metadata(json_path: Path):
    """
    Extract protein, drug, mutation, and chain from the filename/path.

    Example:

        KRAS_Adagrasib_KRAS-G12A_A.log.json

    becomes:

        protein  = KRAS
        drug     = Adagrasib
        mutation = KRAS-G12A
        chain    = A

    The prefix before the mutation is split using the expected protein/drug
    directory structure when possible.
    """

    filename = json_path.name

    match = FILENAME_PATTERN.match(filename)

    if not match:
        raise ValueError(
            f"Filename does not match expected pattern: {filename}"
        )

    prefix = match.group("prefix")
    mutation = match.group("mutation")
    chain = match.group("chain")

    # -------------------------------------------------------------
    # Prefer information from the directory structure.
    #
    # Example:
    #
    # .../vina/KRAS_Adagrasib/
    #     KRAS-G12A/
    #     KRAS-G12A_A/
    #
    # -------------------------------------------------------------

    parts = json_path.parts

    protein = None
    drug = None

    # Find the directory immediately containing the mutation directory.
    #
    # We search backwards because the structure is expected to contain:
    #
    #   <protein>_<drug>/<mutation>/<mutation>_<chain>/
    #
    mutation_dir = None

    for part in reversed(parts):
        if part == mutation:
            mutation_dir = part
            break

    if mutation_dir is not None:
        try:
            mutation_index = parts.index(mutation_dir)

            if mutation_index >= 1:
                protein_drug_dir = parts[mutation_index - 1]

                # Expected:
                # KRAS_Adagrasib
                #
                # Use the first underscore as separator.
                if "_" in protein_drug_dir:
                    protein, drug = protein_drug_dir.split("_", 1)

        except (ValueError, IndexError):
            pass

    # -------------------------------------------------------------
    # Fallback: derive protein/drug from filename prefix.
    # -------------------------------------------------------------

    if protein is None or drug is None:

        # Try to identify the protein from the mutation.
        mutation_protein = mutation.split("-", 1)[0]

        if prefix.startswith(mutation_protein + "_"):
            remainder = prefix[len(mutation_protein) + 1:]

            protein = mutation_protein
            drug = remainder

        else:
            # Last-resort interpretation:
            # first underscore-separated component = protein
            # remaining components = drug
            prefix_parts = prefix.split("_")

            if len(prefix_parts) >= 2:
                protein = prefix_parts[0]
                drug = "_".join(prefix_parts[1:])
            else:
                protein = prefix
                drug = "UNKNOWN"

    return {
        "protein": protein,
        "drug": drug,
        "mutation": mutation,
        "chain": chain,
    }


# -------------------------------------------------------------------------
# JSON validation
# -------------------------------------------------------------------------

def load_vina_json(json_path: Path):
    """
    Load and validate a Vina log JSON.

    Expected structure:

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

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("Top-level JSON object is not a dictionary.")

    if not data:
        raise ValueError("JSON contains no docking modes.")

    valid_modes = {}

    for key, value in data.items():

        if not isinstance(value, dict):
            continue

        # Require at least energy and mode.
        if "energy" not in value:
            continue

        if "mode" not in value:
            continue

        try:
            energy = float(value["energy"])
        except (TypeError, ValueError):
            continue

        try:
            mode = int(value["mode"])
        except (TypeError, ValueError):
            continue

        # Preserve all available information.
        mode_data = {
            "mode": mode,
            "energy": energy,
        }

        # Optional fields.
        if "rmsd_lb" in value:
            mode_data["rmsd_lb"] = value["rmsd_lb"]

        if "rmsd_ub" in value:
            mode_data["rmsd_ub"] = value["rmsd_ub"]

        if "pKd" in value:
            mode_data["pKd"] = value["pKd"]

        valid_modes[str(mode)] = mode_data

    if not valid_modes:
        raise ValueError(
            "No valid docking modes containing 'mode' and 'energy' were found."
        )

    return valid_modes


# -------------------------------------------------------------------------
# Find input JSON files
# -------------------------------------------------------------------------

def find_vina_jsons(root_dir: Path):
    """
    Recursively find all *.log.json files.
    """

    files = sorted(root_dir.rglob("*.log.json"))

    return files


# -------------------------------------------------------------------------
# Determine best docking mode
# -------------------------------------------------------------------------

def find_best_mode(modes):
    """
    Return the mode with the lowest (most negative) Vina energy.
    """

    best_mode = min(
        modes.values(),
        key=lambda x: x["energy"]
    )

    return best_mode


# -------------------------------------------------------------------------
# Create standardized record
# -------------------------------------------------------------------------

def build_record(json_path: Path, root_dir: Path):
    """
    Build the standardized output representation.
    """

    metadata = extract_metadata(json_path)

    modes = load_vina_json(json_path)

    best_mode = find_best_mode(modes)

    try:
        relative_path = json_path.relative_to(root_dir)
    except ValueError:
        relative_path = json_path

    record = {
        "schema_version": "1.0",

        "metadata": {
            "protein": metadata["protein"],
            "drug": metadata["drug"],
            "mutation": metadata["mutation"],
            "chain": metadata["chain"],
        },

        "source": {
            "file": str(relative_path),
            "absolute_path": str(json_path.resolve()),
            "merged_at": datetime.now().isoformat(timespec="seconds"),
        },

        "docking": {
            "engine": "AutoDock Vina",
            "number_of_modes": len(modes),

            "best_mode": {
                "mode": best_mode["mode"],
                "energy": best_mode["energy"],
            },

            "modes": modes,
        }
    }

    # Preserve pKd of best mode when available.
    if "pKd" in best_mode:
        record["docking"]["best_mode"]["pKd"] = best_mode["pKd"]

    return record


# -------------------------------------------------------------------------
# Output filename
# -------------------------------------------------------------------------

def output_filename(record):
    """
    Generate standardized output filename.

    Example:

        KRAS_Adagrasib_KRAS-G12A_A.json
    """

    metadata = record["metadata"]

    return (
        f"{metadata['protein']}_"
        f"{metadata['drug']}_"
        f"{metadata['mutation']}_"
        f"{metadata['chain']}.json"
    )


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Recursively find AutoDock Vina *.log.json files and "
            "create standardized provenance-preserving JSON files."
        )
    )

    parser.add_argument(
        "root",
        type=Path,
        help="Root directory containing the docking results."
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: <root>/merged_vina"
        )
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing merged JSON files."
    )

    args = parser.parse_args()

    root_dir = args.root.resolve()

    if not root_dir.exists():
        print(f"ERROR: Root directory does not exist: {root_dir}")
        sys.exit(1)

    if not root_dir.is_dir():
        print(f"ERROR: Root path is not a directory: {root_dir}")
        sys.exit(1)

    if args.output is None:
        output_dir = root_dir / "merged_vina"
    else:
        output_dir = args.output.resolve()

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    print("=" * 70)
    print("Vina JSON Merge")
    print("=" * 70)

    print(f"Input directory : {root_dir}")
    print(f"Output directory: {output_dir}")
    print()

    # -------------------------------------------------------------
    # Find files
    # -------------------------------------------------------------

    json_files = find_vina_jsons(root_dir)

    # Don't accidentally process files we previously generated.
    json_files = [
        path
        for path in json_files
        if output_dir not in path.parents
    ]

    print(f"Found {len(json_files)} *.log.json files.")
    print()

    if not json_files:
        print("No Vina JSON files found.")
        sys.exit(0)

    # -------------------------------------------------------------
    # Process files
    # -------------------------------------------------------------

    summary = []

    successful = 0
    failed = 0
    skipped = 0

    for index, json_path in enumerate(json_files, start=1):

        print(
            f"[{index}/{len(json_files)}] "
            f"{json_path.name}"
        )

        try:

            record = build_record(
                json_path=json_path,
                root_dir=root_dir
            )

            output_name = output_filename(record)
            output_path = output_dir / output_name

            if output_path.exists() and not args.overwrite:

                print(
                    f"    SKIPPED: output already exists: "
                    f"{output_path.name}"
                )

                skipped += 1

                summary.append({
                    "status": "skipped",
                    "protein": record["metadata"]["protein"],
                    "drug": record["metadata"]["drug"],
                    "mutation": record["metadata"]["mutation"],
                    "chain": record["metadata"]["chain"],
                    "source_file": record["source"]["file"],
                    "output_file": output_name,
                    "number_of_modes": record["docking"]["number_of_modes"],
                    "best_mode": record["docking"]["best_mode"]["mode"],
                    "best_energy": record["docking"]["best_mode"]["energy"],
                })

                continue

            # Write JSON.
            with open(
                output_path,
                "w",
                encoding="utf-8"
            ) as f:

                json.dump(
                    record,
                    f,
                    indent=4
                )

                f.write("\n")

            best = record["docking"]["best_mode"]

            print(
                f"    OK"
                f" | {record['metadata']['protein']}"
                f" | {record['metadata']['drug']}"
                f" | {record['metadata']['mutation']}"
                f" | chain {record['metadata']['chain']}"
            )

            print(
                f"    Modes: {record['docking']['number_of_modes']}"
                f" | Best mode: {best['mode']}"
                f" | Energy: {best['energy']}"
            )

            successful += 1

            summary.append({
                "status": "merged",
                "protein": record["metadata"]["protein"],
                "drug": record["metadata"]["drug"],
                "mutation": record["metadata"]["mutation"],
                "chain": record["metadata"]["chain"],
                "source_file": record["source"]["file"],
                "output_file": output_name,
                "number_of_modes": record["docking"]["number_of_modes"],
                "best_mode": best["mode"],
                "best_energy": best["energy"],
            })

        except Exception as e:

            print(
                f"    ERROR: {e}"
            )

            failed += 1

            summary.append({
                "status": "error",
                "protein": "",
                "drug": "",
                "mutation": "",
                "chain": "",
                "source_file": str(json_path.relative_to(root_dir)),
                "output_file": "",
                "number_of_modes": "",
                "best_mode": "",
                "best_energy": "",
            })

    # -------------------------------------------------------------
    # Summary CSV
    # -------------------------------------------------------------

    summary_path = output_dir / "vina_merge_summary.csv"

    fieldnames = [
        "status",
        "protein",
        "drug",
        "mutation",
        "chain",
        "source_file",
        "output_file",
        "number_of_modes",
        "best_mode",
        "best_energy",
    ]

    with open(
        summary_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )

        writer.writeheader()
        writer.writerows(summary)

    # -------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------

    print()
    print("=" * 70)
    print("MERGE COMPLETE")
    print("=" * 70)

    print(f"Found      : {len(json_files)}")
    print(f"Merged     : {successful}")
    print(f"Skipped    : {skipped}")
    print(f"Errors     : {failed}")
    print()
    print(f"Output     : {output_dir}")
    print(f"Summary CSV: {summary_path}")
    print("=" * 70)

    if failed > 0:
        print()
        print(
            "WARNING: Some files could not be processed. "
            "Check vina_merge_summary.csv."
        )


if __name__ == "__main__":
    main()