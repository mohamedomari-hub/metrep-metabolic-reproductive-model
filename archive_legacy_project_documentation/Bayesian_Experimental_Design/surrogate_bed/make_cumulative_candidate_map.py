"""Create thesis-style cumulative BED candidate maps from existing outputs.

Use this after ``prepare_surrogate_inputs.py`` has already generated
``ode_output_features.csv``. It does not rerun the ODE model.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from prepare_surrogate_inputs import DEFAULT_SPECIES, build_cumulative_candidate_map  # noqa: E402


DAY_PATTERN = re.compile(r"^(?P<species>.+)_day_(?P<day>[0-9]+(?:p[0-9]+)?)$")


def parse_args() -> argparse.Namespace:
    base = Path("Project_Documentation/Bayesian_Experimental_Design/surrogate_bed/input_tables_large")
    parser = argparse.ArgumentParser(description="Build a cumulative BED candidate map from existing output columns.")
    parser.add_argument("--outputs-csv", type=Path, default=base / "ode_output_features.csv")
    parser.add_argument("--output-csv", type=Path, default=base / "candidate_map_cumulative.csv")
    parser.add_argument("--species", nargs="+", default=DEFAULT_SPECIES)
    parser.add_argument("--start-day", type=float, help="Optional first day to include.")
    parser.add_argument("--end-day", type=float, help="Optional last day to include.")
    return parser.parse_args()


def infer_days(outputs: pd.DataFrame, species: list[str]) -> list[float]:
    days = set()
    for column in outputs.columns:
        match = DAY_PATTERN.match(column)
        if not match:
            continue
        if match.group("species") not in species:
            continue
        days.add(float(match.group("day").replace("p", ".")))
    if not days:
        raise ValueError("No species_day_* columns found in the output table.")
    return sorted(days)


def main() -> None:
    args = parse_args()
    outputs = pd.read_csv(args.outputs_csv, nrows=1)
    days = infer_days(outputs, args.species)
    if args.start_day is not None:
        days = [day for day in days if day >= args.start_day]
    if args.end_day is not None:
        days = [day for day in days if day <= args.end_day]
    if not days:
        raise ValueError("No sampling days remain after applying day filters.")

    candidate_map = build_cumulative_candidate_map(args.species, days)
    output_columns = set(outputs.columns)
    unknown_columns = sorted(
        {
            column
            for columns in candidate_map["columns"]
            for column in str(columns).split("|")
            if column not in output_columns
        }
    )
    if unknown_columns:
        raise ValueError(f"Cumulative candidate map uses unknown output columns: {unknown_columns[:10]}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    candidate_map.to_csv(args.output_csv, index=False)
    print(f"Saved cumulative candidate map: {args.output_csv}")
    print(f"Candidates: {len(candidate_map)}")
    print(f"Days: {days[0]:g}..{days[-1]:g} ({len(days)} days)")


if __name__ == "__main__":
    main()
