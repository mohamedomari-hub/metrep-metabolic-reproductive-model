"""Prepare CSV inputs for surrogate-assisted BED from Python ODE simulations.

This script creates the input tables consumed by ``surrogate_bed_pipeline.py``:

- prior_parameter_samples.csv
- ode_output_features.csv
- candidate_map.csv
- nominal_output.csv
- admissibility.csv

It runs the translated MetRep Python ODE model. Start with a modest sample
count for testing, then increase ``--n-samples`` for a reportable analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = PROJECT_ROOT / "MetRep_Python"
sys.path.insert(0, str(PYTHON_ROOT))

from metrep.admissibility import AdmissibilityThresholds, trajectory_admissibility
from metrep.initial_conditions import STATE_INDEX
from metrep.parameters import PARAMETERS, default_parameters
from metrep.scenarios import constant_non_lactating
from metrep.simulate import run_simulation


DEFAULT_SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]
DEFAULT_DAYS = list(range(55, 91))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate surrogate BED input CSVs from MetRep ODE simulations.")
    parser.add_argument("--output-dir", type=Path, default=Path("analyses/bayesian_experimental_design/surrogate_bed/input_tables"))
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--days", type=float, default=90.0)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--sample-low", type=float, default=0.995, help="Lower multiplicative prior bound.")
    parser.add_argument("--sample-high", type=float, default=1.005, help="Upper multiplicative prior bound.")
    parser.add_argument("--parameter-names", nargs="+", default=[parameter.name for parameter in PARAMETERS])
    parser.add_argument("--species", nargs="+", default=DEFAULT_SPECIES)
    parser.add_argument("--sample-days", nargs="+", type=float, default=DEFAULT_DAYS)
    parser.add_argument(
        "--candidate-mode",
        choices=["single-day", "cumulative", "both"],
        default="single-day",
        help=(
            "Candidate-map style. 'single-day' matches the initial pilot; "
            "'cumulative' stacks all sampled days up to each candidate day; "
            "'both' writes both maps."
        ),
    )
    parser.add_argument("--method", default="BDF")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-9)
    return parser.parse_args()


def sample_parameters(args: argparse.Namespace) -> pd.DataFrame:
    nominal = default_parameters()
    unknown = sorted(set(args.parameter_names).difference(nominal))
    if unknown:
        raise ValueError(f"Unknown parameter names: {unknown}")
    rng = np.random.default_rng(args.seed)
    multipliers = rng.uniform(args.sample_low, args.sample_high, size=(args.n_samples, len(args.parameter_names)))
    rows = []
    for sample_index in range(args.n_samples):
        row = {
            name: float(nominal[name] * multipliers[sample_index, parameter_index])
            for parameter_index, name in enumerate(args.parameter_names)
        }
        rows.append(row)
    return pd.DataFrame(rows)


def extract_output_features(simulation, species: list[str], sample_days: list[float]) -> dict[str, float]:
    features: dict[str, float] = {}
    for day in sample_days:
        day_label = f"{day:g}".replace(".", "p")
        for state in species:
            if state not in STATE_INDEX:
                raise ValueError(f"Unknown species/state name: {state}")
            values = simulation.y[:, STATE_INDEX[state]]
            value = float(np.interp(float(day), simulation.t, values))
            features[f"{state}_day_{day_label}"] = value
    return features


def build_candidate_map(species: list[str], sample_days: list[float]) -> pd.DataFrame:
    rows = []
    for day in sample_days:
        day_label = f"{day:g}".replace(".", "p")
        day_columns = [f"{state}_day_{day_label}" for state in species]
        rows.append({"candidate": f"all_species_day_{day_label}", "columns": "|".join(day_columns)})
        for state in species:
            column = f"{state}_day_{day_label}"
            rows.append({"candidate": column, "columns": column})
    return pd.DataFrame(rows)


def build_cumulative_candidate_map(species: list[str], sample_days: list[float]) -> pd.DataFrame:
    """Build thesis-style cumulative measurement candidates.

    A candidate ending at day ``d`` contains every selected sampling day up to
    and including ``d``. For example, ``all_species_cumulative_to_day_68``
    stacks all species measured from the first sampled day through day 68.
    """

    rows = []
    ordered_days = sorted(float(day) for day in sample_days)
    subset_definitions = [
        ("PGF_E2_FSH_INH", ["PGF", "E2", "FSH", "INH"]),
        ("PGF_E2_FSH", ["PGF", "E2", "FSH"]),
        ("PGF_E2", ["PGF", "E2"]),
    ]
    for day in ordered_days:
        day_label = f"{day:g}".replace(".", "p")
        cumulative_days = [past_day for past_day in ordered_days if past_day <= day]
        all_columns = [
            f"{state}_day_{f'{past_day:g}'.replace('.', 'p')}"
            for past_day in cumulative_days
            for state in species
        ]
        rows.append(
            {
                "candidate": f"all_species_cumulative_to_day_{day_label}",
                "columns": "|".join(all_columns),
            }
        )
        for state in species:
            columns = [f"{state}_day_{f'{past_day:g}'.replace('.', 'p')}" for past_day in cumulative_days]
            rows.append({"candidate": f"{state}_cumulative_to_day_{day_label}", "columns": "|".join(columns)})
        for label, subset_species in subset_definitions:
            available_species = [state for state in subset_species if state in species]
            if not available_species:
                continue
            columns = [
                f"{state}_day_{f'{past_day:g}'.replace('.', 'p')}"
                for past_day in cumulative_days
                for state in available_species
            ]
            rows.append({"candidate": f"{label}_cumulative_to_day_{day_label}", "columns": "|".join(columns)})
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    scenario = constant_non_lactating(days=args.days, step=args.dt)
    nominal_params = default_parameters()
    reference = run_simulation(scenario, parameters=nominal_params, method=args.method, rtol=args.rtol, atol=args.atol)
    nominal_output = extract_output_features(reference, args.species, args.sample_days)

    parameter_samples = sample_parameters(args)
    output_rows = []
    admissibility_rows = []
    thresholds = AdmissibilityThresholds()

    for sample_number, row in enumerate(parameter_samples.itertuples(index=False), start=1):
        params = dict(nominal_params)
        params.update({name: float(value) for name, value in zip(parameter_samples.columns, row)})
        try:
            simulation = run_simulation(scenario, parameters=params, method=args.method, rtol=args.rtol, atol=args.atol)
            output_rows.append(extract_output_features(simulation, args.species, args.sample_days))
            _, summary = trajectory_admissibility(reference, simulation, states=args.species, thresholds=thresholds)
            admissibility_rows.append(
                {
                    "admissible": bool(summary["admissible"]),
                    "penalty": float(summary["penalty"]),
                    "min_correlation": float(summary["min_correlation"]),
                    "max_average_difference": float(summary["max_average_difference"]),
                    "max_norm_difference": float(summary["max_norm_difference"]),
                    "worst_species": str(summary["worst_species"]),
                }
            )
        except Exception as exc:
            output_rows.append({key: np.nan for key in nominal_output})
            admissibility_rows.append(
                {
                    "admissible": False,
                    "penalty": np.inf,
                    "min_correlation": np.nan,
                    "max_average_difference": np.nan,
                    "max_norm_difference": np.nan,
                    "worst_species": f"simulation_failed: {exc}",
                }
            )

        if sample_number == 1 or sample_number % 25 == 0 or sample_number == args.n_samples:
            print(f"Simulated {sample_number}/{args.n_samples} parameter samples")

    outputs = pd.DataFrame(output_rows)
    admissibility = pd.DataFrame(admissibility_rows)
    finite_rows = outputs.notna().all(axis=1)
    if not finite_rows.all():
        parameter_samples = parameter_samples.loc[finite_rows].reset_index(drop=True)
        outputs = outputs.loc[finite_rows].reset_index(drop=True)
        admissibility = admissibility.loc[finite_rows].reset_index(drop=True)

    parameter_samples.to_csv(args.output_dir / "prior_parameter_samples.csv", index=False)
    outputs.to_csv(args.output_dir / "ode_output_features.csv", index=False)
    if args.candidate_mode in {"single-day", "both"}:
        build_candidate_map(args.species, args.sample_days).to_csv(args.output_dir / "candidate_map.csv", index=False)
    if args.candidate_mode in {"cumulative", "both"}:
        build_cumulative_candidate_map(args.species, args.sample_days).to_csv(
            args.output_dir / "candidate_map_cumulative.csv",
            index=False,
        )
    pd.DataFrame([nominal_output]).to_csv(args.output_dir / "nominal_output.csv", index=False)
    admissibility.to_csv(args.output_dir / "admissibility.csv", index=False)

    print(f"Saved surrogate BED input tables to {args.output_dir}")
    print(f"Rows retained: {len(parameter_samples)}")
    print(f"Admissible rows: {int(admissibility['admissible'].sum())}/{len(admissibility)}")


if __name__ == "__main__":
    main()
