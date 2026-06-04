"""Generate and cache the faithful +/-5% non-lactating PhD BED simulation bank.

This is the only script in the parameter-BED workflow that runs ODE
simulations. It checkpoints results so interrupted runs can be resumed.
"""

from __future__ import annotations

import argparse
import json
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


ORIGINAL_SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose"]
EXTENDED_SPECIES = ORIGINAL_SPECIES + ["Glucagon"]
DAYS = list(range(54, 90))
PRIOR_HALF_RANGE = 0.05
RELATIVE_NOISE = 0.10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate the cached +/-5% PhD BED ODE bank.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-samples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--checkpoint-every", type=int, default=100)
    parser.add_argument(
        "--species",
        nargs="+",
        default=ORIGINAL_SPECIES,
        help="Biomarkers stored for BED. Add Glucagon for the extended analysis.",
    )
    parser.add_argument(
        "--admissibility-species",
        nargs="+",
        default=ORIGINAL_SPECIES,
        help="States used for biological filtering; defaults to the original eight-marker PhD filter.",
    )
    parser.add_argument("--method", default="BDF")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-9)
    return parser.parse_args()


def extract(simulation, species: list[str]) -> dict[str, float]:
    return {
        f"{state}_day_{day}": float(np.interp(day, simulation.t, simulation.y[:, STATE_INDEX[state]]))
        for day in DAYS
        for state in species
    }


def write_checkpoint(
    output_dir: Path,
    parameters: pd.DataFrame,
    outputs: list[dict[str, float]],
    admissibility: list[dict[str, object]],
) -> None:
    parameters.iloc[: len(outputs)].to_csv(output_dir / "prior_parameter_samples.csv", index=False)
    pd.DataFrame(outputs).to_csv(output_dir / "ode_output_features.csv", index=False)
    pd.DataFrame(admissibility).to_csv(output_dir / "admissibility.csv", index=False)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    unknown = sorted(set(args.species + args.admissibility_species).difference(STATE_INDEX))
    if unknown:
        raise ValueError(f"Unknown species: {unknown}")
    nominal = default_parameters()
    names = [parameter.name for parameter in PARAMETERS]
    rng = np.random.default_rng(args.seed)
    multipliers = rng.uniform(
        1.0 - PRIOR_HALF_RANGE,
        1.0 + PRIOR_HALF_RANGE,
        size=(args.n_samples, len(names)),
    )
    parameters = pd.DataFrame(
        {name: nominal[name] * multipliers[:, index] for index, name in enumerate(names)}
    )

    # Match the original MATLAB non-lactating BED temporal grid: STEP = 1/6 day.
    scenario = constant_non_lactating(days=100.0, step=1.0 / 6.0)
    reference = run_simulation(
        scenario, parameters=nominal, method=args.method, rtol=args.rtol, atol=args.atol
    )
    nominal_output = extract(reference, args.species)
    pd.DataFrame([nominal_output]).to_csv(args.output_dir / "nominal_output.csv", index=False)
    pd.DataFrame(
        {
            "biomarker": args.species,
            "abs_error": [
                RELATIVE_NOISE * float(np.mean(reference.y[:, STATE_INDEX[state]]))
                for state in args.species
            ],
        }
    ).to_csv(args.output_dir / "measurement_noise.csv", index=False)

    output_path = args.output_dir / "ode_output_features.csv"
    admissibility_path = args.output_dir / "admissibility.csv"
    outputs = pd.read_csv(output_path).to_dict("records") if output_path.exists() else []
    admissibility = (
        pd.read_csv(admissibility_path).to_dict("records") if admissibility_path.exists() else []
    )
    if len(outputs) != len(admissibility):
        raise ValueError("Existing output and admissibility checkpoints have different lengths.")
    if outputs and set(outputs[0]) != set(nominal_output):
        raise ValueError(
            "Existing checkpoint biomarker columns do not match --species. "
            "Use a new output directory for a different biomarker panel."
        )
    start = len(outputs)
    if start > args.n_samples:
        raise ValueError(
            f"Existing checkpoint has {start} simulations, exceeding requested {args.n_samples}."
        )
    thresholds = AdmissibilityThresholds()

    for index in range(start, args.n_samples):
        params = dict(nominal)
        params.update(parameters.iloc[index].to_dict())
        try:
            simulation = run_simulation(
                scenario, parameters=params, method=args.method, rtol=args.rtol, atol=args.atol
            )
            outputs.append(extract(simulation, args.species))
            _, summary = trajectory_admissibility(
                reference, simulation, states=args.admissibility_species, thresholds=thresholds
            )
            admissibility.append(
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
            outputs.append({key: np.nan for key in nominal_output})
            admissibility.append(
                {
                    "admissible": False,
                    "penalty": np.inf,
                    "min_correlation": np.nan,
                    "max_average_difference": np.nan,
                    "max_norm_difference": np.nan,
                    "worst_species": f"simulation_failed: {exc}",
                }
            )

        completed = index + 1
        if completed % args.checkpoint_every == 0 or completed == args.n_samples:
            write_checkpoint(args.output_dir, parameters, outputs, admissibility)
            print(f"Checkpointed {completed}/{args.n_samples} simulations")

    metadata = {
        "workflow": "non-lactating PhD BED bank with configurable biomarker panel",
        "n_samples": args.n_samples,
        "prior": "independent uniform nominal +/-5%",
        "prior_half_range": PRIOR_HALF_RANGE,
        "relative_measurement_noise": RELATIVE_NOISE,
        "species": args.species,
        "admissibility_species": args.admissibility_species,
        "candidate_days": DAYS,
        "admissible_rows": int(pd.DataFrame(admissibility)["admissible"].astype(bool).sum()),
    }
    (args.output_dir / "bank_metadata.json").write_text(json.dumps(metadata, indent=2))
    print(f"Saved faithful PhD BED simulation bank to {args.output_dir}")


if __name__ == "__main__":
    main()
