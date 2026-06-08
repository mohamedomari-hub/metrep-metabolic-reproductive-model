"""Run all-parameter uncertainty propagation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from model_analysis.analysis import uncertainty_trajectories
from model_definition.parameters import PARAMETERS
from model_definition.plotting import plot_uncertainty_band
from model_definition.scenarios import constant_non_lactating


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--days", type=float, default=90.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--state", default="Glucose")
    parser.add_argument("--prefix", default="uncertainty_allparams")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scenario = constant_non_lactating(days=args.days, step=args.dt)
    parameter_names = [parameter.name for parameter in PARAMETERS]
    state = args.state
    t, samples = uncertainty_trajectories(
        scenario,
        parameter_names,
        state=state,
        n_samples=args.samples,
        relative_bounds=(0.95, 1.05),
        seed=args.seed,
    )

    sim_path = PROJECT_ROOT / "results" / "simulations" / f"{args.prefix}_{state.lower()}.npz"
    sim_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        sim_path,
        t=t,
        samples=samples,
        state=state,
        parameters=np.array(parameter_names, dtype=object),
    )
    fig_path = plot_uncertainty_band(
        t,
        samples,
        state,
        PROJECT_ROOT / "results" / "figures" / f"{args.prefix}_{state.lower()}.png",
    )
    print(f"Saved uncertainty ensemble: {sim_path}")
    print(f"Saved uncertainty figure: {fig_path}")


if __name__ == "__main__":
    main()
