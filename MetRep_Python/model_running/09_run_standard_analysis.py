"""Run 90-day scenario simulations plus all-parameter standard-scenario analyses.

Run from project root:
    python model_running/09_run_standard_analysis.py
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import os
from pathlib import Path
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))
sys.path.insert(0, str(PROJECT_ROOT))

from model_definition.analysis import (
    local_sensitivity,
    structural_identifiability_svd,
    uncertainty_trajectories,
)
from model_definition.parameters import PARAMETERS
from model_definition.plotting import (
    plot_compensation_network,
    plot_all_states_grid,
    plot_scenario_comparison,
    plot_selected_states,
    plot_uncertainty_band,
)
from model_definition.scenarios import available_scenarios, built_in_scenario
from model_definition.simulate import run_simulation, save_result


DEFAULT_OUTPUTS = ["Glucose", "Insulin", "IGF1", "P4", "E2", "Follicle", "CL"]
LACTATING_C0_SCENARIOS = {
    "lactating_c0_20",
    "lactating_c0_22_5",
    "lactating_c0_25",
    "lactating_c0_30",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scenario",
        nargs="+",
        default=["all"],
        help="Built-in scenario name(s) to simulate for 90 days, or 'all'.",
    )
    parser.add_argument(
        "--analysis-scenario",
        default="baseline_non_lactating",
        help="Built-in scenario used for sensitivity, identifiability, and uncertainty.",
    )
    parser.add_argument(
        "--days",
        type=float,
        default=None,
        help="Override scenario simulation horizon in days. Default: each scenario's MetRep horizon.",
    )
    parser.add_argument("--analysis-days", type=float, default=90.0, help="Analysis horizon for the standard scenario.")
    parser.add_argument("--dt", type=float, default=0.25, help="Output spacing in days.")
    parser.add_argument("--prefix", default="metrep", help="Output filename prefix.")
    parser.add_argument("--outputs", nargs="+", default=DEFAULT_OUTPUTS, help="Observable states.")
    parser.add_argument("--rel-step", type=float, default=1e-3, help="Identifiability relative step.")
    parser.add_argument("--sensitivity-step", type=float, default=0.01)
    parser.add_argument("--uncertainty-samples", type=int, default=40)
    parser.add_argument("--uncertainty-state", default="Glucose")
    parser.add_argument("--network-all-nodes", action="store_true", help="Show all analyzed parameters in the compensation network.")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def scenario_names_from_args(args: argparse.Namespace) -> list[str]:
    if "all" in args.scenario:
        return available_scenarios()
    return args.scenario


def run_one_scenario(
    scenario_name: str,
    args: argparse.Namespace,
    figure_dir: Path,
    simulation_dir: Path,
):
    scenario = built_in_scenario(scenario_name, days=args.days, step=args.dt)
    prefix = f"{args.prefix}_{scenario.name}"

    baseline_scenario = replace(scenario, name=f"{prefix}_baseline")
    baseline = run_simulation(baseline_scenario)
    csv_path, npz_path = save_result(baseline, simulation_dir)
    baseline.scenario_name = scenario.name
    baseline_figure = plot_selected_states(
        baseline,
        figure_dir / f"{prefix}_baseline_selected_states.png",
        states=args.outputs,
    )
    all_states_figure = None
    if scenario.mode == "non_lactating":
        all_states_figure = plot_all_states_grid(
            baseline,
            figure_dir / f"{prefix}_all_species.png",
        )

    print(f"Simulated scenario: {scenario.name}")
    print(f"  duration: {float(scenario.t_eval[-1]):g} days, dt={args.dt:g} days")
    print(f"  CSV: {csv_path}")
    print(f"  NPZ: {npz_path}")
    print(f"  figure: {baseline_figure}")
    if all_states_figure is not None:
        print(f"  all species figure: {all_states_figure}")
    print()

    return baseline


def run_standard_analyses(
    args: argparse.Namespace,
    parameter_names: list[str],
    table_dir: Path,
    figure_dir: Path,
    simulation_dir: Path,
) -> None:
    scenario = built_in_scenario(args.analysis_scenario, days=args.analysis_days, step=args.dt)
    prefix = f"{args.prefix}_{scenario.name}"

    sensitivity = local_sensitivity(
        scenario,
        parameter_names,
        args.outputs,
        relative_step=args.sensitivity_step,
    )
    sensitivity_path = table_dir / f"{prefix}_local_sensitivity_allparams.csv"
    sensitivity.to_csv(sensitivity_path, index=False)

    identifiability = structural_identifiability_svd(
        scenario,
        parameter_names=parameter_names,
        outputs=args.outputs,
        relative_step=args.rel_step,
        tolerance_value=1e-8,
        ranking_metric="rel2_colnorm",
    )
    ident_prefix = f"{prefix}_structid_svd_allparams"
    np.savez_compressed(
        simulation_dir / f"{ident_prefix}.npz",
        params=np.array(parameter_names, dtype=object),
        outputs=np.array(args.outputs, dtype=object),
        t=scenario.t_eval,
        nominal_outputs=identifiability["nominal_outputs"],
        sensitivity_matrix=identifiability["sensitivity_matrix"],
        singular_values=identifiability["singular_values"],
        vt=identifiability["vt"],
        rank=np.array([identifiability["rank"]]),
        nullity=np.array([identifiability["nullity"]]),
        threshold=np.array([identifiability["threshold"]]),
        ranking_scores=identifiability["ranking_scores"],
        nullspace_score=identifiability["nullspace_score"],
        identifiable_score=identifiability["identifiable_score"],
    )
    singular_values_path = table_dir / f"{ident_prefix}_singular_values.csv"
    np.savetxt(
        singular_values_path,
        identifiability["singular_values"],
        delimiter=",",
        header="singular_value",
        comments="",
    )
    holistic_path = table_dir / f"{ident_prefix}_holistic_table.csv"
    identifiability["holistic_table"].to_csv(holistic_path, index=False)
    compensation_edges_path = table_dir / f"{ident_prefix}_compensation_edges.csv"
    identifiability["compensation_edges"].to_csv(compensation_edges_path, index=False)
    network_path = plot_compensation_network(
        identifiability["compensation_edges"],
        identifiability["nullspace"],
        parameter_names,
        figure_dir / f"{ident_prefix}_compensation_network.png",
        include_all_nodes=args.network_all_nodes,
    )

    t, samples = uncertainty_trajectories(
        scenario,
        parameter_names,
        state=args.uncertainty_state,
        n_samples=args.uncertainty_samples,
        relative_bounds=(0.95, 1.05),
        seed=args.seed,
    )
    uncertainty_path = simulation_dir / f"{prefix}_uncertainty_{args.uncertainty_state.lower()}_allparams.npz"
    np.savez_compressed(
        uncertainty_path,
        t=t,
        samples=samples,
        state=args.uncertainty_state,
        parameters=np.array(parameter_names, dtype=object),
    )
    uncertainty_figure = plot_uncertainty_band(
        t,
        samples,
        args.uncertainty_state,
        figure_dir / f"{prefix}_uncertainty_{args.uncertainty_state.lower()}_allparams.png",
    )

    print(f"Analysis scenario: {scenario.name}")
    print(f"  duration: {args.analysis_days:g} days, dt={args.dt:g} days")
    print(f"  parameters analyzed: {len(parameter_names)}")
    print(f"  sensitivity table: {sensitivity_path}")
    print(f"  identifiability NPZ: {simulation_dir / f'{ident_prefix}.npz'}")
    print(f"  singular values: {singular_values_path}")
    print(f"  holistic table: {holistic_path}")
    print(f"  compensation edges: {compensation_edges_path}")
    if network_path is not None:
        print(f"  compensation network: {network_path}")
    print(f"  uncertainty ensemble: {uncertainty_path}")
    print(f"  uncertainty figure: {uncertainty_figure}")
    print()


def main() -> None:
    args = parse_args()

    table_dir = PROJECT_ROOT / "results" / "tables"
    figure_dir = PROJECT_ROOT / "results" / "figures"
    simulation_dir = PROJECT_ROOT / "results" / "simulations"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    simulation_dir.mkdir(parents=True, exist_ok=True)

    parameter_names = [parameter.name for parameter in PARAMETERS]
    baselines = [
        run_one_scenario(name, args, figure_dir, simulation_dir)
        for name in scenario_names_from_args(args)
    ]
    lactating_results = [result for result in baselines if result.scenario_name in LACTATING_C0_SCENARIOS]
    if len(lactating_results) > 1:
        lactating_path = plot_scenario_comparison(
            lactating_results,
            figure_dir / f"{args.prefix}_lactating_c0_overlay.png",
            states=args.outputs,
            title="Lactating c0 Scenarios",
        )
        print(f"Lactating c0 overlay figure: {lactating_path}")

    run_standard_analyses(args, parameter_names, table_dir, figure_dir, simulation_dir)


if __name__ == "__main__":
    main()
