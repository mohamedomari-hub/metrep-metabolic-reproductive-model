"""Run configurable MetRep model scenarios.

Examples
--------
List built-in scenarios:
    python model_analysis/08_run_model_scenarios.py --list

Run all built-in non-Dexa scenarios:
    python model_analysis/08_run_model_scenarios.py --scenario all

Run an exact forcing schedule from CSV:
    python model_analysis/08_run_model_scenarios.py --forcing-csv data/my_scenario.csv --mode lactating

The forcing CSV must contain columns:
    time_days,DMI,Milk
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))
sys.path.insert(0, str(PROJECT_ROOT))

from model_definition.scenarios import (
    available_scenarios,
    built_in_scenario,
    custom_step_feeding,
    scenario_default_days,
    scenario_from_csv,
)
from model_definition.simulate import run_simulation, save_result


DEFAULT_PLOT_STATES = ["Glucose", "Insulin", "IGF1", "P4", "E2", "Follicle", "CL"]
LACTATING_C0_SCENARIOS = {
    "lactating_c0_20",
    "lactating_c0_22_5",
    "lactating_c0_25",
    "lactating_c0_30",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Simulate built-in or CSV-defined BovSys feeding/lactation scenarios."
    )
    parser.add_argument("--list", action="store_true", help="List built-in scenarios and exit.")
    parser.add_argument(
        "--scenario",
        nargs="+",
        default=["baseline_non_lactating"],
        help="Built-in scenario name(s), or 'all'. Ignored when --forcing-csv is used.",
    )
    parser.add_argument("--forcing-csv", type=Path, help="CSV with time_days,DMI,Milk columns.")
    parser.add_argument("--name", help="Output scenario name for --forcing-csv or --step-feeding.")
    parser.add_argument("--mode", choices=["non_lactating", "lactating"], default="lactating")
    parser.add_argument("--c0", type=float, help="Carbohydrate availability coefficient.")
    parser.add_argument(
        "--days",
        type=float,
        default=None,
        help="Override simulation horizon in days. By default, each built-in scenario uses its MetRep horizon.",
    )
    parser.add_argument("--dt", type=float, default=1.0 / 48.0, help="Output spacing in days.")
    parser.add_argument(
        "--states",
        nargs="+",
        default=DEFAULT_PLOT_STATES,
        help="State names to plot.",
    )
    parser.add_argument(
        "--step-feeding",
        action="store_true",
        help="Build a custom DMI step-change scenario from the options below.",
    )
    parser.add_argument("--baseline-dmi", type=float, default=11700.0, help="Baseline DMI, g/day.")
    parser.add_argument("--baseline-milk", type=float, default=0.0, help="Baseline milk, L/day.")
    parser.add_argument("--dmi-multiplier", type=float, default=1.0, help="DMI multiplier during step.")
    parser.add_argument("--start-day", type=float, default=46.0, help="Step start day.")
    parser.add_argument("--end-day", type=float, help="Step end day. Omit for permanent step.")
    return parser.parse_args()


def build_scenarios(args: argparse.Namespace):
    if args.forcing_csv:
        return [
            scenario_from_csv(
                str(args.forcing_csv),
                step=args.dt,
                mode=args.mode,
                c0=args.c0,
                name=args.name,
            )
        ]

    if args.step_feeding:
        return [
            custom_step_feeding(
                days=90.0 if args.days is None else args.days,
                step=args.dt,
                mode=args.mode,
                baseline_dmi_g_per_day=args.baseline_dmi,
                baseline_milk_l_per_day=args.baseline_milk,
                dmi_multiplier=args.dmi_multiplier,
                start_day=args.start_day,
                end_day=args.end_day,
                c0=args.c0,
                name=args.name or "custom_step_feeding",
            )
        ]

    scenario_names = available_scenarios() if "all" in args.scenario else args.scenario
    return [built_in_scenario(name, days=args.days, step=args.dt) for name in scenario_names]


def main() -> None:
    args = parse_args()

    if args.list:
        print("Built-in scenarios:")
        for name in available_scenarios():
            print(f"- {name} ({scenario_default_days(name):g} days)")
        return

    from model_analysis.plotting import plot_all_states_grid, plot_scenario_comparison, plot_selected_states

    simulation_dir = PROJECT_ROOT / "results" / "simulations"
    figure_dir = PROJECT_ROOT / "results" / "figures"
    table_dir = PROJECT_ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    scenarios = build_scenarios(args)
    results = []
    metadata = []

    for scenario in scenarios:
        result = run_simulation(scenario)
        csv_path, npz_path = save_result(result, simulation_dir)
        result.scenario_name = scenario.name
        fig_path = plot_selected_states(
            result,
            figure_dir / f"{scenario.name}_selected_states.png",
            states=args.states,
        )
        all_states_path = None
        if scenario.mode == "non_lactating":
            all_states_path = plot_all_states_grid(
                result,
                figure_dir / f"{scenario.name}_all_species.png",
            )
        results.append(result)
        metadata.append(
            {
                "scenario": scenario.name,
                "mode": scenario.mode,
                "c0": scenario.c0,
                "days": float(scenario.t_eval[-1]),
                "points": len(scenario.t_eval),
                "dmi_min": float(scenario.dmi.min()),
                "dmi_max": float(scenario.dmi.max()),
                "milk_min": float(scenario.milk.min()),
                "milk_max": float(scenario.milk.max()),
                "description": scenario.description,
                "csv": str(csv_path.relative_to(PROJECT_ROOT)),
                "npz": str(npz_path.relative_to(PROJECT_ROOT)),
                "figure": str(fig_path.relative_to(PROJECT_ROOT)),
                "all_species_figure": "" if all_states_path is None else str(all_states_path.relative_to(PROJECT_ROOT)),
            }
        )
        print(f"Saved {scenario.name}: {csv_path}")
        print(f"Saved {scenario.name}: {npz_path}")
        print(f"Saved {scenario.name}: {fig_path}")
        if all_states_path is not None:
            print(f"Saved {scenario.name}: {all_states_path}")

    metadata_path = table_dir / "scenario_run_metadata.csv"
    pd.DataFrame(metadata).to_csv(metadata_path, index=False)
    print(f"Saved metadata: {metadata_path}")

    lactating_results = [result for result in results if result.scenario_name in LACTATING_C0_SCENARIOS]
    if len(lactating_results) > 1:
        comparison_path = plot_scenario_comparison(
            lactating_results,
            figure_dir / "lactating_c0_overlay.png",
            states=args.states,
            title="Lactating c0 Scenarios",
        )
        print(f"Saved lactating c0 overlay figure: {comparison_path}")


if __name__ == "__main__":
    main()
