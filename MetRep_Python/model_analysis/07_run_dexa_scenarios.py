"""Run optional Dexa perturbation scenarios.

The default Python model runners remain non-Dexa. This script explicitly uses
the 25-state Dexa extension translated from the MATLAB v3 Dexa files.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt

from model_definition.dexa import DexaConfig, run_dexa_simulation
from model_analysis.plotting import _apply_plot_style
from model_definition.scenarios import available_scenarios, built_in_scenario
from model_definition.simulate import run_simulation, save_result


DEFAULT_STATES = ["Glucose", "Insulin", "Glucagon", "P4", "IGF1", "A_dep", "A_cent", "C_e"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the optional Dexa PK/PD extension.")
    parser.add_argument(
        "--scenario",
        default="baseline_non_lactating",
        choices=available_scenarios(),
        help="Built-in scenario to perturb with Dexa.",
    )
    parser.add_argument("--days", type=float, default=None, help="Simulation horizon in days.")
    parser.add_argument("--dt", type=float, default=1.0 / 48.0, help="Output spacing in days.")
    parser.add_argument("--dose-day", type=float, default=0.0, help="Day of IM Dexa dose.")
    parser.add_argument(
        "--dose-mg-per-kg",
        type=float,
        default=0.02,
        help="IM Dexa dose in mg/kg; MATLAB reference default is 0.02.",
    )
    parser.add_argument(
        "--body-weight-kg",
        type=float,
        default=600.0,
        help="Body weight used to convert mg/kg dose to ng.",
    )
    parser.add_argument(
        "--states",
        nargs="+",
        default=DEFAULT_STATES,
        help="State names to plot. Dexa states are A_dep, A_cent, and C_e.",
    )
    parser.add_argument(
        "--simulation-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "simulations",
        help="Directory for optional trajectory CSV/NPZ files.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "figures",
        help="Directory for Dexa summary figures.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Directory for summary tables.",
    )
    parser.add_argument(
        "--prefix",
        help="Output prefix. Default is '<scenario>_dexa_day_<dose-day>'.",
    )
    parser.add_argument(
        "--save-trajectories",
        action="store_true",
        help="Save full Dexa and no-Dexa trajectories as CSV/NPZ.",
    )
    return parser.parse_args()


def _state_series(result, state: str) -> np.ndarray | None:
    if state not in result.state_names:
        return None
    return result.y[:, result.state_names.index(state)]


def plot_dexa_summary(dexa_result, core_result, states: list[str], dose_day: float, output_path: Path) -> Path:
    """Plot Dexa response against the matching no-Dexa baseline."""

    _apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(states),
        1,
        figsize=(12, max(5.0, 2.25 * len(states))),
        sharex=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)

    for ax, state in zip(axes, states):
        dexa_values = _state_series(dexa_result, state)
        if dexa_values is None:
            ax.axis("off")
            continue
        core_values = _state_series(core_result, state)
        if core_values is not None:
            ax.plot(core_result.t, core_values, color="0.35", linestyle="--", linewidth=1.7, label="No Dexa")
        ax.plot(dexa_result.t, dexa_values, color="tab:blue", linewidth=2.0, label="Dexa")
        ax.axvline(dose_day, color="tab:red", linestyle=":", linewidth=1.4, label="Dose day")
        ax.set_ylabel(state.replace("_", " "))
        ax.grid(alpha=0.25, linewidth=0.8)
        ax.margins(x=0.01)

    axes[0].legend(loc="upper right", frameon=True)
    axes[-1].set_xlabel("Time (days)")
    fig.suptitle("Dexa perturbation response compared with no-Dexa baseline", y=1.01)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def summarize_response(dexa_result, core_result, states: list[str], dose_day: float) -> pd.DataFrame:
    """Create a compact response table for plotted states."""

    rows = []
    for state in states:
        dexa_values = _state_series(dexa_result, state)
        if dexa_values is None:
            continue
        core_values = _state_series(core_result, state)
        if core_values is None:
            core_interp = np.full_like(dexa_values, np.nan)
            no_dexa_min = np.nan
            no_dexa_max = np.nan
            max_abs_difference = np.nan
            peak_idx = int(np.nanargmax(dexa_values))
        else:
            core_interp = np.interp(dexa_result.t, core_result.t, core_values)
            delta = dexa_values - core_interp
            no_dexa_min = float(np.nanmin(core_interp))
            no_dexa_max = float(np.nanmax(core_interp))
            max_abs_difference = float(np.nanmax(np.abs(delta)))
            peak_idx = int(np.nanargmax(np.abs(delta))) if np.any(np.isfinite(delta)) else 0
        rows.append(
            {
                "state": state,
                "dose_day": dose_day,
                "dexa_min": float(np.nanmin(dexa_values)),
                "dexa_max": float(np.nanmax(dexa_values)),
                "no_dexa_min": no_dexa_min,
                "no_dexa_max": no_dexa_max,
                "max_abs_difference": max_abs_difference,
                "time_of_max_abs_difference": float(dexa_result.t[peak_idx]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    scenario = built_in_scenario(args.scenario, days=args.days, step=args.dt)
    prefix = args.prefix or f"{scenario.name}_dexa_day_{args.dose_day:g}".replace(".", "p")

    config = DexaConfig(
        enabled=True,
        dose_mg_per_kg=args.dose_mg_per_kg,
        body_weight_kg=args.body_weight_kg,
        dose_day=args.dose_day,
    )
    core_result = run_simulation(scenario)
    dexa_result = run_dexa_simulation(scenario, config=config)

    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    figure_path = plot_dexa_summary(
        dexa_result,
        core_result,
        args.states,
        args.dose_day,
        args.figure_dir / f"{prefix}_summary.png",
    )
    summary_path = args.table_dir / f"{prefix}_response_summary.csv"
    summarize_response(dexa_result, core_result, args.states, args.dose_day).to_csv(summary_path, index=False)

    if args.save_trajectories:
        save_result(dexa_result, args.simulation_dir)
        save_result(core_result, args.simulation_dir)

    print(f"Scenario: {scenario.name}")
    print(f"Dose: {config.dose_mg_per_kg:g} mg/kg at day {config.dose_day:g}")
    print(f"Dose amount: {config.dose_ng:g} ng")
    print(f"Saved figure: {figure_path}")
    print(f"Saved response table: {summary_path}")
    if args.save_trajectories:
        print(f"Saved full trajectories in: {args.simulation_dir}")


if __name__ == "__main__":
    main()
