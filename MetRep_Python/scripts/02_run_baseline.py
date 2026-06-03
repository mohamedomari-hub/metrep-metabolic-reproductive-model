"""Run the baseline non-Dexa metabolic-reproductive simulation."""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from metrep.plotting import plot_all_states_grid, plot_selected_states
from metrep.scenarios import constant_non_lactating
from metrep.simulate import run_simulation, save_result


def main() -> None:
    scenario = constant_non_lactating(days=60.0)
    result = run_simulation(scenario)
    csv_path, npz_path = save_result(result, PROJECT_ROOT / "results" / "simulations")
    fig_path = plot_selected_states(
        result,
        PROJECT_ROOT / "results" / "figures" / "baseline_selected_states.png",
    )
    all_states_path = plot_all_states_grid(
        result,
        PROJECT_ROOT / "results" / "figures" / "baseline_all_states.png",
    )
    print(f"Saved simulation: {csv_path}")
    print(f"Saved simulation: {npz_path}")
    print(f"Saved figure: {fig_path}")
    print(f"Saved figure: {all_states_path}")


if __name__ == "__main__":
    main()
