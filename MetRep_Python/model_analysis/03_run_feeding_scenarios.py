"""Run MetRep non-lactating feeding scenarios with their default horizons."""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from model_analysis.plotting import plot_all_states_grid
from model_definition.scenarios import feeding_scenarios
from model_definition.simulate import run_simulation, save_result


def main() -> None:
    for scenario in feeding_scenarios():
        result = run_simulation(scenario)
        save_result(result, PROJECT_ROOT / "results" / "simulations")
        fig_path = plot_all_states_grid(
            result,
            PROJECT_ROOT / "results" / "figures" / f"{scenario.name}_all_species.png",
        )
        print(f"Saved all-species figure: {fig_path}")


if __name__ == "__main__":
    main()
