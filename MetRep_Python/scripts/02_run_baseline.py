"""Run the baseline non-Dexa metabolic-reproductive simulation."""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from bovsys.plotting import plot_selected_states
from bovsys.scenarios import constant_non_lactating
from bovsys.simulate import run_simulation, save_result


def main() -> None:
    scenario = constant_non_lactating(days=60.0)
    result = run_simulation(scenario)
    csv_path, npz_path = save_result(result, PROJECT_ROOT / "results" / "simulations")
    fig_path = plot_selected_states(
        result,
        PROJECT_ROOT / "results" / "figures" / "baseline_selected_states.png",
    )
    print(f"Saved simulation: {csv_path}")
    print(f"Saved simulation: {npz_path}")
    print(f"Saved figure: {fig_path}")


if __name__ == "__main__":
    main()
