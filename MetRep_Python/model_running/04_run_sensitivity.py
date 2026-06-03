"""Run all-parameter local sensitivity analysis for the 22-state core model.

This script deliberately analyzes ``model_definition.parameters.PARAMETERS``:
the 98 metabolic-reproductive parameters from the MATLAB model. Dexa PK/PD
constants are not included here, so the sensitivity workflow remains
reproducible for the non-Dexa baseline model.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from model_definition.analysis import local_sensitivity
from model_definition.parameters import PARAMETERS
from model_definition.scenarios import constant_non_lactating


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run local sensitivity for the 22-state non-Dexa MetRep core. "
            "Dexa PK/PD constants are intentionally excluded."
        )
    )
    parser.add_argument("--days", type=float, default=90.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--rel-step", type=float, default=0.01)
    parser.add_argument("--prefix", default="local_sensitivity_allparams")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scenario = constant_non_lactating(days=args.days, step=args.dt)
    # Keep this analysis on the core parameter set only. Dexa constants live in
    # DexaConfig and are used only by model_running/07_run_dexa_scenarios.py.
    parameters = [parameter.name for parameter in PARAMETERS]
    outputs = ["Glucose", "Insulin", "IGF1", "P4", "E2", "Follicle", "CL"]
    table = local_sensitivity(scenario, parameters, outputs, relative_step=args.rel_step)
    out_path = PROJECT_ROOT / "results" / "tables" / f"{args.prefix}.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path, index=False)
    print(f"Saved sensitivity table: {out_path}")
    print(table.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
