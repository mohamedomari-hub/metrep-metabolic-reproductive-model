"""Run all-parameter local sensitivity analysis."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from metrep.analysis import local_sensitivity
from metrep.parameters import PARAMETERS
from metrep.scenarios import constant_non_lactating


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--days", type=float, default=90.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--rel-step", type=float, default=0.01)
    parser.add_argument("--prefix", default="local_sensitivity_allparams")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scenario = constant_non_lactating(days=args.days, step=args.dt)
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
