"""Validate Python constants against the MATLAB reference metadata.

Run from project root:
    python scripts/01_validate_against_matlab.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from metrep.initial_conditions import (
    CORE_INITIAL_CONDITIONS,
    MATLAB_INITIAL_CONDITIONS_25,
    STATE_NAMES,
    matlab_core_initial_conditions,
)
from metrep.parameters import matlab_parameter_vector, parameter_table


RESULTS_DIR = PROJECT_ROOT / "results" / "tables"


MATLAB_Y0_25 = np.array(
    [
        0.667,
        0.551,
        0.316,
        0.395,
        1.000,
        0.642,
        1.000,
        0.00506,
        0.0,
        0.004,
        0.89,
        0.826,
        0.0,
        0.0183,
        0.35,
        0.48,
        15.5,
        0.48,
        15e4,
        110.0,
        535.0,
        105.0,
        0.0,
        0.0,
        0.0,
    ],
    dtype=float,
)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    parameters = pd.DataFrame(parameter_table())
    parameters.to_csv(RESULTS_DIR / "validated_parameters_from_matlab.csv", index=False)

    state_table = pd.DataFrame(
        {
            "python_index_0_based": range(len(STATE_NAMES)),
            "matlab_index_1_based": range(1, len(STATE_NAMES) + 1),
            "state_name": STATE_NAMES,
            "initial_value": CORE_INITIAL_CONDITIONS,
            "matlab_reference_initial_value": matlab_core_initial_conditions(),
        }
    )
    state_table.to_csv(RESULTS_DIR / "validated_state_order_from_matlab.csv", index=False)

    checks = {
        "parameter_count_is_98": len(matlab_parameter_vector()) == 98,
        "matlab_y0_25_matches_documented_reference": np.allclose(
            MATLAB_INITIAL_CONDITIONS_25, MATLAB_Y0_25
        ),
        "core_state_count_is_22": len(STATE_NAMES) == 22,
        "raw_core_y0_matches_first_22_matlab_states": np.allclose(
            matlab_core_initial_conditions(), MATLAB_Y0_25[:22]
        ),
    }
    report = pd.DataFrame(
        [{"check": name, "passed": passed} for name, passed in checks.items()]
    )
    report.to_csv(RESULTS_DIR / "validation_report.csv", index=False)

    print("Validation summary")
    print(report.to_string(index=False))
    print()
    print("Simulation-output validation TODO:")
    print("- Export MATLAB no-Dexa baseline output as CSV or MAT.")
    print("- Expected columns: time plus states 1..22 in BovSys_Equa_dexa_v3 order.")
    print("- Place it under data/matlab_reference/ and compare trajectories here.")


if __name__ == "__main__":
    main()
