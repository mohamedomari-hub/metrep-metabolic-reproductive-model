"""State definitions and MATLAB-derived initial conditions."""

from __future__ import annotations

import numpy as np


STATE_NAMES: list[str] = [
    "GnRH_H",
    "GnRH",
    "FSH_pituitary",
    "FSH",
    "LH_pituitary",
    "LH",
    "Follicle",
    "PGF",
    "CL",
    "P4",
    "E2",
    "INH",
    "ENZ",
    "OXT",
    "IOF",
    "IGF1",
    "Insulin",
    "Glucose",
    "Fat",
    "Liver_Glucose",
    "Glucose_Store",
    "Glucagon",
]

STATE_INDEX: dict[str, int] = {name: index for index, name in enumerate(STATE_NAMES)}

DEXA_STATE_NAMES: list[str] = ["A_dep", "A_cent", "C_e"]
EXTENDED_STATE_NAMES: list[str] = STATE_NAMES + DEXA_STATE_NAMES
EXTENDED_STATE_INDEX: dict[str, int] = {
    name: index for index, name in enumerate(EXTENDED_STATE_NAMES)
}

MATLAB_INITIAL_CONDITIONS_25 = np.array(
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

CORE_INITIAL_CONDITIONS = MATLAB_INITIAL_CONDITIONS_25[:22].copy()
CORE_INITIAL_CONDITIONS[STATE_INDEX["Insulin"]] = 15.20754495
CORE_INITIAL_CONDITIONS[STATE_INDEX["Glucose"]] = 0.47614015


def matlab_core_initial_conditions() -> np.ndarray:
    """Return the unadjusted first 22 MATLAB initial values."""

    return MATLAB_INITIAL_CONDITIONS_25[:22].copy()


def initial_conditions(mode: str = "non_lactating") -> np.ndarray:
    """Return initial conditions for the 22-state core model.

    Parameters
    ----------
    mode:
        ``"non_lactating"`` uses the MATLAB baseline values. ``"lactating"``
        applies the MATLAB run-file adjustment ``y0(14)=2.5`` for OXT.
    """

    y0 = CORE_INITIAL_CONDITIONS.copy()

    if mode == "lactating":
        y0[STATE_INDEX["OXT"]] = 2.5
    elif mode != "non_lactating":
        raise ValueError("mode must be 'non_lactating' or 'lactating'")

    return y0


def dexa_initial_conditions(mode: str = "non_lactating") -> np.ndarray:
    """Return 25-state initial conditions for the optional Dexa extension."""

    y0 = np.zeros(len(EXTENDED_STATE_NAMES), dtype=float)
    y0[: len(STATE_NAMES)] = initial_conditions(mode=mode)
    return y0
