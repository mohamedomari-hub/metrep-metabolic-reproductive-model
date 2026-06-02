"""Simulation utilities for BovSys."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from bovsys.initial_conditions import STATE_NAMES, initial_conditions
from bovsys.ode_model import bovsys_rhs
from bovsys.parameters import default_parameters
from bovsys.scenarios import Scenario


@dataclass
class SimulationResult:
    """Container for one solved BovSys trajectory."""

    scenario_name: str
    t: np.ndarray
    y: np.ndarray
    state_names: list[str]
    dmi: np.ndarray
    milk: np.ndarray
    parameters: dict[str, float]

    def to_dataframe(self) -> pd.DataFrame:
        """Return a tidy time-indexed DataFrame with states and forcings."""

        data = {"time_days": self.t, "DMI": self.dmi, "Milk": self.milk}
        for index, name in enumerate(self.state_names):
            data[name] = self.y[:, index]
        return pd.DataFrame(data)


def run_simulation(
    scenario: Scenario,
    parameters: dict[str, float] | None = None,
    method: str = "BDF",
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> SimulationResult:
    """Solve one BovSys scenario with SciPy ``solve_ivp``."""

    params = default_parameters() if parameters is None else dict(parameters)
    params["c0"] = scenario.c0
    y0 = initial_conditions(mode=scenario.mode)
    t_eval = np.asarray(scenario.t_eval, dtype=float)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return bovsys_rhs(t, y, params, t_eval, scenario.dmi, scenario.milk, scenario.mode)

    sol = solve_ivp(
        rhs,
        (float(t_eval[0]), float(t_eval[-1])),
        y0,
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solver failed for {scenario.name}: {sol.message}")
    if not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f"ODE solver produced non-finite values for {scenario.name}")

    y = np.maximum(sol.y.T, 0.0)
    return SimulationResult(
        scenario_name=scenario.name,
        t=sol.t,
        y=y,
        state_names=list(STATE_NAMES),
        dmi=np.asarray(scenario.dmi, dtype=float),
        milk=np.asarray(scenario.milk, dtype=float),
        parameters=params,
    )


def save_result(result: SimulationResult, output_dir: Path) -> tuple[Path, Path]:
    """Save a simulation as CSV and compressed NPZ."""

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{result.scenario_name}.csv"
    npz_path = output_dir / f"{result.scenario_name}.npz"
    result.to_dataframe().to_csv(csv_path, index=False)
    np.savez_compressed(
        npz_path,
        t=result.t,
        y=result.y,
        state_names=np.array(result.state_names, dtype=object),
        dmi=result.dmi,
        milk=result.milk,
        scenario_name=result.scenario_name,
    )
    return csv_path, npz_path
