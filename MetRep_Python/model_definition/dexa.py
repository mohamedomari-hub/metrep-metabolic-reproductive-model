"""Optional dexamethasone PK/PD extension for MetRep.

The implementation follows the MATLAB reference files
``BovSys_Equa_dexa_v3.m`` and ``BovSys_run_dexa_v3.m``. Ordinary MetRep
scenarios should continue to use ``run_simulation``; this module is only used
when a Dexa dose is explicitly requested.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from model_definition.initial_conditions import (
    EXTENDED_STATE_NAMES,
    dexa_initial_conditions,
)
from model_definition.ode_model import metrep_rhs
from model_definition.parameters import default_parameters
from model_definition.scenarios import Scenario
from model_definition.simulate import SimulationResult


@dataclass(frozen=True)
class DexaConfig:
    """Dexamethasone dose and PK/PD constants from the MATLAB v3 runner."""

    enabled: bool = True
    dose_mg_per_kg: float = 0.02
    body_weight_kg: float = 600.0
    dose_day: float = 0.0
    ka: float = 13.4352
    ke: float = 2.7086
    bioavailability: float = 0.72
    keo: float = 0.7
    vd_l_per_kg: float = 1.105
    emax_glucagon: float = 3.0
    ca_ng_per_ml: float = 1.8
    cb_ng_per_ml: float = 1.8

    @property
    def dose_ng(self) -> float:
        """Injected IM dose in ng."""

        return self.dose_mg_per_kg * self.body_weight_kg * 1e6

    @property
    def vd_ml(self) -> float:
        """Central volume of distribution in mL."""

        return self.vd_l_per_kg * self.body_weight_kg * 1e3


def dexa_effects(ce_ng_per_ml: float, config: DexaConfig) -> tuple[float, float]:
    """Return glucagon-stimulation and metabolic-inhibition multipliers."""

    ce = max(float(ce_ng_per_ml), 0.0)
    ca10 = config.ca_ng_per_ml**10
    ce10 = ce**10
    cb7 = config.cb_ng_per_ml**7
    ce7 = ce**7

    effect_gluca = 1.0
    if ca10 + ce10 > 0.0:
        effect_gluca += config.emax_glucagon * ce10 / (ca10 + ce10)

    effect_bt = 1.0
    if cb7 + ce7 > 0.0:
        effect_bt -= ce7 / (cb7 + ce7)

    return effect_gluca, max(effect_bt, 0.0)


def dexa_rhs(
    t: float,
    y: np.ndarray,
    params: dict[str, float],
    at: np.ndarray,
    dmi: np.ndarray,
    milk: np.ndarray,
    mode: str,
    config: DexaConfig,
) -> np.ndarray:
    """Return derivatives for the 25-state Dexa-extended model."""

    y = np.maximum(np.asarray(y, dtype=float), 0.0)
    f = np.zeros(len(EXTENDED_STATE_NAMES), dtype=float)

    if config.enabled:
        effect_gluca, effect_bt = dexa_effects(y[24], config)
    else:
        effect_gluca, effect_bt = 1.0, 1.0

    f[:22] = metrep_rhs(
        t,
        y[:22],
        params,
        at,
        dmi,
        milk,
        mode,
        effect_gluca=effect_gluca,
        effect_bt=effect_bt,
    )

    if not config.enabled:
        return f

    adep = y[22]
    acent = y[23]
    ce = y[24]
    central_concentration = acent / config.vd_ml

    f[22] = -config.ka * adep
    f[23] = config.bioavailability * config.ka * adep - config.ke * acent
    f[24] = config.keo * (central_concentration - ce)
    return f


def _solve_segment(
    y0: np.ndarray,
    t_eval: np.ndarray,
    params: dict[str, float],
    scenario: Scenario,
    config: DexaConfig,
    method: str,
    rtol: float,
    atol: float,
) -> tuple[np.ndarray, np.ndarray]:
    if t_eval.size == 1:
        return t_eval, y0.reshape(1, -1)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return dexa_rhs(
            t,
            y,
            params,
            scenario.t_eval,
            scenario.dmi,
            scenario.milk,
            scenario.mode,
            config,
        )

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
        raise RuntimeError(f"Dexa ODE solver failed for {scenario.name}: {sol.message}")
    if not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f"Dexa ODE solver produced non-finite values for {scenario.name}")
    return sol.t, np.maximum(sol.y.T, 0.0)


def run_dexa_simulation(
    scenario: Scenario,
    config: DexaConfig | None = None,
    parameters: dict[str, float] | None = None,
    method: str = "BDF",
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> SimulationResult:
    """Solve a scenario with the optional Dexa dose and PK/PD states."""

    config = DexaConfig() if config is None else config
    params = default_parameters() if parameters is None else dict(parameters)
    params["c0"] = scenario.c0
    t_eval = np.asarray(scenario.t_eval, dtype=float)
    y0 = dexa_initial_conditions(mode=scenario.mode)

    if not config.enabled:
        t, y = _solve_segment(y0, t_eval, params, scenario, config, method, rtol, atol)
    elif config.dose_day <= float(t_eval[0]):
        y0[22] += config.dose_ng
        t, y = _solve_segment(y0, t_eval, params, scenario, config, method, rtol, atol)
    elif config.dose_day >= float(t_eval[-1]):
        raise ValueError("dose_day must fall inside the simulation horizon.")
    else:
        dose_day = float(config.dose_day)
        t_before = t_eval[t_eval < dose_day]
        if t_before.size == 0 or not np.isclose(t_before[-1], dose_day):
            t_before = np.append(t_before, dose_day)
        t_after = t_eval[t_eval > dose_day]
        t_after = np.insert(t_after, 0, dose_day)

        t1, y1 = _solve_segment(y0, t_before, params, scenario, config, method, rtol, atol)
        ydose = y1[-1].copy()
        ydose[22] += config.dose_ng
        t2, y2 = _solve_segment(ydose, t_after, params, scenario, config, method, rtol, atol)
        t = np.concatenate([t1, t2[1:]])
        y = np.vstack([y1, y2[1:]])

    dmi = np.interp(t, scenario.t_eval, scenario.dmi)
    milk = np.interp(t, scenario.t_eval, scenario.milk)
    name = f"{scenario.name}_dexa" if config.enabled else f"{scenario.name}_no_dexa_25state"
    return SimulationResult(
        scenario_name=name,
        t=t,
        y=y,
        state_names=list(EXTENDED_STATE_NAMES),
        dmi=dmi,
        milk=milk,
        parameters=params,
    )
