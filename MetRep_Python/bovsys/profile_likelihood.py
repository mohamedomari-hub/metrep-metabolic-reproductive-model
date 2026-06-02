"""Profile-likelihood utilities for synthetic BovSys experiments."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from bovsys.admissibility import (
    AdmissibilityThresholds,
    trajectory_admissibility,
)
from bovsys.initial_conditions import STATE_INDEX
from bovsys.parameters import default_parameters
from bovsys.scenarios import Scenario
from bovsys.simulate import SimulationResult, run_simulation


@dataclass(frozen=True)
class SyntheticData:
    """Synthetic observations sampled from a nominal simulation."""

    times: np.ndarray
    outputs: list[str]
    values: np.ndarray
    sigma: np.ndarray
    noise_relative: float


def make_synthetic_data(
    reference: SimulationResult,
    outputs: list[str],
    sample_every_days: float = 2.0,
    noise_relative: float = 0.05,
    seed: int = 42,
) -> SyntheticData:
    """Sample measurable outputs from a reference trajectory."""

    rng = np.random.default_rng(seed)
    times = np.arange(float(reference.t[0]), float(reference.t[-1]) + 1e-9, sample_every_days)
    values = np.zeros((times.size, len(outputs)), dtype=float)
    for col, output in enumerate(outputs):
        values[:, col] = np.interp(times, reference.t, reference.y[:, STATE_INDEX[output]])

    scale = np.maximum(np.abs(values), np.maximum(1e-6, 0.05 * np.nanmax(np.abs(values), axis=0)))
    sigma = np.maximum(noise_relative, 1e-6) * scale
    if noise_relative > 0:
        values = np.maximum(values + rng.normal(0.0, sigma), 0.0)

    return SyntheticData(
        times=times,
        outputs=list(outputs),
        values=values,
        sigma=sigma,
        noise_relative=float(noise_relative),
    )


def observation_loss(result: SimulationResult, data: SyntheticData) -> float:
    """Weighted least-squares loss against synthetic observations."""

    prediction = np.zeros_like(data.values)
    for col, output in enumerate(data.outputs):
        prediction[:, col] = np.interp(data.times, result.t, result.y[:, STATE_INDEX[output]])
    residual = (prediction - data.values) / data.sigma
    return 0.5 * float(np.sum(residual * residual))


def _scaled_params(
    base_params: dict[str, float],
    names: list[str],
    log_multipliers: np.ndarray,
) -> dict[str, float]:
    params = dict(base_params)
    for name, log_multiplier in zip(names, log_multipliers):
        params[name] = float(base_params[name]) * float(np.exp(log_multiplier))
    return params


def profile_one_parameter(
    scenario: Scenario,
    reference: SimulationResult,
    data: SyntheticData,
    profile_parameter: str,
    nuisance_parameters: list[str],
    grid_multipliers: np.ndarray,
    admissibility_states: list[str],
    thresholds: AdmissibilityThresholds,
    nuisance_bounds: tuple[float, float] = (0.5, 1.5),
    maxiter: int = 60,
) -> pd.DataFrame:
    """Profile one parameter while re-optimizing selected nuisance parameters."""

    base_params = default_parameters()
    nuisance_parameters = [name for name in nuisance_parameters if name != profile_parameter]
    bounds = [(np.log(nuisance_bounds[0]), np.log(nuisance_bounds[1])) for _ in nuisance_parameters]
    rows = []

    for multiplier in grid_multipliers:
        fixed_value = float(base_params[profile_parameter]) * float(multiplier)

        def objective(log_multipliers: np.ndarray) -> float:
            params = _scaled_params(base_params, nuisance_parameters, log_multipliers)
            params[profile_parameter] = fixed_value
            try:
                result = run_simulation(scenario, params)
                fit_loss = observation_loss(result, data)
                _, admissibility = trajectory_admissibility(
                    reference,
                    result,
                    admissibility_states,
                    thresholds,
                )
                return fit_loss + float(admissibility["penalty"])
            except Exception:
                return 1e12

        if nuisance_parameters:
            opt = minimize(
                objective,
                np.zeros(len(nuisance_parameters), dtype=float),
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": int(maxiter), "ftol": 1e-6},
            )
            best_log = opt.x
            optimizer_success = bool(opt.success)
            optimizer_message = str(opt.message)
        else:
            best_log = np.zeros(0, dtype=float)
            optimizer_success = True
            optimizer_message = "no nuisance parameters"

        best_params = _scaled_params(base_params, nuisance_parameters, best_log)
        best_params[profile_parameter] = fixed_value
        try:
            best_result = run_simulation(scenario, best_params)
            fit_loss = observation_loss(best_result, data)
            _, admissibility = trajectory_admissibility(
                reference,
                best_result,
                admissibility_states,
                thresholds,
            )
            admissibility_penalty = float(admissibility["penalty"])
            total_loss = fit_loss + admissibility_penalty
        except Exception as exc:
            fit_loss = 1e12
            admissibility_penalty = 1e12
            total_loss = 2e12
            admissibility = {
                "admissible": False,
                "min_correlation": np.nan,
                "max_average_difference": np.nan,
                "max_norm_difference": np.nan,
                "worst_species": "solver_failure",
            }
            optimizer_success = False
            optimizer_message = str(exc)

        rows.append(
            {
                "profile_parameter": profile_parameter,
                "fixed_multiplier": float(multiplier),
                "fixed_value": fixed_value,
                "fit_loss": fit_loss,
                "admissibility_penalty": admissibility_penalty,
                "total_loss": total_loss,
                "delta_loss": np.nan,
                "admissible": bool(admissibility["admissible"]),
                "min_correlation": admissibility["min_correlation"],
                "max_average_difference": admissibility["max_average_difference"],
                "max_norm_difference": admissibility["max_norm_difference"],
                "worst_species": admissibility["worst_species"],
                "optimizer_success": optimizer_success,
                "optimizer_message": optimizer_message,
            }
        )

    table = pd.DataFrame(rows)
    table["delta_loss"] = table["total_loss"] - float(table["total_loss"].min())
    return table


def classify_profile(
    profile: pd.DataFrame,
    threshold: float = 1.92,
    admissible_only: bool = False,
) -> dict[str, object]:
    """Classify a profile curve with lightweight diagnostics."""

    profile_for_class = profile[profile["admissible"]].copy() if admissible_only else profile.copy()
    if profile_for_class.empty:
        parameter = str(profile["profile_parameter"].iloc[0])
        return {
            "parameter": parameter,
            "classification": "no admissible profile points",
            "best_multiplier": np.nan,
            "best_value": np.nan,
            "min_total_loss": np.nan,
            "profile_span": np.nan,
            "admissible_fraction": 0.0,
            "worst_species": str(profile.loc[profile["admissibility_penalty"].idxmax()]["worst_species"]),
        }

    best = profile_for_class.loc[profile_for_class["total_loss"].idxmin()]
    admissible_fraction = float(profile["admissible"].mean())
    edge_best = bool(best.name == profile_for_class.index.min() or best.name == profile_for_class.index.max())
    delta = profile_for_class["total_loss"] - float(profile_for_class["total_loss"].min())
    span = float(delta.max() - delta.min())
    below = profile_for_class[delta <= threshold]

    if span < threshold:
        classification = "flat/non-identifiable"
    elif edge_best:
        classification = "boundary-limited"
    elif len(below) >= 0.8 * len(profile):
        classification = "weakly identifiable"
    else:
        classification = "practically identifiable"

    return {
        "parameter": str(best["profile_parameter"]),
        "classification": classification,
        "best_multiplier": float(best["fixed_multiplier"]),
        "best_value": float(best["fixed_value"]),
        "min_total_loss": float(best["total_loss"]),
        "profile_span": span,
        "admissible_fraction": admissible_fraction,
        "worst_species": str(profile.loc[profile["admissibility_penalty"].idxmax()]["worst_species"]),
    }
