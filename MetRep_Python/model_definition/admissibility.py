"""Biological admissibility checks for MetRep trajectories."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from model_definition.initial_conditions import STATE_INDEX
from model_definition.simulate import SimulationResult


@dataclass(frozen=True)
class AdmissibilityThresholds:
    """Thresholds for trajectory similarity to the nominal reference."""

    min_correlation: float = 0.75
    max_average_difference: float = 0.30
    max_norm_difference: float = 0.30
    max_shift_days: float = 5.0
    penalty_weight: float = 100.0


def _safe_denominator(value: float, epsilon: float = 1e-10) -> float:
    return float(value) if abs(float(value)) > epsilon else epsilon


def _max_shifted_correlation(
    reference: np.ndarray,
    candidate: np.ndarray,
    dt: float,
    max_shift_days: float,
) -> float:
    """Return max normalized inner-product similarity over integer time shifts."""

    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    max_lag = int(round(max_shift_days / max(dt, 1e-12)))
    max_lag = min(max_lag, max(0, reference.size - 2))
    best = -np.inf

    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            ref_slice = reference[-lag:]
            cand_slice = candidate[: candidate.size + lag]
        elif lag > 0:
            ref_slice = reference[: reference.size - lag]
            cand_slice = candidate[lag:]
        else:
            ref_slice = reference
            cand_slice = candidate

        if ref_slice.size < 2 or cand_slice.size < 2:
            continue
        numerator = float(np.trapz(ref_slice * cand_slice, dx=dt))
        denominator = np.sqrt(float(np.trapz(ref_slice * ref_slice, dx=dt))) * np.sqrt(
            float(np.trapz(cand_slice * cand_slice, dx=dt))
        )
        best = max(best, numerator / _safe_denominator(denominator))

    return float(best if np.isfinite(best) else 0.0)


def trajectory_admissibility(
    reference: SimulationResult,
    candidate: SimulationResult,
    states: list[str],
    thresholds: AdmissibilityThresholds | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Compare a candidate trajectory with a reference trajectory.

    The metrics follow the thesis-style admissibility checks:
    shifted normalized cross-correlation, normalized average absolute
    difference, and normalized squared-norm difference.
    """

    thresholds = AdmissibilityThresholds() if thresholds is None else thresholds
    t = np.asarray(reference.t, dtype=float)
    dt = float(np.median(np.diff(t))) if t.size > 1 else 1.0
    rows = []

    for state in states:
        index = STATE_INDEX[state]
        ref = np.asarray(reference.y[:, index], dtype=float)
        cand = np.interp(t, candidate.t, candidate.y[:, index])

        correlation = _max_shifted_correlation(ref, cand, dt, thresholds.max_shift_days)
        average_difference = float(
            np.trapz(np.abs(ref - cand), t) / _safe_denominator(np.trapz(np.abs(ref), t))
        )
        ref_norm = np.sqrt(float(np.trapz(ref * ref, t)))
        cand_norm = np.sqrt(float(np.trapz(cand * cand, t)))
        norm_difference = abs(ref_norm - cand_norm) / _safe_denominator(ref_norm)

        correlation_violation = max(0.0, thresholds.min_correlation - correlation)
        average_violation = max(0.0, average_difference - thresholds.max_average_difference)
        norm_violation = max(0.0, norm_difference - thresholds.max_norm_difference)
        penalty = thresholds.penalty_weight * (
            correlation_violation**2 + average_violation**2 + norm_violation**2
        )

        rows.append(
            {
                "state": state,
                "max_correlation": correlation,
                "average_difference": average_difference,
                "norm_difference": norm_difference,
                "admissible": (
                    correlation >= thresholds.min_correlation
                    and average_difference <= thresholds.max_average_difference
                    and norm_difference <= thresholds.max_norm_difference
                ),
                "penalty": penalty,
            }
        )

    metrics = pd.DataFrame(rows)
    if metrics.empty:
        summary = {
            "admissible": True,
            "penalty": 0.0,
            "min_correlation": np.nan,
            "max_average_difference": np.nan,
            "max_norm_difference": np.nan,
            "worst_species": "",
        }
        return metrics, summary

    worst_index = metrics["penalty"].to_numpy(dtype=float).argmax()
    summary = {
        "admissible": bool(metrics["admissible"].all()),
        "penalty": float(metrics["penalty"].sum()),
        "min_correlation": float(metrics["max_correlation"].min()),
        "max_average_difference": float(metrics["average_difference"].max()),
        "max_norm_difference": float(metrics["norm_difference"].max()),
        "worst_species": str(metrics.iloc[worst_index]["state"]),
    }
    return metrics, summary
