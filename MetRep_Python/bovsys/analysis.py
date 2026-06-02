"""Sensitivity, identifiability, and uncertainty helpers."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from bovsys.initial_conditions import STATE_INDEX
from bovsys.parameters import default_parameters
from bovsys.scenarios import Scenario
from bovsys.simulate import run_simulation


def state_auc(result, state: str) -> float:
    """Compute trapezoidal AUC for a state trajectory."""

    idx = STATE_INDEX[state]
    return float(np.trapz(result.y[:, idx], result.t))


def local_sensitivity(
    scenario: Scenario,
    parameter_names: list[str],
    outputs: list[str],
    relative_step: float = 0.01,
) -> pd.DataFrame:
    """One-at-a-time local sensitivity of output AUCs to parameters."""

    base_params = default_parameters()
    base = run_simulation(scenario, base_params)
    base_metrics = {output: state_auc(base, output) for output in outputs}
    rows = []

    for name in parameter_names:
        if name not in base_params:
            raise KeyError(f"Unknown parameter: {name}")
        nominal = base_params[name]
        if nominal == 0:
            continue

        perturbed_params = dict(base_params)
        perturbed_params[name] = nominal * (1.0 + relative_step)
        perturbed = run_simulation(scenario, perturbed_params)

        for output in outputs:
            metric = state_auc(perturbed, output)
            sensitivity = ((metric - base_metrics[output]) / base_metrics[output]) / relative_step
            rows.append(
                {
                    "parameter": name,
                    "output": output,
                    "metric": "AUC",
                    "nominal_value": nominal,
                    "base_metric": base_metrics[output],
                    "perturbed_metric": metric,
                    "relative_sensitivity": sensitivity,
                }
            )

    return pd.DataFrame(rows).sort_values("relative_sensitivity", key=np.abs, ascending=False)


def sensitivity_matrix(
    scenario: Scenario,
    parameter_names: list[str],
    outputs: list[str],
    sample_times: np.ndarray,
    relative_step: float = 0.01,
) -> tuple[np.ndarray, list[str]]:
    """Build a normalized local sensitivity matrix for identifiability analysis."""

    base_params = default_parameters()
    base = run_simulation(scenario, base_params)
    rows = []
    row_labels = []

    for output in outputs:
        idx = STATE_INDEX[output]
        base_values = np.interp(sample_times, base.t, base.y[:, idx])
        scale = np.maximum(np.abs(base_values), 1e-8)
        rows.append((base_values, scale))
        row_labels.extend([f"{output}@{time:.1f}d" for time in sample_times])

    matrix = np.zeros((len(outputs) * len(sample_times), len(parameter_names)), dtype=float)

    for col, name in enumerate(parameter_names):
        nominal = base_params[name]
        perturbed_params = dict(base_params)
        perturbed_params[name] = nominal * (1.0 + relative_step)
        perturbed = run_simulation(scenario, perturbed_params)

        offset = 0
        for output_index, output in enumerate(outputs):
            idx = STATE_INDEX[output]
            perturbed_values = np.interp(sample_times, perturbed.t, perturbed.y[:, idx])
            base_values, scale = rows[output_index]
            matrix[offset : offset + len(sample_times), col] = (
                (perturbed_values - base_values) / scale / relative_step
            )
            offset += len(sample_times)

    return matrix, row_labels


def build_output_matrix(states: np.ndarray, outputs: list[str]) -> np.ndarray:
    """Select observable outputs from a full state trajectory."""

    indexes = [STATE_INDEX[output] for output in outputs]
    return states[:, indexes]


def sensitivity_matrix_central_diff(
    scenario: Scenario,
    parameter_names: list[str],
    outputs: list[str],
    relative_step: float = 1e-3,
    absolute_step_min: float = 1e-8,
    fail_penalty: float = 1e3,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the stacked central-difference sensitivity matrix.

    This follows the approach in the older
    ``bovsys_structid_svd_merged_robust_updated.py`` script:

    ``S = d vec(Y_outputs) / d theta``

    where rows are all selected outputs stacked over all sampled time points and
    columns are parameters. Central differences are used for each parameter.
    """

    base_params = default_parameters()
    for name in parameter_names:
        if name not in base_params:
            raise KeyError(f"Unknown parameter: {name}")

    nominal = run_simulation(scenario, base_params)
    y0_outputs = build_output_matrix(nominal.y, outputs)
    y0_vector = y0_outputs.reshape(-1)

    n_rows = y0_vector.size
    n_params = len(parameter_names)
    sensitivity = np.zeros((n_rows, n_params), dtype=float)

    for col, name in enumerate(parameter_names):
        value = float(base_params[name])
        step = max(absolute_step_min, relative_step * (abs(value) if value != 0 else 1.0))

        high_params = dict(base_params)
        low_params = dict(base_params)
        high_params[name] = value + step
        low_params[name] = value - step

        try:
            high = run_simulation(scenario, high_params)
            high_vector = build_output_matrix(high.y, outputs).reshape(-1)
        except Exception:
            high_vector = y0_vector + fail_penalty

        try:
            low = run_simulation(scenario, low_params)
            low_vector = build_output_matrix(low.y, outputs).reshape(-1)
        except Exception:
            low_vector = y0_vector - fail_penalty

        sensitivity[:, col] = (high_vector - low_vector) / (2.0 * step)

    return sensitivity, y0_outputs


def estimate_rank(
    singular_values: np.ndarray,
    tolerance_mode: str = "relative",
    tolerance_value: float = 1e-8,
) -> tuple[int, float]:
    """Estimate numerical rank from singular values."""

    if singular_values.size == 0:
        return 0, 0.0
    threshold = (
        tolerance_value * singular_values[0]
        if tolerance_mode == "relative"
        else tolerance_value
    )
    rank = int(np.sum(singular_values > threshold))
    return rank, float(threshold)


def ranking_scores(
    sensitivity: np.ndarray,
    nominal_outputs: np.ndarray,
    parameters: dict[str, float],
    parameter_names: list[str],
    metric: str = "rel2_colnorm",
    output_epsilon: float = 1e-8,
) -> np.ndarray:
    """Compute local sensitivity ranking metrics from the old SVD workflow."""

    y_vector = nominal_outputs.reshape(-1)
    y_safe = np.maximum(np.abs(y_vector), output_epsilon)
    scores = np.zeros(len(parameter_names), dtype=float)

    for col, name in enumerate(parameter_names):
        parameter_value = float(parameters[name])
        sensitivity_col = sensitivity[:, col]

        if metric == "colnorm":
            scores[col] = np.linalg.norm(sensitivity_col, ord=2)
        elif metric == "rel_colnorm":
            scores[col] = np.linalg.norm(parameter_value * sensitivity_col, ord=2)
        elif metric == "rel2_colnorm":
            scores[col] = np.linalg.norm((parameter_value / y_safe) * sensitivity_col, ord=2)
        else:
            raise ValueError("metric must be 'colnorm', 'rel_colnorm', or 'rel2_colnorm'")

    return scores


def nullspace_participation(nullspace: np.ndarray, n_parameters: int) -> np.ndarray:
    """Aggregate how strongly each parameter appears in nullspace directions."""

    if nullspace.size == 0:
        return np.zeros(n_parameters, dtype=float)
    return np.linalg.norm(nullspace, axis=0)


def identifiable_participation(vt: np.ndarray, rank: int) -> np.ndarray:
    """Aggregate parameter participation in identifiable right-singular vectors."""

    if vt.size == 0 or rank <= 0:
        return np.zeros(vt.shape[1] if vt.ndim == 2 else 0, dtype=float)
    identifiable = vt[:rank, :]
    return np.sqrt(np.sum(identifiable * identifiable, axis=0))


def extract_compensation_edges(
    nullspace: np.ndarray,
    parameter_names: list[str],
    threshold_relative: float = 0.35,
    top_k: int = 10,
) -> pd.DataFrame:
    """Extract aggregate positive-vs-negative compensation pairs."""

    edges: dict[tuple[str, str], float] = {}
    if nullspace.size == 0:
        return pd.DataFrame(columns=["increase_param", "decrease_param", "weight"])

    for row in nullspace:
        values = row.astype(float)
        abs_values = np.abs(values)
        max_abs = abs_values.max()
        if max_abs <= 0:
            continue

        keep = np.where(abs_values >= threshold_relative * max_abs)[0]
        if keep.size > top_k:
            keep = np.argsort(abs_values)[::-1][:top_k]

        positive = [index for index in keep if values[index] > 0]
        negative = [index for index in keep if values[index] < 0]

        for pos in positive:
            for neg in negative:
                key = (parameter_names[pos], parameter_names[neg])
                edges[key] = edges.get(key, 0.0) + float(abs(values[pos] * values[neg]))

    rows = [
        {"increase_param": key[0], "decrease_param": key[1], "weight": weight}
        for key, weight in sorted(edges.items(), key=lambda item: item[1], reverse=True)
    ]
    return pd.DataFrame(rows)


def _normalize_0_1(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    low = float(np.min(values))
    high = float(np.max(values))
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return np.zeros_like(values)
    return (values - low) / (high - low)


def recommend_fix_estimate(
    parameter_names: list[str],
    sensitivity_score: np.ndarray,
    nullspace_score: np.ndarray,
    identifiable_score: np.ndarray,
    sensitivity_high_quantile: float = 0.75,
    sensitivity_low_quantile: float = 0.25,
    null_high_quantile: float = 0.75,
    null_low_quantile: float = 0.25,
) -> pd.DataFrame:
    """Classify parameters using sensitivity plus nullspace involvement."""

    sensitivity_01 = _normalize_0_1(sensitivity_score)
    nullspace_01 = _normalize_0_1(nullspace_score)
    identifiable_01 = _normalize_0_1(identifiable_score)

    sens_hi = float(np.quantile(sensitivity_01, sensitivity_high_quantile))
    sens_lo = float(np.quantile(sensitivity_01, sensitivity_low_quantile))
    null_hi = float(np.quantile(nullspace_01, null_high_quantile))
    null_lo = float(np.quantile(nullspace_01, null_low_quantile))

    recommendations = []
    for sens, null in zip(sensitivity_01, nullspace_01):
        if sens <= sens_lo:
            recommendations.append("FIX_low_sensitivity")
        elif null >= null_hi:
            recommendations.append(
                "COMPENSATORY_high_impact" if sens >= sens_hi else "FIX_anchor_compensation"
            )
        else:
            recommendations.append("ESTIMATE_candidate" if sens >= sens_hi or null <= null_lo else "ESTIMATE_or_constrain")

    collapsed = [collapse_recommendation(recommendation) for recommendation in recommendations]
    table = pd.DataFrame(
        {
            "param": parameter_names,
            "parameter": parameter_names,
            "sensitivity_raw": sensitivity_score,
            "sensitivity_0to1": sensitivity_01,
            "nullspace_raw": nullspace_score,
            "nullspace_0to1": nullspace_01,
            "identifiable_raw": identifiable_score,
            "identifiable_0to1": identifiable_01,
            "recommendation": recommendations,
            "recommendation_3class": collapsed,
        }
    ).sort_values(["sensitivity_0to1", "nullspace_0to1"], ascending=[False, False])
    table.attrs["thresholds"] = {
        "sens_hi": sens_hi,
        "sens_lo": sens_lo,
        "null_hi": null_hi,
        "null_lo": null_lo,
        "sens_q_hi": sensitivity_high_quantile,
        "sens_q_lo": sensitivity_low_quantile,
        "null_q_hi": null_high_quantile,
        "null_q_lo": null_low_quantile,
    }
    return table


def collapse_recommendation(recommendation: str) -> str:
    """Collapse detailed recommendations into the old robust script's 3 classes."""

    if recommendation in {"FIX_low_sensitivity", "FIX_low_sens", "low_sensitivity", "FIX_ignore"}:
        return "Fix (irrelevant)"
    if recommendation in {
        "FIX_anchor_compensation",
        "COMPENSATORY_high_impact",
        "COMPENSATORY",
        "high_nullspace",
    }:
        return "Fix (anchor)"
    return "Estimate"


def structural_identifiability_svd(
    scenario: Scenario,
    parameter_names: list[str],
    outputs: list[str],
    relative_step: float = 1e-3,
    tolerance_mode: str = "relative",
    tolerance_value: float = 1e-8,
    ranking_metric: str = "rel2_colnorm",
    null_threshold_relative: float = 0.35,
    null_top_k: int = 10,
) -> dict[str, object]:
    """Run the old script's central-difference/SVD identifiability workflow."""

    base_params = default_parameters()
    sensitivity, nominal_outputs = sensitivity_matrix_central_diff(
        scenario,
        parameter_names=parameter_names,
        outputs=outputs,
        relative_step=relative_step,
    )
    _, singular_values, vt = np.linalg.svd(sensitivity, full_matrices=False)
    rank, threshold = estimate_rank(singular_values, tolerance_mode, tolerance_value)
    nullspace = vt[rank:, :]
    nullity = int(nullspace.shape[0])

    scores = ranking_scores(
        sensitivity,
        nominal_outputs,
        base_params,
        parameter_names,
        metric=ranking_metric,
    )
    null_score = nullspace_participation(nullspace, len(parameter_names))
    id_score = identifiable_participation(vt, rank)
    holistic = recommend_fix_estimate(parameter_names, scores, null_score, id_score)
    edges = extract_compensation_edges(
        nullspace,
        parameter_names,
        threshold_relative=null_threshold_relative,
        top_k=null_top_k,
    )

    return {
        "sensitivity_matrix": sensitivity,
        "nominal_outputs": nominal_outputs,
        "singular_values": singular_values,
        "vt": vt,
        "rank": rank,
        "nullity": nullity,
        "threshold": threshold,
        "nullspace": nullspace,
        "ranking_scores": scores,
        "nullspace_score": null_score,
        "identifiable_score": id_score,
        "holistic_table": holistic,
        "holistic_thresholds": dict(holistic.attrs.get("thresholds", {})),
        "compensation_edges": edges,
    }


def identifiability_svd(
    scenario: Scenario,
    parameter_names: list[str],
    outputs: list[str],
    sample_times: np.ndarray,
    relative_step: float = 0.01,
) -> dict[str, object]:
    """Practical identifiability screen using local sensitivity SVD."""

    matrix, row_labels = sensitivity_matrix(
        scenario,
        parameter_names=parameter_names,
        outputs=outputs,
        sample_times=sample_times,
        relative_step=relative_step,
    )
    u, singular_values, vt = np.linalg.svd(matrix, full_matrices=False)
    condition_number = float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0 else np.inf
    participation = pd.DataFrame(
        np.abs(vt[-1, :]),
        index=parameter_names,
        columns=["least_identifiable_direction_abs_loading"],
    ).sort_values("least_identifiable_direction_abs_loading", ascending=False)

    return {
        "sensitivity_matrix": matrix,
        "row_labels": row_labels,
        "singular_values": singular_values,
        "condition_number": condition_number,
        "least_identifiable_participation": participation,
    }


def sample_parameters(
    base_params: dict[str, float],
    parameter_names: list[str],
    n_samples: int,
    relative_bounds: tuple[float, float] = (0.9, 1.1),
    seed: int = 0,
) -> list[dict[str, float]]:
    """Sample parameters independently in log-uniform relative bounds."""

    rng = np.random.default_rng(seed)
    lo, hi = np.log(relative_bounds[0]), np.log(relative_bounds[1])
    samples = []

    for _ in range(n_samples):
        params = dict(base_params)
        for name in parameter_names:
            params[name] = params[name] * float(np.exp(rng.uniform(lo, hi)))
        samples.append(params)

    return samples


def uncertainty_trajectories(
    scenario: Scenario,
    parameter_names: list[str],
    state: str,
    n_samples: int = 50,
    relative_bounds: tuple[float, float] = (0.95, 1.05),
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Run an ensemble and return trajectories for one state."""

    base_params = default_parameters()
    sampled = sample_parameters(base_params, parameter_names, n_samples, relative_bounds, seed)
    idx = STATE_INDEX[state]
    trajectories = []

    for params in sampled:
        result = run_simulation(replace(scenario, name=f"{scenario.name}_sample"), params)
        trajectories.append(result.y[:, idx])

    return scenario.t_eval, np.asarray(trajectories)
