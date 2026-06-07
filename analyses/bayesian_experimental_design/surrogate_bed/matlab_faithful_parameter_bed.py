"""MATLAB-faithful parameter BED from existing admissible MetRep simulations.

This script never runs the ODE model and never trains a surrogate. It reuses
the stored parameter matrix and stored ODE outputs to reproduce the
non-lactating MATLAB KDE/Gaussian-mixture BED logic:

1. global single-day ranking I(Theta; Y_day), where Theta is the full
   parameter vector and Y_day contains all selected biomarkers;
2. individual biomarker ranking I(Theta; Y_i) at the best global day;
3. parameter-specific I(theta_i; Y_day) rankings;
4. cumulative biomarker panels ordered by parameter-specific information;
5. likelihood-weighted prior/posterior density plots.

The density-ratio calculation is evaluated blockwise. Optional evaluation and
reference subsampling control runtime without changing the estimator form.
Set both limits to 0 to use every admissible sample.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import gaussian_kde


DEFAULT_OBSERVABLE_ORDER = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]
DEFAULT_TARGETS = [
    "insulin_glucose_threshold",
    "inhibin_clearance",
    "hp_p4_follicle_scale",
    "blood_to_liver_glucose_threshold",
    "gnrh_clearance",
    "hp_iof_threshold",
    "insulin_igf_threshold",
    "feed_direct_blood_fraction",
    "lh_basal_release",
]
TARGET_CLASSES = {
    "insulin_glucose_threshold": "Practically identifiable",
    "inhibin_clearance": "Practically identifiable",
    "hp_p4_follicle_scale": "Practically identifiable",
    "blood_to_liver_glucose_threshold": "Boundary-limited",
    "gnrh_clearance": "Boundary-limited",
    "hp_iof_threshold": "Boundary-limited",
    "insulin_igf_threshold": "Weak/flat",
    "feed_direct_blood_fraction": "Weak/flat",
    "lh_basal_release": "Weak/flat",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce MATLAB-style non-lactating parameter BED from existing "
            "biologically admissible ODE simulations. No ODE simulations are run."
        )
    )
    parser.add_argument("--parameters-csv", type=Path, required=True)
    parser.add_argument("--outputs-csv", type=Path, required=True)
    parser.add_argument("--nominal-output-csv", type=Path, required=True)
    parser.add_argument(
        "--measurement-noise-csv",
        type=Path,
        required=True,
        help="MATLAB Abs_Error_NL equivalent computed from the full nominal trajectory.",
    )
    parser.add_argument("--admissibility-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--species",
        nargs="+",
        help="Candidate observable biomarkers. Omit to use all stored observable *_day_* outputs.",
    )
    parser.add_argument("--target-parameters", nargs="+", default=DEFAULT_TARGETS)
    parser.add_argument("--day-start", type=int, default=54)
    parser.add_argument("--day-end", type=int, default=89)
    parser.add_argument("--noise-floor", type=float, default=1e-10)
    parser.add_argument("--mi-repeats", type=int, default=1)
    parser.add_argument(
        "--evaluation-samples",
        type=int,
        default=0,
        help="Monte Carlo evaluation rows per MI estimate; 0 uses all admissible rows.",
    )
    parser.add_argument(
        "--reference-samples",
        type=int,
        default=0,
        help="Gaussian-mixture reference rows per MI estimate; 0 uses all admissible rows.",
    )
    parser.add_argument("--block-size", type=int, default=48)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--posterior-grid-size", type=int, default=500)
    parser.add_argument(
        "--posterior-grid-margin",
        type=float,
        default=0.20,
        help="Fraction of the sampled parameter span added to both posterior KDE plot limits.",
    )
    parser.add_argument(
        "--analysis-label",
        default="Faithful PhD +/-5% bank",
        help="Short bank label shown on figures and recorded in the input audit.",
    )
    parser.add_argument(
        "--posterior-observation-index",
        type=int,
        default=0,
        help="Accepted simulation row used as the fixed synthetic observation z*.",
    )
    parser.add_argument(
        "--max-global-days",
        type=int,
        help="Optional diagnostic limit. Omit for every available requested day.",
    )
    parser.add_argument(
        "--skip-global",
        action="store_true",
        help="Skip expensive full-vector global MI and use --best-day.",
    )
    parser.add_argument("--best-day", type=int, help="Required when --skip-global is used.")
    parser.add_argument(
        "--required-prior-half-range",
        type=float,
        default=0.05,
        help="Required relative prior half-range. The PhD workflow uses 0.05 (+/-5%%).",
    )
    parser.add_argument(
        "--prior-range-tolerance",
        type=float,
        default=0.005,
        help="Absolute tolerance when validating the stored prior half-range.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    table = pd.read_csv(path)
    if table.empty:
        raise ValueError(f"Empty table: {path}")
    return table


def available_days(outputs: pd.DataFrame, species: list[str]) -> list[int]:
    days = set()
    pattern = re.compile(r"^(.+)_day_(\d+)$")
    for column in outputs.columns:
        match = pattern.match(column)
        if match and match.group(1) in species:
            days.add(int(match.group(2)))
    return sorted(days)


def infer_observable_species(outputs: pd.DataFrame) -> list[str]:
    found = set()
    pattern = re.compile(r"^(.+)_day_(\d+)$")
    for column in outputs.columns:
        match = pattern.match(column)
        if match:
            found.add(match.group(1))
    ordered = [name for name in DEFAULT_OBSERVABLE_ORDER if name in found]
    ordered.extend(name for name in sorted(found) if name not in ordered)
    return ordered


def day_columns(day: int, species: list[str]) -> list[str]:
    return [f"{name}_day_{day}" for name in species]


def load_inputs(args: argparse.Namespace):
    parameters = read_csv(args.parameters_csv)
    outputs = read_csv(args.outputs_csv)
    nominal = read_csv(args.nominal_output_csv)
    admissibility = read_csv(args.admissibility_csv)
    if not (len(parameters) == len(outputs) == len(admissibility)):
        raise ValueError("Parameters, outputs, and admissibility tables must have matching rows.")
    if "admissible" not in admissibility:
        raise ValueError("Admissibility table must contain an 'admissible' column.")
    if args.species is None:
        args.species = infer_observable_species(outputs)
    if not args.species:
        raise ValueError("No observable biomarker day columns were found.")

    keep = admissibility["admissible"].astype(bool).to_numpy()
    parameters = parameters.loc[keep].reset_index(drop=True)
    outputs = outputs.loc[keep].reset_index(drop=True)
    if len(parameters) < 100:
        raise ValueError("Fewer than 100 biologically admissible rows remain.")

    unknown_targets = sorted(set(args.target_parameters).difference(parameters.columns))
    if unknown_targets:
        raise ValueError(f"Unknown target parameters: {unknown_targets}")

    stored_days = available_days(outputs, args.species)
    requested = list(range(args.day_start, args.day_end + 1))
    days = [day for day in requested if day in stored_days]
    missing = [day for day in requested if day not in stored_days]
    if not days:
        raise ValueError("No requested days are available in the stored output table.")
    for day in days:
        missing_columns = sorted(set(day_columns(day, args.species)).difference(outputs.columns))
        if missing_columns:
            raise ValueError(f"Day {day} is missing output columns: {missing_columns}")
    if args.max_global_days:
        days = days[: args.max_global_days]
    return parameters, outputs, nominal, days, missing


def robust_bandwidth(values: np.ndarray, floor: float = 1e-12) -> np.ndarray:
    """Approximate MATLAB univariate ksdensity bandwidth per dimension."""

    values = np.asarray(values, dtype=float)
    n = values.shape[0]
    std = np.std(values, axis=0, ddof=1)
    q25, q75 = np.percentile(values, [25, 75], axis=0)
    robust_scale = np.minimum(std, (q75 - q25) / 1.34)
    robust_scale = np.where(robust_scale > floor, robust_scale, std)
    robust_scale = np.where(robust_scale > floor, robust_scale, 1.0)
    return np.maximum(1.06 * robust_scale * n ** (-1.0 / 5.0), floor)


def load_measurement_noise(path: Path, species: list[str], noise_floor: float) -> dict[str, float]:
    """Load the MATLAB Abs_Error_NL equivalent from the full nominal trajectory."""

    table = read_csv(path)
    if not {"biomarker", "abs_error"}.issubset(table.columns):
        raise ValueError("Measurement-noise CSV must contain biomarker and abs_error columns.")
    noise = dict(zip(table["biomarker"], table["abs_error"]))
    missing = sorted(set(species).difference(noise))
    if missing:
        raise ValueError(f"Measurement-noise CSV is missing biomarkers: {missing}")
    return {name: max(float(noise[name]), noise_floor) for name in species}


def choose_indices(n: int, limit: int, rng: np.random.Generator) -> np.ndarray:
    if limit <= 0 or limit >= n:
        return np.arange(n)
    return np.sort(rng.choice(n, size=limit, replace=False))


def prior_range_audit(parameters: pd.DataFrame) -> dict[str, float]:
    """Summarize the relative half-range represented by the stored ensemble."""

    values = parameters.to_numpy(dtype=float)
    centers = np.mean(values, axis=0)
    half_ranges = 0.5 * (np.max(values, axis=0) - np.min(values, axis=0))
    relative = np.divide(
        half_ranges,
        np.abs(centers),
        out=np.full_like(half_ranges, np.nan),
        where=np.abs(centers) > 1e-12,
    )
    finite = relative[np.isfinite(relative)]
    return {
        "relative_half_range_min": float(np.min(finite)),
        "relative_half_range_median": float(np.median(finite)),
        "relative_half_range_max": float(np.max(finite)),
        "matlab_non_lactating_relative_half_range": 0.05,
    }


def log_gaussian_kernel(eval_values: np.ndarray, centers: np.ndarray, scales: np.ndarray) -> np.ndarray:
    standardized = (eval_values[:, None, :] - centers[None, :, :]) / scales[None, None, :]
    normalizer = np.sum(np.log(scales)) + 0.5 * eval_values.shape[1] * np.log(2.0 * np.pi)
    return -0.5 * np.sum(standardized * standardized, axis=2) - normalizer


def density_ratio_mi(
    theta: np.ndarray,
    measurements: np.ndarray,
    measurement_sigma: np.ndarray,
    repeats: int,
    evaluation_samples: int,
    reference_samples: int,
    block_size: int,
    rng: np.random.Generator,
) -> tuple[float, float, int, int]:
    """Estimate MATLAB-style E[log p(theta,z)-log p(theta)-log p(z)]."""

    theta = np.asarray(theta, dtype=float)
    measurements = np.asarray(measurements, dtype=float)
    theta_bw = robust_bandwidth(theta)
    values = []
    ref_count = min(len(theta), reference_samples) if reference_samples > 0 else len(theta)
    eval_count = min(ref_count, evaluation_samples) if evaluation_samples > 0 else ref_count

    for _ in range(repeats):
        ref_idx = choose_indices(len(theta), reference_samples, rng)
        # MATLAB evaluates samples against a mixture containing those same
        # samples. Draw evaluation rows from the reference ensemble to preserve
        # that estimator behavior when runtime-controlling subsamples are used.
        if evaluation_samples <= 0 or evaluation_samples >= len(ref_idx):
            eval_idx = ref_idx
        else:
            eval_idx = np.sort(rng.choice(ref_idx, size=evaluation_samples, replace=False))
        theta_eval = theta[eval_idx]
        theta_ref = theta[ref_idx]
        y_eval = measurements[eval_idx]
        y_ref = measurements[ref_idx]
        z_eval = np.maximum(
            y_eval + rng.normal(0.0, measurement_sigma, size=y_eval.shape),
            0.0,
        )
        log_ratios = []

        for start in range(0, len(eval_idx), block_size):
            stop = min(start + block_size, len(eval_idx))
            log_theta = log_gaussian_kernel(theta_eval[start:stop], theta_ref, theta_bw)
            log_z = log_gaussian_kernel(z_eval[start:stop], y_ref, measurement_sigma)
            log_p_theta = logsumexp(log_theta, axis=1) - np.log(len(ref_idx))
            log_p_z = logsumexp(log_z, axis=1) - np.log(len(ref_idx))
            log_p_joint = logsumexp(log_theta + log_z, axis=1) - np.log(len(ref_idx))
            log_ratios.append(log_p_joint - log_p_theta - log_p_z)

        values.append(float(np.mean(np.concatenate(log_ratios))))

    return float(np.mean(values)), float(np.std(values, ddof=1)) if repeats > 1 else 0.0, eval_count, ref_count


def global_day_ranking(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    species: list[str],
    days: list[int],
    noise: dict[str, float],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    theta = parameters.to_numpy(dtype=float)
    sigma = np.array([noise[name] for name in species], dtype=float)
    for position, day in enumerate(days, start=1):
        y = outputs[day_columns(day, species)].to_numpy(dtype=float)
        mean, sd, n_eval, n_ref = density_ratio_mi(
            theta, y, sigma, args.mi_repeats, args.evaluation_samples,
            args.reference_samples, args.block_size, rng,
        )
        rows.append(
            {
                "day": day,
                "mi_mean_nats": mean,
                "mi_sd_nats": sd,
                "n_parameters": theta.shape[1],
                "n_biomarkers": len(species),
                "evaluation_samples": n_eval,
                "reference_samples": n_ref,
            }
        )
        print(f"Global I(Theta;Y_day): {position}/{len(days)} day={day} MI={mean:.6g}")
    return pd.DataFrame(rows).sort_values("mi_mean_nats", ascending=False).reset_index(drop=True)


def species_ranking(
    theta: np.ndarray,
    outputs: pd.DataFrame,
    species: list[str],
    day: int,
    noise: dict[str, float],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    for name in species:
        y = outputs[[f"{name}_day_{day}"]].to_numpy(dtype=float)
        sigma = np.array([noise[name]], dtype=float)
        mean, sd, n_eval, n_ref = density_ratio_mi(
            theta, y, sigma, args.mi_repeats, args.evaluation_samples,
            args.reference_samples, args.block_size, rng,
        )
        rows.append(
            {
                "biomarker": name,
                "day": day,
                "mi_mean_nats": mean,
                "mi_sd_nats": sd,
                "evaluation_samples": n_eval,
                "reference_samples": n_ref,
            }
        )
    return pd.DataFrame(rows).sort_values("mi_mean_nats", ascending=False).reset_index(drop=True)


def parameter_day_rankings(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    targets: list[str],
    species: list[str],
    days: list[int],
    noise: dict[str, float],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    sigma = np.array([noise[name] for name in species], dtype=float)
    for target in targets:
        theta = parameters[[target]].to_numpy(dtype=float)
        for day in days:
            y = outputs[day_columns(day, species)].to_numpy(dtype=float)
            mean, sd, n_eval, n_ref = density_ratio_mi(
                theta, y, sigma, args.mi_repeats, args.evaluation_samples,
                args.reference_samples, args.block_size, rng,
            )
            rows.append(
                {
                    "parameter": target,
                    "identifiability_class": TARGET_CLASSES.get(target, "Selected"),
                    "day": day,
                    "mi_mean_nats": mean,
                    "mi_sd_nats": sd,
                    "evaluation_samples": n_eval,
                    "reference_samples": n_ref,
                }
            )
        best = max(rows[-len(days):], key=lambda row: row["mi_mean_nats"])
        print(f"Parameter I(theta;Y_day): {target} best day={best['day']} MI={best['mi_mean_nats']:.6g}")
    return pd.DataFrame(rows)


def cumulative_biomarker_rankings(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    targets: list[str],
    biomarker_order: list[str],
    selected_day: int,
    noise: dict[str, float],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Evaluate parameter-specific MI using the global thesis acquisition order."""

    rows = []
    for target in targets:
        for count in range(1, len(biomarker_order) + 1):
            selected = biomarker_order[:count]
            y = outputs[day_columns(selected_day, selected)].to_numpy(dtype=float)
            sigma = np.array([noise[name] for name in selected], dtype=float)
            mean, sd, n_eval, n_ref = density_ratio_mi(
                parameters[[target]].to_numpy(dtype=float), y, sigma,
                args.mi_repeats, args.evaluation_samples, args.reference_samples,
                args.block_size, rng,
            )
            rows.append(
                {
                    "parameter": target,
                    "identifiability_class": TARGET_CLASSES.get(target, "Selected"),
                    "global_best_day": selected_day,
                    "n_biomarkers": count,
                    "biomarkers": "|".join(selected),
                    "mi_mean_nats": mean,
                    "mi_sd_nats": sd,
                    "evaluation_samples": n_eval,
                    "reference_samples": n_ref,
                }
            )
    return pd.DataFrame(rows)


def weighted_kde_curve(values: np.ndarray, weights: np.ndarray, grid: np.ndarray) -> np.ndarray:
    kde = gaussian_kde(values, weights=weights)
    density = kde(grid)
    area = np.trapezoid(density, grid)
    return density / area if area > 0 else density


def posterior_curves(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    biomarker_order: list[str],
    best_day: int,
    lowest_day: int | None,
    noise: dict[str, float],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    curve_rows = []
    diagnostic_rows = []
    if not 0 <= args.posterior_observation_index < len(outputs):
        raise ValueError(
            f"--posterior-observation-index must be between 0 and {len(outputs) - 1}."
        )
    day_contexts = [("highest_information_day", best_day)]
    if lowest_day is not None:
        day_contexts.append(("lowest_information_day", lowest_day))
    for day_role, selected_day in day_contexts:
        for target in args.target_parameters:
            values = parameters[target].to_numpy(dtype=float)
            value_span = float(values.max() - values.min())
            grid_margin = max(args.posterior_grid_margin * value_span, 1e-12)
            grid = np.linspace(
                values.min() - grid_margin,
                values.max() + grid_margin,
                args.posterior_grid_size,
            )
            prior_weights = np.full(len(values), 1.0 / len(values))
            prior_density = weighted_kde_curve(values, prior_weights, grid)
            prior_variance = float(np.var(values, ddof=1))
            for x, density in zip(grid, prior_density):
                curve_rows.append(
                    {
                        "parameter": target, "day_role": day_role, "day": selected_day,
                        "n_biomarkers": 0, "biomarkers": "Prior", "grid": x, "density": density,
                    }
                )
            diagnostic_rows.append(
                {
                    "parameter": target, "identifiability_class": TARGET_CLASSES.get(target, "Selected"),
                    "day_role": day_role, "day": selected_day, "n_biomarkers": 0,
                    "biomarkers": "Prior", "ess": len(values), "ess_fraction": 1.0,
                    "normalization_integral": float(np.trapezoid(prior_density, grid)),
                    "posterior_to_prior_variance": 1.0,
                    "observation_index": args.posterior_observation_index,
                }
            )

            for count in range(1, len(biomarker_order) + 1):
                selected = biomarker_order[:count]
                columns = day_columns(selected_day, selected)
                y = outputs[columns].to_numpy(dtype=float)
                z = y[args.posterior_observation_index]
                sigma = np.array([noise[name] for name in selected], dtype=float)
                # Discrete form of p(theta_i,z*) followed by division by
                # integral p(theta_i,z*) dtheta_i.
                log_joint_weights = -0.5 * np.sum(((z[None, :] - y) / sigma[None, :]) ** 2, axis=1)
                log_normalizer = logsumexp(log_joint_weights)
                weights = np.exp(log_joint_weights - log_normalizer)
                density = weighted_kde_curve(values, weights, grid)
                normalization_integral = float(np.trapezoid(density, grid))
                ess = float(1.0 / np.sum(weights * weights))
                weighted_mean = float(np.sum(weights * values))
                posterior_variance = float(np.sum(weights * (values - weighted_mean) ** 2))
                label = "|".join(selected)
                for x, density_value in zip(grid, density):
                    curve_rows.append(
                        {
                            "parameter": target, "day_role": day_role, "day": selected_day,
                            "n_biomarkers": count, "biomarkers": label,
                            "grid": x, "density": density_value,
                        }
                    )
                diagnostic_rows.append(
                    {
                        "parameter": target,
                        "identifiability_class": TARGET_CLASSES.get(target, "Selected"),
                        "day_role": day_role, "day": selected_day, "n_biomarkers": count,
                        "biomarkers": label, "ess": ess, "ess_fraction": ess / len(values),
                        "normalization_integral": normalization_integral,
                        "posterior_to_prior_variance": posterior_variance / prior_variance,
                        "observation_index": args.posterior_observation_index,
                    }
                )
    return pd.DataFrame(curve_rows), pd.DataFrame(diagnostic_rows)


def plot_outputs(
    output_dir: Path,
    global_ranking: pd.DataFrame | None,
    global_species: pd.DataFrame,
    parameter_days: pd.DataFrame,
    cumulative: pd.DataFrame,
    curves: pd.DataFrame,
    n_admissible: int,
    analysis_label: str,
) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    sample_source = global_ranking if global_ranking is not None else global_species
    n_evaluated = int(sample_source["evaluation_samples"].iloc[0])
    n_reference = int(sample_source["reference_samples"].iloc[0])
    sample_note = (
        f"{analysis_label} | admissible bank: {n_admissible:,}; "
        f"MI evaluated: {n_evaluated:,}; mixture reference: {n_reference:,}"
    )

    if global_ranking is not None:
        ordered = global_ranking.sort_values("day")
        fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
        ax.bar(ordered["day"], ordered["mi_mean_nats"], color="tab:blue")
        ax.errorbar(ordered["day"], ordered["mi_mean_nats"], yerr=ordered["mi_sd_nats"], fmt="none", color="0.2", capsize=2)
        best = global_ranking.iloc[0]
        worst = global_ranking.iloc[-1]
        ax.bar(int(best["day"]), float(best["mi_mean_nats"]), color="black", label=f"Highest: day {int(best['day'])}")
        ax.bar(int(worst["day"]), float(worst["mi_mean_nats"]), color="tab:red", label=f"Lowest: day {int(worst['day'])}")
        day_grid = np.linspace(float(ordered["day"].min()), float(ordered["day"].max()), 500)
        nonnegative_mi = np.maximum(ordered["mi_mean_nats"].to_numpy(dtype=float), 0.0)
        if np.any(nonnegative_mi > 0):
            kde_guide = gaussian_kde(ordered["day"].to_numpy(dtype=float), weights=nonnegative_mi)
            smoothed = kde_guide(day_grid)
            smoothed *= float(nonnegative_mi.max()) / float(smoothed.max())
            ax.plot(
                day_grid,
                smoothed,
                color="tab:orange",
                linewidth=2.0,
                label="KDE-smoothed MI-by-day guide",
            )
        ax.set(
            title="Global BED: information between full parameter vector and all biomarkers",
            xlabel="Candidate sampling day",
            ylabel="Estimated I(Theta; Y_day), nats",
        )
        ax.text(0.01, 0.99, sample_note, transform=ax.transAxes, va="top", fontsize=8, color="0.25")
        ax.legend()
        ax.grid(axis="y", alpha=0.25)
        fig.savefig(figure_dir / "global_full_parameter_mi_by_day.png", dpi=220)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ranked = global_species.sort_values("mi_mean_nats")
    ax.barh(ranked["biomarker"], ranked["mi_mean_nats"], xerr=ranked["mi_sd_nats"], color="tab:blue")
    ax.set(title=f"Biomarker information about full parameter vector at day {int(global_species['day'].iloc[0])}", xlabel="Estimated I(Theta; Y_i), nats", ylabel="Biomarker")
    ax.text(0.99, 0.01, sample_note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="0.25")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(figure_dir / "global_species_ranking_at_best_day.png", dpi=220)
    plt.close(fig)

    targets = parameter_days["parameter"].unique().tolist()
    fig, axes = plt.subplots(3, 3, figsize=(16, 12), constrained_layout=True)
    for ax, target in zip(axes.ravel(), targets):
        subset = parameter_days[parameter_days["parameter"] == target].sort_values("day")
        ax.plot(subset["day"], subset["mi_mean_nats"], marker="o", linewidth=1.5)
        ax.fill_between(subset["day"], subset["mi_mean_nats"] - subset["mi_sd_nats"], subset["mi_mean_nats"] + subset["mi_sd_nats"], alpha=0.2)
        ax.set_title(f"{target}\n{TARGET_CLASSES.get(target, 'Selected')}", fontsize=10)
        ax.set_xlabel("Day")
        ax.set_ylabel("I(theta; Y_day), nats")
        ax.grid(alpha=0.25)
    for ax in axes.ravel()[len(targets):]:
        ax.axis("off")
    fig.suptitle("Parameter-specific BED: all biomarkers at each single day", fontsize=15)
    fig.savefig(figure_dir / "parameter_specific_mi_by_day_3x3.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(16, 12), constrained_layout=True)
    for ax, target in zip(axes.ravel(), targets):
        subset = cumulative[cumulative["parameter"] == target].sort_values("n_biomarkers")
        ax.plot(subset["n_biomarkers"], subset["mi_mean_nats"], marker="o")
        cumulative_labels = []
        previous: list[str] = []
        for names in subset["biomarkers"]:
            current = str(names).split("|")
            added = current[-1] if len(current) > len(previous) else current[-1]
            cumulative_labels.append(f"{len(current)}: +{added}" if len(current) > 1 else f"1: {added}")
            previous = current
        ax.set_title(f"{target}\nglobal best day {int(subset['global_best_day'].iloc[0])}", fontsize=10)
        ax.set_xlabel("Cumulative biomarkers in global MI order")
        ax.set_ylabel("I(theta; Y_1:k), nats")
        ax.set_xticks(subset["n_biomarkers"])
        ax.set_xticklabels(cumulative_labels, rotation=35, ha="right", fontsize=7)
        ax.grid(alpha=0.25)
    for ax in axes.ravel()[len(targets):]:
        ax.axis("off")
    fig.suptitle(
        "Cumulative biomarker information using global full-vector MI ranking\n" + sample_note,
        fontsize=15,
    )
    fig.savefig(figure_dir / "cumulative_biomarker_mi_3x3.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(16, 12), constrained_layout=True)
    for ax, target in zip(axes.ravel(), targets):
        subset = curves[
            (curves["parameter"] == target)
            & (curves["day_role"] == "highest_information_day")
        ]
        counts = sorted(subset["n_biomarkers"].unique())
        colors = plt.cm.tab10(np.linspace(0, 1, len(counts)))
        line_styles = ["-", "--", "-.", ":", "--", "-.", ":", "--", "-"]
        markers = [None, "o", "s", "^", "D", "v", "P", "X", "*"]
        for position, (color, count) in enumerate(zip(colors, counts)):
            curve = subset[subset["n_biomarkers"] == count]
            if count == 0:
                label = "Prior"
            else:
                biomarkers = str(curve["biomarkers"].iloc[0]).replace("|", " + ")
                label = f"Best {count}: {biomarkers}"
            ax.plot(
                curve["grid"],
                curve["density"],
                linewidth=1.8 if count == 0 else 1.2,
                color="black" if count == 0 else color,
                linestyle=line_styles[position % len(line_styles)],
                marker=markers[position % len(markers)],
                markevery=(position * 4, 55) if count > 0 else None,
                markersize=3.2,
                markerfacecolor="none",
                zorder=len(counts) - position,
                label=label,
            )
        ax.set_title(target, fontsize=10)
        ax.set_xlabel("Parameter value")
        ax.set_ylabel("Density")
        ax.grid(alpha=0.25)
    axes.ravel()[0].legend(fontsize=5.5)
    for ax in axes.ravel()[len(targets):]:
        ax.axis("off")
    fig.suptitle(
        "Prior-to-posterior narrowing with cumulative biomarkers\n"
        + sample_note
        + "\nMarkers distinguish overlapping cumulative posteriors",
        fontsize=15,
    )
    fig.savefig(figure_dir / "posterior_narrowing_cumulative_biomarkers_3x3.png", dpi=220)
    plt.close(fig)

    if "lowest_information_day" in set(curves["day_role"]):
        fig, axes = plt.subplots(3, 3, figsize=(16, 12), constrained_layout=True)
        for ax, target in zip(axes.ravel(), targets):
            subset = curves[curves["parameter"] == target]
            for role, style in [("highest_information_day", "-"), ("lowest_information_day", "--")]:
                role_curves = subset[subset["day_role"] == role]
                max_count = int(role_curves["n_biomarkers"].max())
                curve = role_curves[role_curves["n_biomarkers"] == max_count]
                ax.plot(curve["grid"], curve["density"], style, linewidth=1.8, label=role.replace("_", " "))
            prior = subset[
                (subset["day_role"] == "highest_information_day") & (subset["n_biomarkers"] == 0)
            ]
            ax.plot(prior["grid"], prior["density"], color="0.5", linewidth=1.2, label="prior")
            ax.set_title(target, fontsize=10)
            ax.set_xlabel("Parameter value")
            ax.set_ylabel("Density")
            ax.grid(alpha=0.25)
        axes.ravel()[0].legend(fontsize=8)
        fig.suptitle(
            "Posterior comparison: highest- versus lowest-information day\n" + sample_note,
            fontsize=15,
        )
        fig.savefig(figure_dir / "posterior_highest_vs_lowest_information_day_3x3.png", dpi=220)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    parameters, outputs, nominal, days, missing_days = load_inputs(args)
    noise = load_measurement_noise(args.measurement_noise_csv, args.species, args.noise_floor)
    prior_audit = prior_range_audit(parameters)
    observed_half_range = prior_audit["relative_half_range_median"]
    if abs(observed_half_range - args.required_prior_half_range) > args.prior_range_tolerance:
        raise ValueError(
            "Stored simulation bank does not match the required PhD prior: "
            f"observed median relative half-range={observed_half_range:.6g}, "
            f"required={args.required_prior_half_range:.6g}. Generate the +/-5% bank first."
        )
    rng = np.random.default_rng(args.seed)

    audit = {
        "ode_simulations_run": 0,
        "admissible_samples": len(parameters),
        "parameter_count": len(parameters.columns),
        "species": args.species,
        "requested_days": list(range(args.day_start, args.day_end + 1)),
        "used_days": days,
        "missing_requested_days": missing_days,
        "measurement_sigma": noise,
        "evaluation_samples": args.evaluation_samples,
        "reference_samples": args.reference_samples,
        "mi_repeats": args.mi_repeats,
        "estimator": "MATLAB-style Gaussian-mixture density ratio, blockwise",
        "stored_prior_range_audit": prior_audit,
        "required_prior_half_range": args.required_prior_half_range,
        "analysis_label": args.analysis_label,
        "posterior_grid_margin": args.posterior_grid_margin,
    }
    (args.output_dir / "input_audit.json").write_text(json.dumps(audit, indent=2))

    if args.skip_global:
        if args.best_day is None:
            raise ValueError("--best-day is required with --skip-global.")
        if args.best_day not in days:
            raise ValueError(f"--best-day {args.best_day} is not available. Available: {days}")
        global_ranking = None
        best_day = args.best_day
        lowest_day = None
    else:
        global_ranking = global_day_ranking(parameters, outputs, args.species, days, noise, args, rng)
        global_ranking.to_csv(args.output_dir / "global_full_parameter_mi_by_day.csv", index=False)
        best_day = int(global_ranking.iloc[0]["day"])
        lowest_day = int(global_ranking.iloc[-1]["day"])

    global_species = species_ranking(
        parameters.to_numpy(dtype=float), outputs, args.species, best_day, noise, args, rng
    )
    global_species.to_csv(args.output_dir / "global_species_ranking_at_best_day.csv", index=False)
    biomarker_order = global_species["biomarker"].tolist()
    pd.DataFrame(
        {
            "rank": range(1, len(biomarker_order) + 1),
            "biomarker": biomarker_order,
            "global_best_day": best_day,
        }
    ).to_csv(args.output_dir / "global_cumulative_biomarker_order.csv", index=False)

    parameter_days = parameter_day_rankings(
        parameters, outputs, args.target_parameters, args.species, days, noise, args, rng
    )
    parameter_days.to_csv(args.output_dir / "parameter_specific_mi_by_day.csv", index=False)

    cumulative = cumulative_biomarker_rankings(
        parameters, outputs, args.target_parameters, biomarker_order, best_day, noise, args, rng
    )
    cumulative.to_csv(args.output_dir / "cumulative_biomarker_mi.csv", index=False)

    curves, diagnostics = posterior_curves(
        parameters, outputs, biomarker_order, best_day, lowest_day, noise, args
    )
    curves.to_csv(args.output_dir / "posterior_curves.csv", index=False)
    diagnostics.to_csv(args.output_dir / "posterior_diagnostics.csv", index=False)
    final_learning = (
        diagnostics[
            (diagnostics["day_role"] == "highest_information_day")
            & (diagnostics["n_biomarkers"] == len(biomarker_order))
        ][
            [
                "parameter",
                "identifiability_class",
                "day",
                "biomarkers",
                "ess",
                "ess_fraction",
                "normalization_integral",
                "posterior_to_prior_variance",
            ]
        ]
        .sort_values(["identifiability_class", "posterior_to_prior_variance"])
        .reset_index(drop=True)
    )
    final_learning.to_csv(args.output_dir / "posterior_learning_by_identifiability_class.csv", index=False)
    max_normalization_error = float(
        np.max(np.abs(diagnostics["normalization_integral"].to_numpy(dtype=float) - 1.0))
    )
    if max_normalization_error > 1e-3:
        raise ValueError(
            "Posterior normalization verification failed: "
            f"maximum absolute integral error={max_normalization_error:.6g}"
        )
    plot_outputs(
        args.output_dir,
        global_ranking,
        global_species,
        parameter_days,
        cumulative,
        curves,
        len(parameters),
        args.analysis_label,
    )

    print(f"Saved MATLAB-faithful parameter BED outputs to {args.output_dir}")
    print(f"ODE simulations run: 0")
    print(f"Admissible rows used: {len(parameters)}")
    print(f"Available requested days used: {days[0]}-{days[-1]}")
    if missing_days:
        print(f"Unavailable requested days skipped: {missing_days}")
    print(
        "Stored prior median relative half-range: "
        f"{100 * prior_audit['relative_half_range_median']:.3g}% "
        "(MATLAB non-lactating workflow: 5%)"
    )
    print(f"Best global day: {best_day}")
    if lowest_day is not None:
        print(f"Lowest global day: {lowest_day}")
    print(f"Maximum posterior normalization error: {max_normalization_error:.3g}")


if __name__ == "__main__":
    main()
