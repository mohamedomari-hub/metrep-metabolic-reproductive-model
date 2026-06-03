"""Create thesis-style surrogate BED summary figures.

The figure is generated from the documented surrogate BED pilot outputs. It
does not replace the MATLAB thesis result; it gives a Python/surrogate view
with the same scientific structure: nominal trajectories, information by day,
information by species, posterior narrowing, and ODE-versus-surrogate timing.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from sklearn.neighbors import KernelDensity


PROJECT_ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = PROJECT_ROOT / "MetRep_Python"
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from model_definition.initial_conditions import STATE_INDEX
from model_definition.parameters import default_parameters
from model_definition.scenarios import constant_non_lactating
from model_definition.simulate import run_simulation
from surrogate_bed_pipeline import train_surrogate


SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]
SPECIES_LABELS = {
    "FSH": "FSH",
    "PGF": "PGF",
    "P4": "P4",
    "E2": "E2",
    "INH": "INH",
    "IGF1": "IGF",
    "Insulin": "Ins",
    "Glucose": "Glu",
    "Glucagon": "Gluca",
}


def parse_args() -> argparse.Namespace:
    base = Path("docs/Bayesian_Experimental_Design/surrogate_bed")
    parser = argparse.ArgumentParser(description="Create thesis-style surrogate BED plots.")
    parser.add_argument("--input-dir", type=Path, default=base / "input_tables")
    parser.add_argument("--run-dir", type=Path, default=base / "run_outputs")
    parser.add_argument("--candidate-map-csv", type=Path, help="Candidate map used for the BED run.")
    parser.add_argument("--target-column", default="follicle_to_e2_scale")
    parser.add_argument(
        "--measurement-mode",
        choices=["auto", "cumulative", "single-day"],
        default="auto",
        help="Use thesis-style cumulative candidates, single-day candidates, or auto-detect from mi_candidate_ranking.csv.",
    )
    parser.add_argument("--output-dir", type=Path, default=base / "run_outputs" / "figures")
    parser.add_argument("--relative-noise", type=float, default=0.05)
    parser.add_argument("--noise-floor", type=float, default=1e-8)
    parser.add_argument("--synthetic-observation-seed", type=int, default=7)
    parser.add_argument("--posterior-grid-size", type=int, default=500)
    parser.add_argument("--benchmark-samples", type=int, default=12)
    parser.add_argument("--surrogate-benchmark-repeats", type=int, default=200)
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-benchmark", action="store_true")
    return parser.parse_args()


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")
    return pd.read_csv(path)


def load_tables(args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    parameters = read_csv(args.input_dir / "prior_parameter_samples.csv")
    outputs = read_csv(args.input_dir / "ode_output_features.csv")
    admissibility = read_csv(args.input_dir / "admissibility.csv")
    admissible = admissibility["admissible"].astype(bool).to_numpy()
    mi = read_csv(args.run_dir / "mi_candidate_ranking.csv")
    candidate_map_path = args.candidate_map_csv
    if candidate_map_path is None:
        has_cumulative = mi["candidate"].astype(str).str.contains("_cumulative_to_day_").any()
        cumulative_path = args.input_dir / "candidate_map_cumulative.csv"
        single_day_path = args.input_dir / "candidate_map.csv"
        candidate_map_path = cumulative_path if has_cumulative and cumulative_path.exists() else single_day_path
    return {
        "parameters": parameters.loc[admissible].reset_index(drop=True),
        "outputs": outputs.loc[admissible].reset_index(drop=True),
        "nominal": read_csv(args.input_dir / "nominal_output.csv"),
        "candidate_map": read_csv(candidate_map_path),
        "mi": mi,
        "predicted": read_csv(args.run_dir / "surrogate_predicted_outputs.csv"),
    }


def parse_day(candidate: str) -> int | None:
    if "_cumulative_to_day_" in candidate:
        day_text = candidate.rsplit("_cumulative_to_day_", maxsplit=1)[-1]
    elif "_day_" in candidate:
        day_text = candidate.rsplit("_day_", maxsplit=1)[-1]
    else:
        return None
    try:
        return int(float(day_text.replace("p", ".")))
    except ValueError:
        return None


def choose_measurement_mode(mi: pd.DataFrame, requested: str) -> str:
    if requested != "auto":
        return requested
    candidates = mi["candidate"].astype(str)
    if candidates.str.startswith("all_species_cumulative_to_day_").any():
        return "cumulative"
    return "single-day"


def all_species_candidate_name(day: int, mode: str) -> str:
    if mode == "cumulative":
        return f"all_species_cumulative_to_day_{day}"
    return f"all_species_day_{day}"


def species_candidate_name(species: str, day: int, mode: str) -> str:
    if mode == "cumulative":
        return f"{species}_cumulative_to_day_{day}"
    return f"{species}_day_{day}"


def subset_candidate_name(species: list[str], day: int, mode: str) -> str:
    if mode == "cumulative":
        return f"{'_'.join(species)}_cumulative_to_day_{day}"
    return "|".join(f"{state}_day_{day}" for state in species)


def all_species_mi_by_day(mi: pd.DataFrame, mode: str) -> pd.DataFrame:
    rows = []
    for _, row in mi.iterrows():
        candidate = str(row["candidate"])
        if mode == "cumulative" and not candidate.startswith("all_species_cumulative_to_day_"):
            continue
        if mode == "single-day" and not candidate.startswith("all_species_day_"):
            continue
        day = parse_day(candidate)
        if day is not None:
            rows.append({"day": day, "mi_nats": float(row["mi_nats"]), "candidate": candidate})
    return pd.DataFrame(rows).sort_values("day").reset_index(drop=True)


def species_mi_at_day(mi: pd.DataFrame, day: int, mode: str) -> pd.DataFrame:
    rows = []
    for species in SPECIES:
        candidate = species_candidate_name(species, day, mode)
        match = mi[mi["candidate"] == candidate]
        rows.append(
            {
                "species": species,
                "label": SPECIES_LABELS[species],
                "mi_nats": float(match["mi_nats"].iloc[0]) if not match.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def candidate_columns(candidate_map: pd.DataFrame, candidate: str) -> list[str]:
    match = candidate_map[candidate_map["candidate"] == candidate]
    if match.empty:
        if "|" in candidate:
            return [column.strip() for column in candidate.split("|") if column.strip()]
        raise ValueError(f"Candidate not found in candidate map: {candidate}")
    return [column.strip() for column in str(match["columns"].iloc[0]).split("|") if column.strip()]


def nominal_follicle_p4() -> pd.DataFrame:
    scenario = constant_non_lactating(days=90.0, step=0.25)
    result = run_simulation(scenario, method="BDF", rtol=1e-6, atol=1e-9)
    return pd.DataFrame(
        {
            "day": result.t,
            "Follicle": result.y[:, STATE_INDEX["Follicle"]],
            "P4": result.y[:, STATE_INDEX["P4"]],
        }
    )


def kde_curve(values: np.ndarray, weights: np.ndarray | None, grid: np.ndarray) -> np.ndarray:
    value_range = float(np.max(values) - np.min(values))
    bandwidth = max(value_range / 35.0, 1e-8)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    if weights is None:
        kde.fit(values[:, None])
    else:
        kde.fit(values[:, None], sample_weight=weights)
    density = np.exp(kde.score_samples(grid[:, None]))
    area = np.trapezoid(density, grid)
    return density / area if area > 0 else density


def posterior_weights(predicted: pd.DataFrame, observation: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    residual = (predicted.to_numpy(dtype=float) - observation[None, :]) / sigma[None, :]
    log_likelihood = -0.5 * np.sum(residual * residual, axis=1)
    shifted = log_likelihood - float(np.max(log_likelihood))
    weights = np.exp(shifted)
    total = float(np.sum(weights))
    if total <= 0 or not np.isfinite(total):
        raise ValueError("Posterior weights collapsed for one plotted candidate.")
    return weights / total


def synthetic_observation_for_columns(
    nominal: pd.DataFrame,
    columns: list[str],
    relative_noise: float,
    noise_floor: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    all_columns = list(nominal.columns)
    nominal_values = nominal.loc[nominal.index[0], all_columns].to_numpy(dtype=float)
    sigma_all = np.maximum(relative_noise * np.abs(nominal_values), noise_floor)
    rng = np.random.default_rng(seed)
    synthetic_all = nominal_values + sigma_all * rng.normal(size=len(all_columns))
    index = [all_columns.index(column) for column in columns]
    return synthetic_all[index], sigma_all[index]


def posterior_curves(
    parameters: pd.DataFrame,
    predicted: pd.DataFrame,
    nominal: pd.DataFrame,
    candidate_map: pd.DataFrame,
    target_column: str,
    best_day: int,
    low_day: int,
    mode: str,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target = parameters[target_column].to_numpy(dtype=float)
    pad = 0.06 * (float(np.max(target)) - float(np.min(target)))
    grid = np.linspace(float(np.min(target)) - pad, float(np.max(target)) + pad, args.posterior_grid_size)

    day_label = "cumulative to day" if mode == "cumulative" else "day"
    candidate_specs: list[tuple[str, str | None]] = [
        ("Prior", None),
        (f"All species, {day_label} {best_day}", all_species_candidate_name(best_day, mode)),
        (f"All species, {day_label} {low_day}", all_species_candidate_name(low_day, mode)),
        (
            f"PGF+E2+FSH+INH, {day_label} {best_day}",
            subset_candidate_name(["PGF", "E2", "FSH", "INH"], best_day, mode),
        ),
        (
            f"PGF+E2+FSH, {day_label} {best_day}",
            subset_candidate_name(["PGF", "E2", "FSH"], best_day, mode),
        ),
        (f"PGF+E2, {day_label} {best_day}", subset_candidate_name(["PGF", "E2"], best_day, mode)),
        (f"PGF, {day_label} {best_day}", species_candidate_name("PGF", best_day, mode)),
    ]
    designs = [
        (label, [] if candidate is None else candidate_columns(candidate_map, candidate))
        for label, candidate in candidate_specs
    ]

    curve_rows = []
    diagnostic_rows = []
    for label, columns in designs:
        if not columns:
            density = kde_curve(target, None, grid)
            ess = float(len(target))
        else:
            observation, sigma = synthetic_observation_for_columns(
                nominal,
                columns,
                args.relative_noise,
                args.noise_floor,
                args.synthetic_observation_seed,
            )
            weights = posterior_weights(predicted[columns], observation, sigma)
            density = kde_curve(target, weights, grid)
            ess = float(1.0 / np.sum(weights * weights))
        curve_rows.extend({"target_value": x, "density": y, "curve": label} for x, y in zip(grid, density))
        diagnostic_rows.append({"curve": label, "n_columns": len(columns), "effective_sample_size": ess})
    return pd.DataFrame(curve_rows), pd.DataFrame(diagnostic_rows)


def benchmark_ode_vs_surrogate(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    n = min(args.benchmark_samples, len(parameters))
    x_subset = parameters.iloc[:n].copy()

    surrogate_args = SimpleNamespace(
        regressor="extra_trees",
        n_estimators=args.n_estimators,
        min_samples_leaf=2,
        pca_variance=0.995,
        seed=args.seed,
    )
    train_start = time.perf_counter()
    surrogate = train_surrogate(parameters, outputs, surrogate_args)
    train_seconds = time.perf_counter() - train_start

    predict_start = time.perf_counter()
    for _ in range(args.surrogate_benchmark_repeats):
        surrogate.predict(x_subset)
    surrogate_total = time.perf_counter() - predict_start
    surrogate_seconds = surrogate_total / args.surrogate_benchmark_repeats

    scenario = constant_non_lactating(days=90.0, step=1.0)
    nominal = default_parameters()
    ode_start = time.perf_counter()
    for _, row in x_subset.iterrows():
        params = dict(nominal)
        params.update({column: float(row[column]) for column in parameters.columns})
        run_simulation(scenario, parameters=params, method="BDF", rtol=1e-6, atol=1e-9)
    ode_seconds = time.perf_counter() - ode_start

    ode_per_sample = ode_seconds / n
    surrogate_per_sample = surrogate_seconds / n
    speedup = ode_per_sample / surrogate_per_sample if surrogate_per_sample > 0 else np.inf
    return pd.DataFrame(
        [
            {
                "method": "ODE",
                "n_samples": n,
                "seconds_total": ode_seconds,
                "seconds_per_sample": ode_per_sample,
                "speedup_vs_surrogate_predict": 1.0,
            },
            {
                "method": "Surrogate prediction",
                "n_samples": n,
                "seconds_total": surrogate_seconds,
                "seconds_per_sample": surrogate_per_sample,
                "speedup_vs_surrogate_predict": speedup,
            },
            {
                "method": "Surrogate training",
                "n_samples": len(parameters),
                "seconds_total": train_seconds,
                "seconds_per_sample": train_seconds / len(parameters),
                "speedup_vs_surrogate_predict": np.nan,
            },
        ]
    )


def plot_speed(speed: pd.DataFrame, output_path: Path) -> None:
    plot_rows = speed[speed["method"].isin(["ODE", "Surrogate prediction"])].copy()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.bar(plot_rows["method"], plot_rows["seconds_per_sample"], color=["#3b6ea8", "#d95f02"])
    ax.set_yscale("log")
    ax.set_ylabel("Seconds per parameter sample (log scale)")
    ax.set_title("ODE simulation is the bottleneck; surrogate prediction is faster")
    speedup = float(plot_rows.loc[plot_rows["method"] == "Surrogate prediction", "speedup_vs_surrogate_predict"].iloc[0])
    ax.text(
        0.5,
        0.92,
        f"Measured local speedup: {speedup:,.0f}x per sample",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "0.8", "boxstyle": "round,pad=0.3"},
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_composite(
    trajectories: pd.DataFrame,
    mi_day: pd.DataFrame,
    mi_species: pd.DataFrame,
    posterior: pd.DataFrame,
    speed: pd.DataFrame | None,
    output_path: Path,
) -> None:
    best = mi_day.loc[mi_day["mi_nats"].idxmax()]
    low = mi_day.loc[mi_day["mi_nats"].idxmin()]

    fig = plt.figure(figsize=(15, 10.5))
    grid = GridSpec(3, 2, figure=fig, height_ratios=[0.95, 1.05, 1.35], width_ratios=[1.05, 1.25])
    ax_traj = fig.add_subplot(grid[0, :])
    ax_day = fig.add_subplot(grid[1, :])
    ax_species = fig.add_subplot(grid[2, 0])
    ax_post = fig.add_subplot(grid[2, 1])

    window = trajectories[(trajectories["day"] >= 55) & (trajectories["day"] <= 90)].copy()
    follicle = window["Follicle"] / max(float(window["Follicle"].max()), 1e-12)
    p4 = window["P4"] / max(float(window["P4"].max()), 1e-12)
    ax_traj.plot(window["day"], follicle, color="#1f77b4", linewidth=2.0, label="Follicle")
    ax_traj.plot(window["day"], p4, color="black", linestyle="--", linewidth=1.7, label="P4")
    ax_traj.annotate("", xy=(64, 1.08), xytext=(59, 1.08), arrowprops={"arrowstyle": "<->", "color": "black"})
    ax_traj.text(61.5, 1.15, "Follicular phase", ha="center", fontsize=9)
    ax_traj.annotate("", xy=(81, 1.08), xytext=(65, 1.08), arrowprops={"arrowstyle": "<->", "color": "black"})
    ax_traj.text(73, 1.15, "Luteal phase", ha="center", fontsize=9)
    ax_traj.set_xlim(55, 90)
    ax_traj.set_ylim(0, 1.25)
    ax_traj.set_ylabel("Relative level")
    ax_traj.set_title("Nominal reproductive trajectory used to place candidate measurement days")
    ax_traj.legend(loc="upper right")
    ax_traj.grid(alpha=0.25)

    colors = np.full(len(mi_day), "#1f77b4", dtype=object)
    colors[mi_day["day"].to_numpy() == int(best["day"])] = "black"
    colors[mi_day["day"].to_numpy() == int(low["day"])] = "#d62728"
    ax_day.bar(mi_day["day"], mi_day["mi_nats"], color=colors, edgecolor="white", linewidth=0.4)
    ax_day.set_xlim(54.3, 90.7)
    ax_day.set_ylabel("Estimated MI (nats)")
    ax_day.set_xlabel("Candidate sampling day")
    ax_day.set_title("Experiment choice by expected information about the target parameter")
    ax_day.axvline(float(best["day"]), color="black", linewidth=1.2)
    ax_day.axvline(float(low["day"]), color="#d62728", linewidth=1.2, linestyle="--")
    ax_day.text(float(best["day"]) + 0.4, float(best["mi_nats"]), f"best day {int(best['day'])}", va="bottom", fontsize=9)
    ax_day.text(float(low["day"]) + 0.4, float(low["mi_nats"]), f"lowest day {int(low['day'])}", va="bottom", fontsize=9, color="#d62728")
    ax_day.grid(axis="y", alpha=0.25)

    species_colors = ["black" if np.isfinite(v) else "0.7" for v in mi_species["mi_nats"]]
    ax_species.bar(mi_species["label"], mi_species["mi_nats"], color=species_colors)
    ax_species.set_ylabel("Estimated MI (nats)")
    ax_species.set_xlabel("Measured species")
    ax_species.set_title(f"Per-species information at day {int(best['day'])}")
    ax_species.tick_params(axis="x", rotation=35)
    ax_species.grid(axis="y", alpha=0.25)

    line_styles = {
        "Prior": {"color": "#1f77b4", "linestyle": "--", "linewidth": 2.2},
        f"All species, day {int(best['day'])}": {"color": "black", "linestyle": "-", "linewidth": 2.0},
        f"All species, day {int(low['day'])}": {"color": "#d62728", "linestyle": "-", "linewidth": 1.8},
        f"All species, cumulative to day {int(best['day'])}": {"color": "black", "linestyle": "-", "linewidth": 2.0},
        f"All species, cumulative to day {int(low['day'])}": {"color": "#d62728", "linestyle": "-", "linewidth": 1.8},
    }
    fallback_colors = ["#ff7f0e", "#8c564b", "#2ca02c", "#bcbd22", "#9467bd"]
    for index, (curve, subset) in enumerate(posterior.groupby("curve", sort=False)):
        style = line_styles.get(curve, {"color": fallback_colors[index % len(fallback_colors)], "linestyle": "-", "linewidth": 1.6})
        ax_post.plot(subset["target_value"], subset["density"], label=curve, **style)
    ax_post.set_xlabel("Target parameter value")
    ax_post.set_ylabel("Probability density")
    ax_post.set_title("Prior updating by measuring different species/day choices")
    ax_post.legend(fontsize=8, loc="upper left")
    ax_post.grid(alpha=0.25)

    if speed is not None:
        speed_rows = speed[speed["method"].isin(["ODE", "Surrogate prediction"])]
        speedup = float(speed_rows.loc[speed_rows["method"] == "Surrogate prediction", "speedup_vs_surrogate_predict"].iloc[0])
        ax_post.text(
            0.98,
            0.96,
            f"ODE vs surrogate prediction:\n{speedup:,.0f}x faster per sample",
            transform=ax_post.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "0.8", "boxstyle": "round,pad=0.35"},
        )

    fig.suptitle("Surrogate-assisted experiment choice for parameter inference", fontsize=18, fontweight="bold", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables = load_tables(args)

    mode = choose_measurement_mode(tables["mi"], args.measurement_mode)
    mi_day = all_species_mi_by_day(tables["mi"], mode)
    if mi_day.empty:
        raise ValueError(f"No all-species {mode} candidates found in mi_candidate_ranking.csv")
    best_day = int(mi_day.loc[mi_day["mi_nats"].idxmax(), "day"])
    low_day = int(mi_day.loc[mi_day["mi_nats"].idxmin(), "day"])
    mi_species = species_mi_at_day(tables["mi"], best_day, mode)
    posterior, posterior_diag = posterior_curves(
        tables["parameters"],
        tables["predicted"],
        tables["nominal"],
        tables["candidate_map"],
        args.target_column,
        best_day,
        low_day,
        mode,
        args,
    )
    trajectories = nominal_follicle_p4()

    speed = None
    if not args.skip_benchmark:
        speed = benchmark_ode_vs_surrogate(tables["parameters"], tables["outputs"], args)
        speed.to_csv(args.run_dir / "thesis_style_speed_benchmark.csv", index=False)
        plot_speed(speed, args.output_dir / "surrogate_bed_ode_vs_surrogate_speed.png")

    posterior.to_csv(args.run_dir / "thesis_style_posterior_curves.csv", index=False)
    posterior_diag.to_csv(args.run_dir / "thesis_style_posterior_diagnostics.csv", index=False)
    plot_composite(
        trajectories,
        mi_day,
        mi_species,
        posterior,
        speed,
        args.output_dir / "surrogate_bed_thesis_style_summary.png",
    )
    print(f"Saved thesis-style surrogate BED figures to {args.output_dir}")
    print(f"Measurement mode: {mode}")
    print(f"Best all-species day: {best_day}; lowest all-species day: {low_day}")
    if speed is not None:
        speedup = float(speed.loc[speed["method"] == "Surrogate prediction", "speedup_vs_surrogate_predict"].iloc[0])
        print(f"Measured ODE/surrogate prediction speedup: {speedup:,.0f}x per sample")


if __name__ == "__main__":
    main()
