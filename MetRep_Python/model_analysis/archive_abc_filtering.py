"""Archive-based sequential ABC filtering using the saved broad +/-5% ODE archive.

No new ODE simulations are run here. The script implements a proper
likelihood-free sequential thresholding workflow over real ODE archive rows for
the fixed 3x3 representative parameter set and the same BED observation
scenarios used by posterior reweighting and reduced archive posterior updates.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parents[2]
BROAD_BANK = ROOT / "local_data/phd_bed_bank_5pct_50k_glucagon"
DEFAULT_BED = ROOT / "local_outputs/bed_targeted"
DEFAULT_GSA_UNCERTAINTY_SCENARIOS = ROOT / "results_final/tables/bed_targeted_gsa_uncertainty_observation_scenarios.csv"
DEFAULT_BED_GUIDED_SCENARIOS = ROOT / "results_final/tables/bed_guided_selected_observation_scenarios.csv"
DEFAULT_OBSERVABLE_ORDER = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]


@dataclass(frozen=True)
class Scenario:
    name: str
    columns: tuple[str, ...]
    description: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run archive-based sequential ABC filtering for selected 3x3 parameters.")
    parser.add_argument("--target-parameters", type=Path, required=True)
    parser.add_argument("--bank-dir", type=Path, default=BROAD_BANK)
    parser.add_argument("--bed-output-dir", type=Path, default=DEFAULT_BED)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--day", type=int, default=81)
    parser.add_argument("--round-quantiles", nargs="+", type=float, default=[0.50, 0.25, 0.10, 0.05])
    parser.add_argument("--min-final-particles", type=int, default=100)
    return parser.parse_args()


def infer_observables(outputs: pd.DataFrame) -> list[str]:
    found = sorted({m.group(1) for col in outputs.columns if (m := re.match(r"^(.+)_day_(\d+)$", col))})
    ordered = [name for name in DEFAULT_OBSERVABLE_ORDER if name in found]
    ordered.extend(name for name in found if name not in ordered)
    return ordered


def load_targets(path: Path) -> list[str]:
    table = pd.read_csv(path)
    targets = table["parameter"].dropna().astype(str).tolist()
    if len(targets) != 9:
        raise ValueError(f"Expected 9 fixed representative parameters, got {len(targets)}")
    return targets


def load_bank(bank_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    parameters = pd.read_csv(bank_dir / "prior_parameter_samples.csv")
    outputs = pd.read_csv(bank_dir / "ode_output_features.csv")
    finite = parameters.notna().all(axis=1) & outputs.notna().all(axis=1)
    return parameters.loc[finite].reset_index(drop=True), outputs.loc[finite].reset_index(drop=True)


def day_columns(outputs: pd.DataFrame, day: int, biomarkers: list[str]) -> list[str]:
    cols = [f"{b}_day_{day}" for b in biomarkers if f"{b}_day_{day}" in outputs.columns]
    if not cols:
        raise ValueError(f"No observable columns found for day {day}")
    return cols


def column_label(columns: tuple[str, ...], max_items: int = 6) -> str:
    parts = []
    for column in columns:
        if "_day_" in column:
            biomarker, day = column.rsplit("_day_", 1)
            parts.append(f"{biomarker} d{day}")
        else:
            parts.append(column)
    shown = parts[:max_items]
    if len(parts) > max_items:
        shown.append(f"+{len(parts) - max_items} more")
    return ", ".join(shown)


def gsa_uncertainty_scenario(outputs: pd.DataFrame, observables: list[str]) -> Scenario | None:
    if not DEFAULT_GSA_UNCERTAINTY_SCENARIOS.exists():
        return None
    table = pd.read_csv(DEFAULT_GSA_UNCERTAINTY_SCENARIOS)
    required = {"biomarker", "day"}
    if not required.issubset(table.columns):
        return None
    table = table[table["biomarker"].isin(observables)].copy()
    if table.empty:
        return None
    counts = (
        table.assign(pair=table["biomarker"].astype(str) + "_day_" + table["day"].astype(int).astype(str))
        .groupby(["pair", "biomarker", "day"], as_index=False)
        .size()
        .sort_values(["size", "biomarker"], ascending=[False, True])
    )
    columns: list[str] = []
    for _, row in counts.iterrows():
        column = f"{row['biomarker']}_day_{int(row['day'])}"
        if column in outputs.columns and column not in columns:
            columns.append(column)
        if len(columns) >= 9:
            break
    if not columns:
        return None
    return Scenario(
        name="gsa_uncertainty_guided",
        columns=tuple(columns),
        description=f"GSA + uncertainty guided: {column_label(tuple(columns), max_items=9)}",
    )


def parse_biomarker_day_combination(combination: str, outputs: pd.DataFrame, observables: list[str]) -> tuple[str, ...]:
    columns: list[str] = []
    for item in str(combination).split("|"):
        if "_day_" not in item:
            continue
        biomarker, day_text = item.rsplit("_day_", 1)
        biomarker = biomarker.strip()
        try:
            day = int(day_text)
        except ValueError:
            continue
        column = f"{biomarker}_day_{day}"
        if biomarker in observables and column in outputs.columns and column not in columns:
            columns.append(column)
    return tuple(columns)


def parameter_specific_guided_scenarios(
    outputs: pd.DataFrame,
    targets: list[str],
    observables: list[str],
    max_columns: int = 5,
) -> dict[str, Scenario]:
    if DEFAULT_BED_GUIDED_SCENARIOS.exists():
        guided = pd.read_csv(DEFAULT_BED_GUIDED_SCENARIOS)
        required_guided = {"parameter", "scenario", "biomarker_day_combination"}
        if required_guided.issubset(guided.columns):
            scenarios: dict[str, Scenario] = {}
            guided = guided[guided["parameter"].isin(targets)].copy()
            for _, row in guided.iterrows():
                columns = parse_biomarker_day_combination(row["biomarker_day_combination"], outputs, observables)
                if columns:
                    parameter = str(row["parameter"])
                    scenarios[parameter] = Scenario(
                        name=str(row["scenario"]),
                        columns=columns,
                        description=f"{parameter}: {column_label(columns, max_items=max_columns)}",
                    )
            if scenarios:
                return scenarios

    if not DEFAULT_GSA_UNCERTAINTY_SCENARIOS.exists():
        return {}
    table = pd.read_csv(DEFAULT_GSA_UNCERTAINTY_SCENARIOS)
    required = {"parameter", "biomarker", "day"}
    if not required.issubset(table.columns):
        return {}
    table = table[table["parameter"].isin(targets) & table["biomarker"].isin(observables)].copy()
    if table.empty:
        return {}
    if "biomarker_rank" not in table.columns:
        table["biomarker_rank"] = 999
    if "window_rank" not in table.columns:
        table["window_rank"] = 999
    scenarios: dict[str, Scenario] = {}
    for parameter, group in table.groupby("parameter", sort=False):
        columns: list[str] = []
        ordered = group.sort_values(["biomarker_rank", "window_rank", "biomarker", "day"])
        for _, row in ordered.iterrows():
            column = f"{row['biomarker']}_day_{int(row['day'])}"
            if column in outputs.columns and column not in columns:
                columns.append(column)
            if len(columns) >= max_columns:
                break
        if columns:
            scenarios[parameter] = Scenario(
                name=f"gsa_uncertainty_guided_{parameter}",
                columns=tuple(columns),
                description=f"{parameter}: {column_label(tuple(columns), max_items=max_columns)}",
            )
    return scenarios


def scenario_columns_from_bed(bed_dir: Path, outputs: pd.DataFrame, default_day: int, observables: list[str]) -> list[Scenario]:
    scenarios: list[Scenario] = []
    guided = gsa_uncertainty_scenario(outputs, observables)
    if guided is not None:
        scenarios.append(guided)
    selection_path = bed_dir / "target_parameter_selection.csv"
    if selection_path.exists():
        selection = pd.read_csv(selection_path)
        species: list[str] = []
        best_days: list[int] = []
        low_days: list[int] = []
        for _, row in selection.iterrows():
            for biomarker in str(row.get("best_species", "")).split("|"):
                if biomarker and biomarker in observables and biomarker not in species:
                    species.append(biomarker)
            if pd.notna(row.get("best_day", np.nan)):
                best_days.append(int(row["best_day"]))
            if pd.notna(row.get("low_day", np.nan)):
                low_days.append(int(row["low_day"]))
        best_day = int(pd.Series(best_days).mode().iloc[0]) if best_days else default_day
        low_day = int(pd.Series(low_days).mode().iloc[0]) if low_days else default_day
        for idx, biomarker in enumerate(species[:4], start=1):
            scenarios.append(
                Scenario(
                    name=f"independent_{idx}_{biomarker}_day_{best_day}",
                    columns=tuple(day_columns(outputs, best_day, [biomarker])),
                    description=f"Independent observation scenario {idx}: {biomarker}",
                )
            )
        cumulative: list[str] = []
        for biomarker in species:
            cumulative.append(biomarker)
            scenarios.append(
                Scenario(
                    name=f"cumulative_best_{len(cumulative):02d}",
                    columns=tuple(day_columns(outputs, best_day, cumulative)),
                    description=f"Cumulative biomarkers: {' + '.join(cumulative)}",
                )
            )
        scenarios.append(Scenario("high_information_day_all_biomarkers", tuple(day_columns(outputs, best_day, observables)), "All observables at high-information day"))
        scenarios.append(Scenario("low_information_day_all_biomarkers", tuple(day_columns(outputs, low_day, observables)), "All observables at low-information day"))
    if not scenarios:
        scenarios.append(Scenario(f"all_observables_day_{default_day}", tuple(day_columns(outputs, default_day, observables)), "All observables at requested day"))
    return scenarios


def weighted_quantile(values: np.ndarray, weights: np.ndarray, qs: list[float]) -> np.ndarray:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return np.interp(qs, cdf, values)


def scenario_distance(outputs: pd.DataFrame, columns: tuple[str, ...]) -> np.ndarray:
    y = outputs[list(columns)].to_numpy(float)
    observation = y[0]
    scale = np.maximum(np.std(y, axis=0, ddof=1), 1e-12)
    return np.sqrt(np.mean(((y - observation[None, :]) / scale[None, :]) ** 2, axis=1))


def abc_rounds(distance: np.ndarray, round_quantiles: list[float], min_final: int) -> pd.DataFrame:
    rows = []
    n = len(distance)
    for i, q in enumerate(round_quantiles, start=1):
        eps = float(np.quantile(distance, q))
        accepted = int(np.sum(distance <= eps))
        rows.append({"round": i, "epsilon": eps, "quantile": q, "accepted": accepted, "acceptance_fraction": accepted / n})
    if rows and rows[-1]["accepted"] < min_final:
        sorted_distance = np.sort(distance)
        eps = float(sorted_distance[min(min_final - 1, n - 1)])
        accepted = int(np.sum(distance <= eps))
        rows.append({"round": len(rows) + 1, "epsilon": eps, "quantile": np.nan, "accepted": accepted, "acceptance_fraction": accepted / n, "note": "min_final_particles_enforced"})
    return pd.DataFrame(rows)


def summarize_particles(parameters: pd.DataFrame, targets: list[str], mask: np.ndarray, distance: np.ndarray) -> pd.DataFrame:
    weights = np.exp(-0.5 * (distance[mask] / max(np.median(distance[mask]), 1e-12)) ** 2)
    weights /= weights.sum()
    rows = []
    for parameter in targets:
        prior = parameters[parameter].to_numpy(float)
        post = parameters.loc[mask, parameter].to_numpy(float)
        p05, p50, p95 = np.quantile(prior, [0.05, 0.50, 0.95])
        q05, q50, q95 = weighted_quantile(post, weights, [0.05, 0.50, 0.95])
        rows.append(
            {
                "parameter": parameter,
                "prior_median": p50,
                "posterior_median": q50,
                "prior_90_width": p95 - p05,
                "posterior_90_width": q95 - q05,
                "interval_reduction_fraction": 1.0 - (q95 - q05) / max(p95 - p05, 1e-12),
            }
        )
    return pd.DataFrame(rows)


def kde_on_grid(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 5 or len(np.unique(values)) < 3:
        hist, edges = np.histogram(values, bins=min(20, max(3, len(values))), density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        return np.interp(grid, centers, hist, left=0.0, right=0.0)
    try:
        density = gaussian_kde(values)(grid)
    except Exception:
        hist, edges = np.histogram(values, bins=30, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        density = np.interp(grid, centers, hist, left=0.0, right=0.0)
    area = np.trapezoid(density, grid)
    return density / area if area > 0 else density


def plot_prior_posterior(
    parameters: pd.DataFrame,
    targets: list[str],
    mask: np.ndarray,
    path: Path,
    title: str,
    scenario: Scenario,
) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 10), constrained_layout=True)
    for ax, parameter in zip(axes.ravel(), targets):
        prior = parameters[parameter].to_numpy(float)
        posterior = parameters.loc[mask, parameter].to_numpy(float)
        lo = min(np.nanmin(prior), np.nanmin(posterior))
        hi = max(np.nanmax(prior), np.nanmax(posterior))
        span = max(hi - lo, np.finfo(float).eps)
        grid = np.linspace(lo - 0.18 * span, hi + 0.18 * span, 600)
        prior_density = kde_on_grid(prior, grid)
        posterior_density = kde_on_grid(posterior, grid)
        ax.fill_between(grid, prior_density, color="#9ecae1", alpha=0.22)
        ax.plot(grid, prior_density, color="#1f77b4", linewidth=1.8, label="broad +/-5% prior KDE")
        ax.fill_between(grid, posterior_density, color="#fdae6b", alpha=0.24)
        ax.plot(grid, posterior_density, color="#ff7f0e", linewidth=1.9, label="archive ABC posterior KDE")
        ax.set_xlim(grid[0], grid[-1])
        ax.set_title(parameter, fontsize=9)
        ax.set_xlabel("Parameter value")
        ax.set_ylabel("Density")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=14)
    fig.text(
        0.5,
        0.012,
        f"Observation scenario: {scenario.description}; columns: {column_label(scenario.columns, max_items=12)}",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_parameter_specific_guided(
    parameters: pd.DataFrame,
    targets: list[str],
    masks: dict[str, np.ndarray],
    scenarios: dict[str, Scenario],
    path: Path,
) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 10), constrained_layout=True)
    for ax, parameter in zip(axes.ravel(), targets):
        prior = parameters[parameter].to_numpy(float)
        posterior = parameters.loc[masks[parameter], parameter].to_numpy(float)
        lo = min(np.nanmin(prior), np.nanmin(posterior))
        hi = max(np.nanmax(prior), np.nanmax(posterior))
        span = max(hi - lo, np.finfo(float).eps)
        grid = np.linspace(lo - 0.18 * span, hi + 0.18 * span, 600)
        prior_density = kde_on_grid(prior, grid)
        posterior_density = kde_on_grid(posterior, grid)
        ax.fill_between(grid, prior_density, color="#9ecae1", alpha=0.22)
        ax.plot(grid, prior_density, color="#1f77b4", linewidth=1.8, label="broad +/-5% prior KDE")
        ax.fill_between(grid, posterior_density, color="#fdae6b", alpha=0.24)
        ax.plot(grid, posterior_density, color="#ff7f0e", linewidth=1.9, label="archive ABC posterior KDE")
        ax.set_xlim(grid[0], grid[-1])
        ax.set_title(f"{parameter}\n{column_label(scenarios[parameter].columns, max_items=5)}", fontsize=8)
        ax.set_xlabel("Parameter value")
        ax.set_ylabel("Density")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8)
    fig.suptitle("Parameter-specific GSA + uncertainty + MI guided archive ABC", fontsize=14)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    targets = load_targets(args.target_parameters)
    parameters, outputs = load_bank(args.bank_dir)
    observables = infer_observables(outputs)
    scenarios = scenario_columns_from_bed(args.bed_output_dir, outputs, args.day, observables)
    parameter_guided = parameter_specific_guided_scenarios(outputs, targets, observables)
    round_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    final_masks: dict[str, np.ndarray] = {}
    for scenario in scenarios:
        distance = scenario_distance(outputs, scenario.columns)
        rounds = abc_rounds(distance, args.round_quantiles, args.min_final_particles)
        rounds.insert(0, "scenario", scenario.name)
        rounds.insert(1, "description", scenario.description)
        round_frames.append(rounds)
        epsilon = float(rounds.iloc[-1]["epsilon"])
        mask = distance <= epsilon
        final_masks[scenario.name] = mask
        particles = parameters.loc[mask, targets].copy()
        particles["abc_distance"] = distance[mask]
        particles.to_csv(args.output_dir / f"abc_smc_particles_final_{scenario.name}.csv", index=False)
        summary = summarize_particles(parameters, targets, mask, distance)
        summary.insert(0, "scenario", scenario.name)
        summary.insert(1, "description", scenario.description)
        summary["final_epsilon"] = epsilon
        summary["final_particles"] = int(mask.sum())
        summary_frames.append(summary)
    parameter_guided_masks: dict[str, np.ndarray] = {}
    parameter_guided_rows: list[pd.DataFrame] = []
    for parameter, scenario in parameter_guided.items():
        distance = scenario_distance(outputs, scenario.columns)
        rounds = abc_rounds(distance, args.round_quantiles, args.min_final_particles)
        epsilon = float(rounds.iloc[-1]["epsilon"])
        mask = distance <= epsilon
        parameter_guided_masks[parameter] = mask
        summary = summarize_particles(parameters, [parameter], mask, distance)
        summary.insert(0, "scenario", scenario.name)
        summary.insert(1, "description", scenario.description)
        summary["final_epsilon"] = epsilon
        summary["final_particles"] = int(mask.sum())
        parameter_guided_rows.append(summary)

    rounds_all = pd.concat(round_frames, ignore_index=True)
    summaries = pd.concat(summary_frames, ignore_index=True)
    parameter_guided_summary = pd.concat(parameter_guided_rows, ignore_index=True) if parameter_guided_rows else pd.DataFrame()
    rounds_all.to_csv(args.output_dir / "abc_smc_round_summary.csv", index=False)
    rounds_all.to_csv(args.output_dir / "abc_smc_distance_summary.csv", index=False)
    summaries.to_csv(args.output_dir / "abc_smc_posterior_summary.csv", index=False)
    summaries.to_csv(args.output_dir / "abc_smc_posterior_narrowing_summary.csv", index=False)
    parameter_guided_summary.to_csv(args.output_dir / "abc_smc_parameter_specific_guided_posterior_summary.csv", index=False)
    scenario_audit = pd.DataFrame(
        [
            {"scenario": s.name, "description": s.description, "columns": "|".join(s.columns), "final_particles": int(final_masks[s.name].sum())}
            for s in scenarios
        ]
    )
    scenario_audit.to_csv(args.output_dir / "abc_smc_scenarios.csv", index=False)
    if parameter_guided:
        pd.DataFrame(
            [
                {
                    "parameter": parameter,
                    "scenario": scenario.name,
                    "description": scenario.description,
                    "columns": "|".join(scenario.columns),
                    "final_particles": int(parameter_guided_masks[parameter].sum()),
                }
                for parameter, scenario in parameter_guided.items()
            ]
        ).to_csv(args.output_dir / "abc_smc_parameter_specific_guided_scenarios.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
    for scenario, group in rounds_all.groupby("scenario", sort=False):
        ax.plot(group["round"], group["epsilon"], marker="o", label=scenario)
    ax.set_xlabel("ABC filtering round")
    ax.set_ylabel("epsilon")
    ax.set_title("Archive ABC threshold schedule by observation scenario")
    ax.legend(fontsize=6, ncol=2)
    fig.savefig(args.figure_dir / "abc_smc_epsilon_by_round.png", dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
    for scenario, group in rounds_all.groupby("scenario", sort=False):
        ax.plot(group["round"], group["acceptance_fraction"], marker="o", label=scenario)
    ax.set_xlabel("ABC filtering round")
    ax.set_ylabel("Archive acceptance fraction")
    ax.set_title("Archive ABC acceptance by scenario")
    ax.legend(fontsize=6, ncol=2)
    fig.savefig(args.figure_dir / "abc_smc_acceptance_by_round.png", dpi=220)
    plt.close(fig)
    scenario_by_name = {scenario.name: scenario for scenario in scenarios}
    main_scenario = "cumulative_best_09" if "cumulative_best_09" in final_masks else scenarios[-1].name
    guided_scenario = "gsa_uncertainty_guided" if "gsa_uncertainty_guided" in final_masks else main_scenario
    plot_prior_posterior(
        parameters,
        targets,
        final_masks[main_scenario],
        args.figure_dir / "abc_smc_prior_vs_posterior_3x3.png",
        f"Archive ABC prior vs posterior: {main_scenario}",
        scenario_by_name[main_scenario],
    )
    plot_prior_posterior(
        parameters,
        targets,
        final_masks[guided_scenario],
        args.figure_dir / "abc_smc_gsa_uncertainty_guided_posteriors_3x3.png",
        "Common GSA + uncertainty guided archive ABC",
        scenario_by_name[guided_scenario],
    )
    if set(targets).issubset(parameter_guided_masks):
        plot_parameter_specific_guided(
            parameters,
            targets,
            parameter_guided_masks,
            parameter_guided,
            args.figure_dir / "abc_smc_parameter_specific_gsa_uncertainty_guided_posteriors_3x3.png",
        )

    final_tables = ROOT / "results_final/tables"
    final_figures = ROOT / "results_final/figures"
    final_tables.mkdir(parents=True, exist_ok=True)
    final_figures.mkdir(parents=True, exist_ok=True)
    rounds_all.to_csv(final_tables / "abc_smc_round_summary.csv", index=False)
    summaries.to_csv(final_tables / "abc_smc_posterior_narrowing_summary.csv", index=False)
    parameter_guided_summary.to_csv(final_tables / "abc_smc_parameter_specific_guided_posterior_narrowing_summary.csv", index=False)
    scenario_audit.to_csv(final_tables / "abc_smc_observation_scenarios.csv", index=False)
    parameter_guided_scenarios_path = args.output_dir / "abc_smc_parameter_specific_guided_scenarios.csv"
    if parameter_guided_scenarios_path.exists():
        (final_tables / "abc_smc_parameter_specific_guided_scenarios.csv").write_bytes(parameter_guided_scenarios_path.read_bytes())
    for figure in [
        "abc_smc_epsilon_by_round.png",
        "abc_smc_acceptance_by_round.png",
        "abc_smc_prior_vs_posterior_3x3.png",
        "abc_smc_gsa_uncertainty_guided_posteriors_3x3.png",
        "abc_smc_parameter_specific_gsa_uncertainty_guided_posteriors_3x3.png",
    ]:
        src = args.figure_dir / figure
        if src.exists():
            (final_figures / figure).write_bytes(src.read_bytes())

    metadata = {
        "target_parameters": targets,
        "observable_biomarkers": observables,
        "n_archive_rows": int(len(parameters)),
        "n_scenarios": int(len(scenarios)),
        "n_parameter_specific_guided_scenarios": int(len(parameter_guided)),
        "round_quantiles": args.round_quantiles,
        "scientific_caveat": "Archive ABC here is likelihood-free filtering of an existing broad-prior ODE archive; no surrogate-predicted candidates are treated as truth.",
    }
    (args.output_dir / "abc_smc_metadata.json").write_text(json.dumps(metadata, indent=2))
    print(f"archive ABC scenarios evaluated: {len(scenarios)}")
    print(f"archive rows used: {len(parameters)}")
    print(f"main scenario: {main_scenario}; final particles: {int(final_masks[main_scenario].sum())}")
    print(f"guided scenario: {guided_scenario}; {scenario_by_name[guided_scenario].description}")
    print(f"parameter-specific guided scenarios: {len(parameter_guided)}")
    print(f"outputs saved: {args.output_dir}")
    print("caveat: ABC particles are ODE archive rows only; no surrogate-predicted samples are used as scientific truth.")


if __name__ == "__main__":
    main()
