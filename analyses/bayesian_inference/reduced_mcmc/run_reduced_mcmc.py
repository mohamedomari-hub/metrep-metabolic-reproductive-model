"""Reduced 3x3-parameter Bayesian inference from the saved ODE archive.

This script does not run new ODE simulations and does not use a surrogate. It
conditions the fixed profile-likelihood 3x3 representative parameter set on
observation scenarios by likelihood-weighting precomputed broad +/-5% ODE rows,
then resampling posterior rows for plotting and summaries.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "MetRep_Python"))

from model_definition.parameters import default_parameters


BROAD_BANK = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon"
DEFAULT_BED = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/run_outputs_targeted"
DEFAULT_OBSERVABLE_ORDER = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]


@dataclass(frozen=True)
class Scenario:
    name: str
    columns: tuple[str, ...]
    description: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run reduced selected-parameter posterior inference from the ODE archive.")
    parser.add_argument("--target-parameters", type=Path, required=True)
    parser.add_argument("--bank-dir", type=Path, default=BROAD_BANK)
    parser.add_argument("--bed-output-dir", type=Path, default=DEFAULT_BED)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-steps", type=int, default=20000, help="Number of posterior rows to resample from the weighted ODE archive.")
    parser.add_argument("--burn-in", type=int, default=5000, help="Accepted for CLI compatibility; not used by archive reweighting.")
    parser.add_argument("--thin", type=int, default=5, help="Accepted for CLI compatibility; not used by archive reweighting.")
    parser.add_argument("--day", type=int, default=81)
    parser.add_argument("--validation-median-nrmse-threshold", type=float, default=0.12, help="Deprecated; retained for old commands.")
    parser.add_argument("--validation-worst-nrmse-threshold", type=float, default=0.30, help="Deprecated; retained for old commands.")
    return parser.parse_args()


def infer_observables(outputs: pd.DataFrame) -> list[str]:
    found = sorted({m.group(1) for col in outputs.columns if (m := re.match(r"^(.+)_day_(\d+)$", col))})
    ordered = [name for name in DEFAULT_OBSERVABLE_ORDER if name in found]
    ordered.extend(name for name in found if name not in ordered)
    return ordered


def load_targets(path: Path) -> list[str]:
    table = pd.read_csv(path)
    if "parameter" not in table.columns:
        raise ValueError(f"{path} must contain a 'parameter' column")
    targets = table["parameter"].dropna().astype(str).tolist()
    if len(targets) != 9:
        raise ValueError(f"Expected fixed 3x3 representative target set of 9 parameters, got {len(targets)}")
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


def column_label(columns: tuple[str, ...], max_items: int = 8) -> str:
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


def scenario_columns_from_bed(bed_dir: Path, outputs: pd.DataFrame, default_day: int, observables: list[str]) -> list[Scenario]:
    scenarios: list[Scenario] = []
    selection_path = bed_dir / "target_parameter_selection.csv"
    if selection_path.exists():
        selection = pd.read_csv(selection_path)
        all_species: list[str] = []
        best_days: list[int] = []
        low_days: list[int] = []
        for _, row in selection.iterrows():
            for b in str(row.get("best_species", "")).split("|"):
                if b and b in observables and b not in all_species:
                    all_species.append(b)
            if pd.notna(row.get("best_day", np.nan)):
                best_days.append(int(row["best_day"]))
            if pd.notna(row.get("low_day", np.nan)):
                low_days.append(int(row["low_day"]))
        best_day = int(pd.Series(best_days).mode().iloc[0]) if best_days else default_day
        low_day = int(pd.Series(low_days).mode().iloc[0]) if low_days else default_day
        if all_species:
            for idx, biomarker in enumerate(all_species[:4], start=1):
                scenarios.append(
                    Scenario(
                        name=f"independent_{idx}_{biomarker}_day_{best_day}",
                        columns=tuple(day_columns(outputs, best_day, [biomarker])),
                        description=f"Independent observation scenario {idx}: {biomarker} at BED high-information day",
                    )
                )
            cumulative = []
            for biomarker in all_species:
                cumulative.append(biomarker)
                scenarios.append(
                    Scenario(
                        name=f"cumulative_best_{len(cumulative):02d}",
                        columns=tuple(day_columns(outputs, best_day, cumulative)),
                        description=f"Cumulative BED biomarkers: {' + '.join(cumulative)}",
                    )
                )
            scenarios.append(
                Scenario(
                    name="high_information_day_all_biomarkers",
                    columns=tuple(day_columns(outputs, best_day, observables)),
                    description="All observables at BED high-information day",
                )
            )
            scenarios.append(
                Scenario(
                    name="low_information_day_all_biomarkers",
                    columns=tuple(day_columns(outputs, low_day, observables)),
                    description="All observables at BED low-information day",
                )
            )
    if not scenarios:
        scenarios.append(
            Scenario(
                name=f"all_observables_day_{default_day}",
                columns=tuple(day_columns(outputs, default_day, observables)),
                description="All observables at requested day",
            )
        )
    return scenarios


def safe_observation(outputs: pd.DataFrame, columns: tuple[str, ...]) -> np.ndarray:
    return outputs.loc[0, list(columns)].to_numpy(float)


def weighted_archive_posterior(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    targets: list[str],
    columns: tuple[str, ...],
    observation: np.ndarray,
    rng: np.random.Generator,
    n_samples: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    sigma = np.maximum(0.10 * np.abs(observation), 1e-8)
    residual = (outputs.loc[:, list(columns)].to_numpy(float) - observation) / sigma
    log_likelihood = -0.5 * np.sum(residual * residual, axis=1)
    shifted = log_likelihood - np.nanmax(log_likelihood)
    weights = np.exp(shifted)
    weights = weights / np.sum(weights)
    indices = rng.choice(np.arange(len(parameters)), size=n_samples, replace=True, p=weights)
    samples = parameters.iloc[indices][targets].reset_index(drop=True)
    ess = 1.0 / np.sum(weights * weights)
    diagnostics = {
        "posterior_rows": int(len(samples)),
        "effective_sample_size": float(ess),
        "max_weight": float(np.max(weights)),
        "n_archive_rows": int(len(parameters)),
    }
    return samples, diagnostics


def summarize_posterior(samples: pd.DataFrame, lower: np.ndarray, upper: np.ndarray) -> pd.DataFrame:
    rows = []
    for i, parameter in enumerate(samples.columns):
        q05, q50, q95 = samples[parameter].quantile([0.05, 0.50, 0.95])
        prior_width = upper[i] - lower[i]
        post_width = q95 - q05
        rows.append(
            {
                "parameter": parameter,
                "posterior_q05": q05,
                "posterior_median": q50,
                "posterior_q95": q95,
                "prior_90_width": 0.90 * prior_width,
                "posterior_90_width": post_width,
                "interval_reduction_fraction": 1.0 - post_width / max(0.90 * prior_width, 1e-12),
            }
        )
    return pd.DataFrame(rows)


def kde_on_grid(values: np.ndarray, grid: np.ndarray, bw_adjust: float = 1.0) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 5 or len(np.unique(values)) < 3:
        hist, edges = np.histogram(values, bins=min(20, max(3, len(values))), density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        return np.interp(grid, centers, hist, left=0.0, right=0.0)
    try:
        kde = gaussian_kde(values)
        kde.set_bandwidth(kde.factor * bw_adjust)
        density = kde(grid)
    except Exception:
        hist, edges = np.histogram(values, bins=30, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        density = np.interp(grid, centers, hist, left=0.0, right=0.0)
    area = np.trapezoid(density, grid)
    return density / area if area > 0 else density


def plot_prior_vs_posterior(
    targets: list[str],
    samples: pd.DataFrame,
    lower: np.ndarray,
    upper: np.ndarray,
    path: Path,
    title: str,
    scenario_note: str = "",
) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 10), constrained_layout=True)
    rng = np.random.default_rng(123)
    for ax, parameter, lo, hi in zip(axes.ravel(), targets, lower, upper):
        prior = rng.uniform(lo, hi, 12000)
        posterior = samples[parameter].to_numpy(float)
        xlo = min(np.nanmin(prior), np.nanmin(posterior))
        xhi = max(np.nanmax(prior), np.nanmax(posterior))
        span = max(xhi - xlo, np.finfo(float).eps)
        grid = np.linspace(xlo - 0.18 * span, xhi + 0.18 * span, 600)
        prior_density = kde_on_grid(prior, grid, bw_adjust=1.4)
        posterior_density = kde_on_grid(posterior, grid, bw_adjust=3.0)
        ax.fill_between(grid, prior_density, color="#9ecae1", alpha=0.22)
        ax.plot(grid, prior_density, color="#1f77b4", linewidth=1.8, label="broad +/-5% prior KDE")
        ax.fill_between(grid, posterior_density, color="#fdae6b", alpha=0.24)
        ax.plot(grid, posterior_density, color="#ff7f0e", linewidth=1.9, label="ODE-archive posterior KDE")
        ax.set_xlim(grid[0], grid[-1])
        ax.set_title(parameter, fontsize=9)
        ax.set_xlabel("Parameter value")
        ax.set_ylabel("Density")
        ax.grid(alpha=0.20)
    axes.ravel()[0].legend(fontsize=8)
    fig.suptitle(title, fontsize=14)
    if scenario_note:
        fig.text(0.5, 0.012, scenario_note, ha="center", va="bottom", fontsize=9)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    targets = load_targets(args.target_parameters)
    parameters, outputs = load_bank(args.bank_dir)
    missing = [p for p in targets if p not in parameters.columns]
    if missing:
        raise ValueError(f"Missing target parameters in bank: {missing}")
    observables = infer_observables(outputs)
    scenarios = scenario_columns_from_bed(args.bed_output_dir, outputs, args.day, observables)
    nominal = default_parameters()
    lower = np.array([0.95 * nominal[p] for p in targets], dtype=float)
    upper = np.array([1.05 * nominal[p] for p in targets], dtype=float)

    scenario_rows: list[dict[str, object]] = []
    all_summaries: list[pd.DataFrame] = []
    main_samples: pd.DataFrame | None = None
    main_scenario_note = ""
    for scenario in scenarios:
        observation = safe_observation(outputs, scenario.columns)
        samples, chain_diag = weighted_archive_posterior(
            parameters=parameters,
            outputs=outputs,
            targets=targets,
            columns=scenario.columns,
            observation=observation,
            rng=rng,
            n_samples=args.n_steps,
        )
        samples.to_csv(args.output_dir / f"reduced_posterior_samples_{scenario.name}.csv", index=False)
        summary = summarize_posterior(samples, lower, upper)
        summary.insert(0, "scenario", scenario.name)
        summary.insert(1, "description", scenario.description)
        summary["effective_sample_size"] = chain_diag["effective_sample_size"]
        summary["max_weight"] = chain_diag["max_weight"]
        all_summaries.append(summary)
        scenario_rows.append(
            {
                "scenario": scenario.name,
                "description": scenario.description,
                "columns": "|".join(scenario.columns),
                "effective_sample_size": chain_diag["effective_sample_size"],
                "max_weight": chain_diag["max_weight"],
                "posterior_rows": chain_diag["posterior_rows"],
                "inference_status": "ode_archive_likelihood_weighted",
            }
        )
        if main_samples is None or scenario.name.startswith("cumulative_best_09") or scenario.name.startswith("high_information"):
            main_samples = samples
            main_scenario_note = f"Observation scenario: {scenario.description}; columns: {column_label(scenario.columns, max_items=12)}"

    scenario_table = pd.DataFrame(scenario_rows)
    scenario_table.to_csv(args.output_dir / "reduced_posterior_scenario_summary.csv", index=False)
    if all_summaries:
        posterior_summary = pd.concat(all_summaries, ignore_index=True)
        posterior_summary.to_csv(args.output_dir / "reduced_posterior_summary.csv", index=False)
        posterior_summary.to_csv(args.output_dir / "reduced_posterior_narrowing_summary.csv", index=False)
    else:
        posterior_summary = pd.DataFrame()
        posterior_summary.to_csv(args.output_dir / "reduced_posterior_summary.csv", index=False)
        posterior_summary.to_csv(args.output_dir / "reduced_posterior_narrowing_summary.csv", index=False)

    diagnostics = {
        "target_parameters": targets,
        "observable_biomarkers": observables,
        "n_bank_rows_used": int(len(parameters)),
        "n_scenarios": int(len(scenarios)),
        "inference_mode": "ode_archive_likelihood_weighting_no_surrogate",
        "scientific_caveat": "This reduced posterior uses precomputed ODE archive rows only. It does not use surrogate predictions and does not run new ODE simulations.",
    }
    (args.output_dir / "reduced_posterior_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    pd.DataFrame([diagnostics | {"target_parameters": "|".join(targets), "observable_biomarkers": "|".join(observables)}]).to_csv(args.output_dir / "mcmc_diagnostics.csv", index=False)

    if main_samples is not None:
        plot_prior_vs_posterior(targets, main_samples, lower, upper, args.figure_dir / "mcmc_prior_vs_posterior_3x3.png", "Reduced MCMC prior vs posterior", main_scenario_note)
        plot_prior_vs_posterior(targets, main_samples, lower, upper, args.figure_dir / "mcmc_gsa_uncertainty_guided_posteriors_3x3.png", "Reduced ODE-archive posterior", main_scenario_note)
        corr = np.corrcoef(main_samples[targets].to_numpy(float), rowvar=False)
        fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
        im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(targets)), targets, rotation=90, fontsize=7)
        ax.set_yticks(range(len(targets)), targets, fontsize=7)
        fig.colorbar(im, ax=ax, label="posterior correlation")
        fig.savefig(args.figure_dir / "mcmc_pairplot_or_correlation.png", dpi=220)
        plt.close(fig)

    final_tables = ROOT / "results_final/tables"
    final_figures = ROOT / "results_final/figures"
    final_tables.mkdir(parents=True, exist_ok=True)
    final_figures.mkdir(parents=True, exist_ok=True)
    scenario_table.to_csv(final_tables / "mcmc_scenario_validation_summary.csv", index=False)
    scenario_table.to_csv(final_tables / "reduced_posterior_scenario_summary.csv", index=False)
    posterior_summary.to_csv(final_tables / "mcmc_posterior_narrowing_summary.csv", index=False)
    posterior_summary.to_csv(final_tables / "reduced_posterior_narrowing_summary.csv", index=False)
    for figure in ["mcmc_prior_vs_posterior_3x3.png", "mcmc_gsa_uncertainty_guided_posteriors_3x3.png", "mcmc_pairplot_or_correlation.png"]:
        src = args.figure_dir / figure
        if src.exists():
            (final_figures / figure).write_bytes(src.read_bytes())

    print(f"target parameters: {targets}")
    print(f"scenarios evaluated: {len(scenarios)}")
    print("surrogate used: no")
    print("inference mode: ODE-archive likelihood weighting")
    print(f"outputs saved: {args.output_dir}")
    print("caveat: no new ODEs are run; posterior support is limited to the saved broad +/-5% ODE archive.")


if __name__ == "__main__":
    main()
