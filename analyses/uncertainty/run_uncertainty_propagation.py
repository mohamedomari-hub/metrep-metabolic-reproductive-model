"""Uncertainty propagation from the existing admissible 0.5% MetRep bank."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BANK = PROJECT_ROOT / "analyses/bayesian_experimental_design/surrogate_bed/input_tables_large"
BIOMARKERS = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose"]
REPRODUCTIVE = ["FSH", "PGF", "P4", "E2", "INH"]
METABOLIC = ["IGF1", "Insulin", "Glucose"]
README_BIOMARKERS = ["Glucose", "Insulin", "P4", "E2"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Propagate uncertainty using stored admissible simulations only.")
    parser.add_argument("--bank-dir", type=Path, default=DEFAULT_BANK)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "outputs")
    return parser.parse_args()


def stored_days(outputs: pd.DataFrame, biomarker: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(biomarker)}_day_(\d+)$")
    return sorted(int(match.group(1)) for column in outputs.columns if (match := pattern.match(column)))


def trajectory_quantiles(outputs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for biomarker in BIOMARKERS:
        for day in stored_days(outputs, biomarker):
            values = outputs[f"{biomarker}_day_{day}"].to_numpy(dtype=float)
            q05, median, q95 = np.quantile(values, [0.05, 0.5, 0.95])
            rows.append(
                {
                    "biomarker": biomarker,
                    "day": day,
                    "q05": q05,
                    "median": median,
                    "q95": q95,
                    "interval_width": q95 - q05,
                    "samples": len(values),
                }
            )
    return pd.DataFrame(rows)


def plot_group(
    quantiles: pd.DataFrame,
    nominal: pd.DataFrame,
    biomarkers: list[str],
    section_title: str,
    path: Path,
    ncols: int = 3,
) -> None:
    nrows = int(np.ceil(len(biomarkers) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4.2 * nrows), constrained_layout=True, squeeze=False)
    for ax, biomarker in zip(axes.ravel(), biomarkers):
        subset = quantiles[quantiles["biomarker"] == biomarker].sort_values("day")
        day = subset["day"].to_numpy(dtype=float)
        q05 = subset["q05"].to_numpy(dtype=float)
        median = subset["median"].to_numpy(dtype=float)
        q95 = subset["q95"].to_numpy(dtype=float)
        nominal_values = nominal[[f"{biomarker}_day_{int(value)}" for value in day]].iloc[0].to_numpy(dtype=float)
        ax.fill_between(day, q05, q95, color="tab:blue", alpha=0.22, label="5th-95th percentile")
        ax.plot(day, median, color="tab:blue", linewidth=2, label="Median")
        ax.plot(day, nominal_values, color="black", linestyle="--", linewidth=1.5, label="Nominal")
        ax.set_title(biomarker)
        ax.set_xlabel("Day")
        ax.set_ylabel("Model output")
        ax.grid(alpha=0.25)
    for ax in axes.ravel()[len(biomarkers):]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=8)
    fig.suptitle(
        "Uncertainty propagation under biologically admissible parameter variability\n"
        f"{section_title} | n = {int(quantiles['samples'].iloc[0]):,} admissible simulations | "
        "narrow +/-0.5% parameter ensemble",
        fontsize=15,
    )
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    parameters = pd.read_csv(args.bank_dir / "prior_parameter_samples.csv")
    outputs = pd.read_csv(args.bank_dir / "ode_output_features.csv")
    nominal = pd.read_csv(args.bank_dir / "nominal_output.csv")
    admissibility = pd.read_csv(args.bank_dir / "admissibility.csv")
    if not (len(parameters) == len(outputs) == len(admissibility)):
        raise ValueError("Stored bank tables have inconsistent row counts.")
    keep = admissibility["admissible"].astype(bool) & outputs.notna().all(axis=1)
    accepted = outputs.loc[keep].reset_index(drop=True)

    quantiles = trajectory_quantiles(accepted)
    quantiles.to_csv(args.output_dir / "trajectory_quantiles.csv", index=False)
    summary = (
        quantiles.groupby("biomarker", as_index=False)
        .agg(
            first_day=("day", "min"),
            last_day=("day", "max"),
            maximum_interval_width=("interval_width", "max"),
            mean_interval_width=("interval_width", "mean"),
            samples=("samples", "first"),
        )
    )
    summary.to_csv(args.output_dir / "uncertainty_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "stored_simulations": len(outputs),
                "admissible_simulations": len(accepted),
                "rejected_simulations": len(outputs) - len(accepted),
                "prior_relative_half_range": 0.005,
                "ode_simulations_run": 0,
            }
        ]
    ).to_csv(args.output_dir / "accepted_simulation_summary.csv", index=False)

    plot_group(quantiles, nominal, REPRODUCTIVE, "Reproductive biomarkers", args.output_dir / "uncertainty_reproductive.png")
    plot_group(quantiles, nominal, METABOLIC, "Metabolic biomarkers", args.output_dir / "uncertainty_metabolic.png")
    plot_group(quantiles, nominal, BIOMARKERS, "All observable biomarkers", args.output_dir / "uncertainty_all_biomarkers.png")
    plot_group(
        quantiles,
        nominal,
        README_BIOMARKERS,
        "Representative observable biomarkers",
        args.output_dir / "uncertainty_readme_summary.png",
        ncols=2,
    )

    note = (
        "Local sensitivity used +1% one-at-a-time perturbations around the nominal model and AUC endpoints. "
        "Global sensitivity screening uses simulation-bank ensembles and AUC endpoints. The admissible-bank "
        "Spearman/PRCC analysis summarizes parameter-biomarker associations across biologically plausible "
        "simulations. Uncertainty propagation uses the biologically admissible 0.5% simulation bank to generate "
        "stable trajectory bands. Perturbation scales differ because each method asks a different question. "
        "BED is handled separately.\n"
    )
    (args.output_dir / "README.md").write_text(note)
    summary_note = (
        "# Uncertainty Propagation\n\n"
        "Uncertainty propagation was estimated from the biologically admissible Monte Carlo ensemble "
        "generated around the calibrated parameter regime. The shaded region shows the 5th-95th "
        "percentile range, the blue line shows the ensemble median, and the black dashed line shows "
        "the nominal trajectory. Because the ensemble uses a narrow +/-0.5% perturbation range and "
        "biological admissibility filtering, the bands should be interpreted as local robustness "
        "around the calibrated model, not as full population variability.\n"
    )
    (args.output_dir / "uncertainty_summary.md").write_text(summary_note)
    print(f"Saved uncertainty outputs to {args.output_dir}")
    print(f"Admissible rows: {len(accepted)}; ODE simulations run: 0")


if __name__ == "__main__":
    main()
