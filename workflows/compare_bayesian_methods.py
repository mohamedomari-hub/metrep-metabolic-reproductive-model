"""Compare targeted Bayesian update methods."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare posterior reweighting, archive posterior, and archive ABC summaries.")
    parser.add_argument("--bed-dir", type=Path, default=ROOT / "local_outputs/bed_targeted")
    parser.add_argument("--mcmc-dir", type=Path, default=ROOT / "local_outputs/reduced_archive_posterior")
    parser.add_argument("--abc-dir", type=Path, default=ROOT / "local_outputs/archive_abc")
    parser.add_argument("--table-dir", type=Path, default=ROOT / "results_final/tables")
    parser.add_argument("--figure-dir", type=Path, default=ROOT / "results_final/figures")
    return parser.parse_args()


def load_method(path: Path, method: str, source: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    table = pd.read_csv(path)
    if "interval_reduction_fraction" not in table.columns and "interval_width_reduction_fraction" in table.columns:
        table["interval_reduction_fraction"] = table["interval_width_reduction_fraction"]
    if "posterior_90_width" not in table.columns and "posterior_90_width" not in table:
        return pd.DataFrame()
    keep = [c for c in ["parameter", "prior_90_width", "posterior_90_width", "interval_reduction_fraction", "posterior_median", "scenario", "biomarkers"] if c in table.columns]
    out = table[keep].copy()
    out["method"] = method
    out["source_file"] = source
    out["computational_cost"] = {
        "posterior reweighting": "low; saved ODE bank weights",
        "reduced ODE-archive posterior": "low/moderate; saved ODE bank likelihood weights",
        "archive-based sequential ABC filtering": "moderate/high; ODE archive likelihood-free filtering",
    }.get(method, "not recorded")
    out["interpretation"] = {
        "posterior reweighting": "fast BED design comparison",
        "reduced ODE-archive posterior": "likelihood-based selected-parameter inference on precomputed ODE rows",
        "archive-based sequential ABC filtering": "likelihood-free selected-parameter inference over precomputed ODE rows",
    }.get(method, "")
    return out


def main() -> None:
    args = parse_args()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    frames = [
        load_method(args.bed_dir / "posterior_narrowing_cumulative_biomarkers.csv", "posterior reweighting", "posterior_narrowing_cumulative_biomarkers.csv"),
        load_method(args.mcmc_dir / "reduced_posterior_narrowing_summary.csv", "reduced ODE-archive posterior", "reduced_posterior_narrowing_summary.csv"),
        load_method(args.abc_dir / "abc_smc_posterior_narrowing_summary.csv", "archive-based sequential ABC filtering", "abc_smc_posterior_narrowing_summary.csv"),
    ]
    comparison = pd.concat([f for f in frames if not f.empty], ignore_index=True) if any(not f.empty for f in frames) else pd.DataFrame()
    comparison.to_csv(args.table_dir / "bayesian_method_comparison.csv", index=False)
    if not comparison.empty:
        summary = comparison.groupby("method", as_index=False)["interval_reduction_fraction"].mean()
        fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
        ax.bar(summary["method"], summary["interval_reduction_fraction"])
        ax.set_ylabel("Mean interval reduction fraction")
        ax.set_title("Bayesian update method comparison")
        ax.tick_params(axis="x", rotation=20)
        fig.savefig(args.figure_dir / "bayesian_method_comparison_summary.png", dpi=220)
        plt.close(fig)
    print("Bayesian method comparison saved")


if __name__ == "__main__":
    main()
