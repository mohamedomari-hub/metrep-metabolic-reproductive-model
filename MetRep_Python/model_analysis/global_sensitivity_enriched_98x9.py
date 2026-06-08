"""Final 98 x 9 admissible-ensemble PRCC/Spearman sensitivity outputs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _global_sensitivity_core import (
    REPRESENTATIVE_CLASSES,
    admissible_screening,
    auc_endpoints,
    infer_observable_biomarkers,
)


FORMATTED_ENRICHED = ROOT / "local_data/phd_bed_bank_5pct_50k_glucagon_smc_enriched"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute final enriched 98x9 PRCC/Spearman sensitivity.")
    parser.add_argument("--enriched-bank-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    return parser.parse_args()


def resolve_bank(path: Path) -> Path:
    if (path / "prior_parameter_samples.csv").exists():
        return path
    if FORMATTED_ENRICHED.exists():
        return FORMATTED_ENRICHED
    raise FileNotFoundError(f"{path} is not a formatted bank and {FORMATTED_ENRICHED} was not found.")


def heatmap(table: pd.DataFrame, value: str, path: Path, title: str) -> None:
    matrix = table.pivot(index="parameter", columns="biomarker", values=value)
    order = matrix.abs().max(axis=1).sort_values().index
    matrix = matrix.loc[order]
    fig, ax = plt.subplots(figsize=(11, max(10, 0.16 * len(matrix))), constrained_layout=True)
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), matrix.columns)
    ax.set_yticks(range(len(matrix.index)), matrix.index, fontsize=5)
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label=value)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def representative_heatmap(prcc: pd.DataFrame, path: Path) -> None:
    matrix = prcc[prcc["parameter"].isin(REPRESENTATIVE_CLASSES)].pivot(index="parameter", columns="biomarker", values="prcc")
    matrix = matrix.reindex(list(REPRESENTATIVE_CLASSES))
    fig, ax = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), matrix.columns)
    labels = [f"{p}\n{REPRESENTATIVE_CLASSES[p]}" for p in matrix.index]
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_title("Representative profile-likelihood parameters: PRCC with observable AUC")
    fig.colorbar(image, ax=ax, label="PRCC")
    fig.savefig(path, dpi=240)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    bank = resolve_bank(args.enriched_bank_dir)
    parameters = pd.read_csv(bank / "prior_parameter_samples.csv")
    outputs = pd.read_csv(bank / "ode_output_features.csv")
    admissibility = pd.read_csv(bank / "admissibility.csv")
    keep = admissibility["admissible"].astype(bool) & parameters.notna().all(axis=1) & outputs.notna().all(axis=1)
    parameters = parameters.loc[keep].reset_index(drop=True)
    outputs = outputs.loc[keep].reset_index(drop=True)
    biomarkers = infer_observable_biomarkers(outputs)
    auc = auc_endpoints(outputs, biomarkers)
    spearman, prcc = admissible_screening(parameters, auc, biomarkers)

    prcc.to_csv(args.output_dir / "global_sensitivity_98x9_prcc.csv", index=False)
    spearman.to_csv(args.output_dir / "global_sensitivity_98x9_spearman.csv", index=False)
    prcc.to_csv(args.table_dir / "global_sensitivity_98x9_prcc.csv", index=False)
    spearman.to_csv(args.table_dir / "global_sensitivity_98x9_spearman.csv", index=False)

    top_links = prcc.sort_values("abs_prcc", ascending=False).head(100)
    top_links.to_csv(args.table_dir / "global_sensitivity_top_parameter_biomarker_links.csv", index=False)
    targets = list(REPRESENTATIVE_CLASSES)
    top_targets = (
        prcc[prcc["parameter"].isin(targets)]
        .sort_values(["parameter", "abs_prcc"], ascending=[True, False])
        .groupby("parameter")
        .head(5)
        .reset_index(drop=True)
    )
    top_targets.to_csv(args.table_dir / "global_sensitivity_top_biomarkers_for_profile_targets.csv", index=False)

    heatmap(prcc, "prcc", args.figure_dir / "global_sensitivity_98x9_prcc_heatmap.png", "Admissible ensemble PRCC: 98 parameters x 9 observables")
    heatmap(spearman, "spearman_rho", args.figure_dir / "global_sensitivity_98x9_spearman_heatmap.png", "Admissible ensemble Spearman: 98 parameters x 9 observables")
    representative_heatmap(prcc, args.figure_dir / "global_sensitivity_representative_parameters_heatmap.png")

    print(f"enriched bank count: {len(parameters)}")
    print("global sensitivity 98 x 9 outputs saved")


if __name__ == "__main__":
    main()
