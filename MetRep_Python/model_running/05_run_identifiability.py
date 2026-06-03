"""Run SVD-based local structural identifiability diagnostics.

This script follows the approach from the older
``metrep_structid_svd_merged_robust_updated.py`` workflow, adapted to the
refactored MetRep package:

    1. central-difference stacked sensitivity matrix
    2. SVD rank/nullspace analysis
    3. sensitivity ranking metrics
    4. nullspace participation
    5. compensation edge extraction
    6. holistic estimate/fix recommendation

Run from project root:
    python model_running/05_run_identifiability.py

By default this analyzes all 98 MATLAB-derived core parameters for 90 days.
Dexa PK/PD constants are not part of this parameter list; Dexa is an optional
perturbation runner and does not interfere with SVD/sensitivity/profile
likelihood workflows.
"""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PROJECT_ROOT.parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))
sys.path.insert(0, str(WORKSPACE_ROOT))
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt

from model_definition.analysis import structural_identifiability_svd
from model_definition.parameters import PARAMETERS
from model_definition.plotting import (
    parameter_scenario_classes,
    plot_compensation_network,
    plot_parameter_scenario_map,
    plot_sensitivity_ranked_by_class,
)
from model_definition.scenarios import constant_non_lactating


FOCUSED_PARAMETERS = [
    "feed_direct_blood_fraction",
    "insulin_secretion_max",
    "insulin_clearance",
    "glucagon_secretion_max",
    "glucagon_clearance",
    "hepatic_glucose_release",
    "blood_to_liver_max",
    "liver_to_storage_max",
    "storage_to_liver_max",
    "blood_usage_max",
    "insulin_igf_scale",
    "igf_clearance",
    "cl_to_p4_scale",
    "p4_clearance",
    "follicle_to_e2_scale",
    "e2_clearance",
]

DEFAULT_OUTPUTS = [
    "Glucose",
    "Insulin",
    "IGF1",
    "P4",
    "E2",
    "Follicle",
    "CL",
]


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_singular_values_plot(path: Path, singular_values: np.ndarray, threshold: float, rank: int) -> None:
    """Save log singular-value plot."""

    apply_plot_style()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    ax.semilogy(np.arange(1, len(singular_values) + 1), singular_values, marker="o")
    ax.axhline(threshold, linestyle="--", linewidth=1, color="black")
    ax.set_xlabel("Singular value index")
    ax.set_ylabel("Singular value (log scale)")
    ax.set_title(f"SVD of stacked sensitivity matrix\nrank={rank}/{len(singular_values)}")
    ax.grid(alpha=0.25, linewidth=0.8)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_barplot(path: Path, names: list[str], values: np.ndarray, title: str, xlabel: str) -> None:
    """Save horizontal bar plot for all selected parameters."""

    apply_plot_style()
    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values, dtype=float)
    normalized = values / values.max() if values.size and values.max() > 0 else values
    order = np.argsort(normalized)[::-1]

    fig_height = max(7.0, 0.28 * len(names))
    fig, ax = plt.subplots(figsize=(12, fig_height), constrained_layout=True)
    ax.barh([names[i] for i in order][::-1], normalized[order][::-1])
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25, linewidth=0.8)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_nullspace_heatmap(path: Path, nullspace: np.ndarray, parameter_names: list[str], max_labeled_params: int = 101) -> None:
    """Save normalized nullspace basis heatmap."""

    if nullspace.size == 0:
        return

    apply_plot_style()
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix = nullspace.copy()
    for row_index in range(matrix.shape[0]):
        max_abs = np.max(np.abs(matrix[row_index, :]))
        if max_abs > 0:
            matrix[row_index, :] = matrix[row_index, :] / max_abs

    fig, ax = plt.subplots(
        figsize=(max(10, 0.12 * matrix.shape[1]), max(4.5, 0.35 * matrix.shape[0])),
        constrained_layout=True,
    )
    image = ax.imshow(matrix, aspect="auto", interpolation="nearest")
    fig.colorbar(image, ax=ax, label="Coeff (normalized per null vector)")
    ax.set_xlabel("Parameter index")
    ax.set_ylabel("Nullspace basis vector index")
    ax.set_title("Nullspace basis heatmap (normalized rows)")
    if matrix.shape[1] <= max_labeled_params:
        step = max(1, matrix.shape[1] // 20)
        ticks = list(range(0, matrix.shape[1], step))
        ax.set_xticks(ticks)
        ax.set_xticklabels([parameter_names[index] for index in ticks], rotation=90, fontsize=8)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_holistic_scatter(path: Path, holistic: pd.DataFrame, thresholds: dict[str, float] | None = None) -> None:
    """Save presentation-friendly 3-class sensitivity/nullspace decision map."""

    apply_plot_style()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)

    color_map = {
        "Estimate": "#1f77b4",
        "Fix (anchor)": "#d62728",
        "Fix (irrelevant)": "#7f7f7f",
    }
    class_column = "recommendation_3class" if "recommendation_3class" in holistic.columns else "recommendation"
    legend_order = ["Estimate", "Fix (anchor)", "Fix (irrelevant)"]

    for recommendation in legend_order:
        group = holistic[holistic[class_column] == recommendation]
        if group.empty:
            continue
        ax.scatter(
            group["nullspace_0to1"],
            group["sensitivity_0to1"],
            label=recommendation,
            s=35,
            alpha=0.86,
            color=color_map.get(recommendation, "tab:purple"),
        )

    if thresholds:
        ax.axhline(thresholds.get("sens_hi", 0.0), linestyle="--", linewidth=1, color="0.25")
        ax.axhline(thresholds.get("sens_lo", 0.0), linestyle="--", linewidth=1, color="0.45")
        ax.axvline(thresholds.get("null_hi", 0.0), linestyle="--", linewidth=1, color="0.25")
        ax.axvline(thresholds.get("null_lo", 0.0), linestyle="--", linewidth=1, color="0.45")

    ax.set_xlabel("Nullspace involvement (0..1) -> higher = more compensatory")
    ax.set_ylabel("Sensitivity (0..1) -> higher = more output impact")
    ax.set_title("Holistic decision map: Estimate vs Fix (3-class)")
    ax.grid(alpha=0.25, linewidth=0.8)
    ax.legend(loc="best", framealpha=0.9)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_holistic_bar_allparams(path: Path, holistic: pd.DataFrame) -> None:
    """Save all parameters ranked by sensitivity and colored by 3-class recommendation."""

    apply_plot_style()
    path.parent.mkdir(parents=True, exist_ok=True)
    class_column = "recommendation_3class" if "recommendation_3class" in holistic.columns else "recommendation"
    ordered = holistic.sort_values("sensitivity_0to1", ascending=False)
    color_map = {
        "Estimate": "#1f77b4",
        "Fix (anchor)": "#d62728",
        "Fix (irrelevant)": "#7f7f7f",
    }
    names = ordered["parameter"].tolist()
    values = ordered["sensitivity_0to1"].to_numpy(dtype=float)
    classes = ordered[class_column].tolist()

    fig, ax = plt.subplots(figsize=(max(14, 0.25 * len(names)), 8), constrained_layout=True)
    ax.bar(range(len(names)), values, color=[color_map.get(name, "tab:purple") for name in classes], width=0.9)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_ylabel("Sensitivity (0..1)")
    ax.set_title("ALL parameters ranked by sensitivity (3-class FIX/ESTIMATE recommendation)")
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)

    from matplotlib.patches import Patch

    handles = [
        Patch(facecolor=color_map[label], label=label)
        for label in ["Estimate", "Fix (anchor)", "Fix (irrelevant)"]
        if label in set(classes)
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.9)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_nullspace_text_report(
    path: Path,
    nullspace: np.ndarray,
    parameter_names: list[str],
    threshold_relative: float = 0.35,
    top_k: int = 10,
) -> None:
    """Write dominant coefficients and compensation patterns for each null direction."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        if nullspace.size == 0:
            handle.write("Nullspace empty: no local structural non-identifiability directions detected.\n")
            return

        handle.write(
            f"Nullspace basis vectors: {nullspace.shape[0]} directions, "
            f"D={nullspace.shape[1]} parameters\n"
        )
        handle.write(
            f"Reporting coeffs: abs(v) >= {threshold_relative}*max(abs(v)) "
            f"OR top-{top_k}\n\n"
        )

        for row_index, row in enumerate(nullspace):
            values = row.astype(float)
            values_abs = np.abs(values)
            max_abs = values_abs.max()
            if max_abs <= 0:
                continue

            keep = np.where(values_abs >= threshold_relative * max_abs)[0]
            if keep.size > top_k:
                keep = np.argsort(values_abs)[::-1][:top_k]
            keep = keep[np.argsort(values_abs[keep])[::-1]]

            handle.write(f"=== Null direction {row_index + 1}/{nullspace.shape[0]} ===\n")
            for index in keep:
                handle.write(f"  {parameter_names[int(index)]:>32s}: {values[int(index)]: .6e}\n")

            positive = [parameter_names[int(index)] for index in keep if values[int(index)] > 0]
            negative = [parameter_names[int(index)] for index in keep if values[int(index)] < 0]
            if positive and negative:
                handle.write(f"  Compensation pattern: increase {positive} while decreasing {negative}\n")
            handle.write("\n")


def save_holistic_summary(path: Path, holistic: pd.DataFrame, thresholds: dict[str, float]) -> None:
    """Write class counts and quantile thresholds used for the holistic recommendation."""

    path.parent.mkdir(parents=True, exist_ok=True)
    class_column = "recommendation_3class" if "recommendation_3class" in holistic.columns else "recommendation"
    counts = Counter(holistic[class_column].tolist())
    lines = ["Holistic recommendation summary (Sensitivity + Nullspace involvement)", ""]
    for label, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        lines.append(f"{label}: {count}")
    lines.append("")
    lines.append("Thresholds (on 0..1 normalized scores):")
    for key in ["sens_lo", "sens_hi", "null_lo", "null_hi", "sens_q_lo", "sens_q_hi", "null_q_lo", "null_q_hi"]:
        if key in thresholds:
            lines.append(f"  {key}: {thresholds[key]}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_horizons(raw: str) -> list[float]:
    horizons: list[float] = []
    for token in (raw or "").split(","):
        token = token.strip()
        if not token:
            continue
        value = float(token)
        if value not in horizons:
            horizons.append(value)
    return horizons


def save_robust_consensus(holistic_paths: list[Path], horizons: list[float], output_prefix: Path) -> None:
    """Build robust consensus and horizon-stability figures from multiple holistic tables."""

    if len(holistic_paths) < 2:
        print("[robust] Not enough successful runs to build consensus.")
        return

    frames = []
    for path, horizon in zip(holistic_paths, horizons):
        frame = pd.read_csv(path)
        frame["horizon_days"] = float(horizon)
        if "recommendation_3class" not in frame.columns:
            frame["recommendation_3class"] = frame["recommendation"]
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    parameter_column = "param" if "param" in combined.columns else "parameter"
    pivot = combined.pivot_table(
        index=parameter_column,
        columns="horizon_days",
        values="recommendation_3class",
        aggfunc="first",
    )

    preference = {"Estimate": 0, "Fix (anchor)": 1, "Fix (irrelevant)": 2}

    def mode_with_ties(values: list[str]) -> str:
        counts: dict[str, int] = {}
        for value in values:
            counts[value] = counts.get(value, 0) + 1
        return sorted(counts.items(), key=lambda item: (-item[1], preference.get(item[0], 9), item[0]))[0][0]

    rows = []
    for parameter, row in pivot.iterrows():
        recommendations = [value for value in row.values.tolist() if isinstance(value, str) and value.strip()]
        if recommendations:
            consensus = mode_with_ties(recommendations)
            stability = float(sum(value == consensus for value in recommendations)) / float(len(recommendations))
        else:
            consensus = ""
            stability = 0.0
        rows.append((parameter, consensus, stability, ";".join(recommendations)))

    consensus = pd.DataFrame(
        rows,
        columns=["param", "consensus_recommendation", "stability_frac", "per_horizon_recs"],
    )
    mean_scores = combined.groupby(parameter_column)[
        ["sensitivity_0to1", "nullspace_0to1", "identifiable_0to1"]
    ].mean().reset_index()
    mean_scores = mean_scores.rename(columns={parameter_column: "param"})
    consensus = consensus.merge(mean_scores, on="param", how="left")
    consensus = consensus.sort_values(["stability_frac", "sensitivity_0to1"], ascending=[False, False])

    figure_prefix = PROJECT_ROOT / "results" / "figures" / output_prefix.name
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_prefix.with_name(output_prefix.name + "_ROBUST_consensus.csv")
    consensus.to_csv(csv_path, index=False)

    class_to_int = {"Estimate": 0, "Fix (anchor)": 1, "Fix (irrelevant)": 2}
    matrix = pivot.map(lambda value: class_to_int.get(value, np.nan))
    matrix = matrix.loc[consensus["param"].values]

    from matplotlib.colors import BoundaryNorm, ListedColormap

    cmap = ListedColormap(["#1f77b4", "#d62728", "#7f7f7f"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    heatmap_path = figure_prefix.with_name(figure_prefix.name + "_ROBUST_class_heatmap.png")
    fig, ax = plt.subplots(figsize=(max(8.0, 1.0 * len(horizons) + 4.0), max(4.5, 0.16 * len(matrix) + 2.0)))
    image = ax.imshow(matrix.values, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    step = 1 if len(matrix.index) <= 80 else max(2, len(matrix.index) // 80)
    ticks = list(range(0, len(matrix.index), step))
    ax.set_yticks(ticks)
    ax.set_yticklabels([matrix.index[index] for index in ticks], fontsize=6)
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels([str(int(float(column))) for column in matrix.columns], fontsize=10)
    ax.set_xlabel("Horizon (days)")
    ax.set_title("Recommendation across horizons (categorical)")
    cbar = fig.colorbar(image, ax=ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["Estimate", "Fix (anchor)", "Fix (irrelevant)"])
    fig.tight_layout()
    fig.savefig(heatmap_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    classbar_path = figure_prefix.with_name(figure_prefix.name + "_ROBUST_class_fractions.png")
    counts = []
    for column in matrix.columns:
        values = matrix[column].values
        counts.append((np.sum(values == 0), np.sum(values == 1), np.sum(values == 2)))
    counts_array = np.asarray(counts, dtype=float)
    totals = counts_array.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    fractions = counts_array / totals

    fig, ax = plt.subplots(figsize=(max(8.0, 1.2 * len(matrix.columns) + 4.0), 4.8), constrained_layout=True)
    x = np.arange(len(matrix.columns))
    ax.bar(x, fractions[:, 0], label="Estimate", color="#1f77b4")
    ax.bar(x, fractions[:, 1], bottom=fractions[:, 0], label="Fix (anchor)", color="#d62728")
    ax.bar(
        x,
        fractions[:, 2],
        bottom=fractions[:, 0] + fractions[:, 1],
        label="Fix (irrelevant)",
        color="#7f7f7f",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(float(column))) for column in matrix.columns])
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Fraction of parameters")
    ax.set_xlabel("Horizon (days)")
    ax.set_title("Recommendation composition by horizon")
    ax.legend(loc="upper right")
    fig.savefig(classbar_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    stability_path = figure_prefix.with_name(figure_prefix.name + "_ROBUST_stability.png")
    fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)
    ax.plot(consensus["stability_frac"].values)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Parameters (sorted by stability then sensitivity)")
    ax.set_ylabel("Stability fraction")
    ax.set_title("Stability of recommendation across horizons")
    fig.savefig(stability_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print("Saved robust consensus outputs:")
    print(f"  {csv_path}")
    print(f"  {heatmap_path}")
    print(f"  {classbar_path}")
    print(f"  {stability_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.description = (
        "Run SVD identifiability for the 22-state non-Dexa MetRep core. "
        "Dexa PK/PD constants are intentionally excluded."
    )
    parser.add_argument("--all-params", action="store_true", help="Use all 98 MATLAB-derived parameters. This is now the default.")
    parser.add_argument("--focused-params", action="store_true", help="Use the historical focused 16-parameter subset.")
    parser.add_argument("--params", nargs="+", default=None, help="Explicit parameter names to analyze.")
    parser.add_argument("--outputs", nargs="+", default=DEFAULT_OUTPUTS, help="Observable state names.")
    parser.add_argument("--days", type=float, default=90.0, help="Simulation horizon in days.")
    parser.add_argument("--dt", type=float, default=2.0, help="Output spacing in days.")
    parser.add_argument("--rel-step", type=float, default=1e-3, help="Central-difference relative step.")
    parser.add_argument("--tol-value", type=float, default=1e-8, help="Relative SVD rank tolerance.")
    parser.add_argument("--rank-metric", choices=["colnorm", "rel_colnorm", "rel2_colnorm"], default="rel2_colnorm")
    parser.add_argument("--null-thr-rel", type=float, default=0.35, help="Nullspace coefficient threshold relative to max(abs(v)).")
    parser.add_argument("--null-topk", type=int, default=10, help="Top-k coefficients per null vector to report if threshold keeps too many.")
    parser.add_argument("--edges-top", type=int, default=60, help="Maximum compensation edges to show in the network plot.")
    parser.add_argument("--network-all-nodes", action="store_true", help="Show all analyzed parameters in the compensation network.")
    parser.add_argument(
        "--network-layout",
        choices=["filtered", "all-zones", "diagnostic"],
        default="filtered",
        help="Compensation network layout: filtered strongest-edge network, all analyzed parameters in zones, or diagnostic multi-figure output.",
    )
    parser.add_argument("--robust", action="store_true", help="Run robust identifiability consensus across multiple horizons.")
    parser.add_argument("--horizons", default="42,100,224", help="Comma-separated horizons for --robust.")
    parser.add_argument("--prefix", default="structid_svd_refactored", help="Output filename prefix.")
    return parser.parse_args()


def select_parameter_names(args: argparse.Namespace) -> list[str]:
    if args.params:
        return args.params
    if args.focused_params:
        return FOCUSED_PARAMETERS
    return [parameter.name for parameter in PARAMETERS]


def run_identifiability_once(
    args: argparse.Namespace,
    parameter_names: list[str],
    prefix: str,
    days: float,
) -> Path:
    scenario = constant_non_lactating(days=days, step=args.dt)

    result = structural_identifiability_svd(
        scenario,
        parameter_names=parameter_names,
        outputs=args.outputs,
        relative_step=args.rel_step,
        tolerance_value=args.tol_value,
        ranking_metric=args.rank_metric,
        null_threshold_relative=args.null_thr_rel,
        null_top_k=args.null_topk,
    )

    table_dir = PROJECT_ROOT / "results" / "tables"
    figure_dir = PROJECT_ROOT / "results" / "figures"
    simulation_dir = PROJECT_ROOT / "results" / "simulations"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    simulation_dir.mkdir(parents=True, exist_ok=True)

    singular_values = result["singular_values"]
    rank = result["rank"]
    nullity = result["nullity"]
    threshold = result["threshold"]
    nullspace = result["nullspace"]
    holistic = result["holistic_table"]
    thresholds = result.get("holistic_thresholds", {})

    np.savez_compressed(
        simulation_dir / f"{prefix}.npz",
        params=np.array(parameter_names, dtype=object),
        outputs=np.array(args.outputs, dtype=object),
        t=scenario.t_eval,
        nominal_outputs=result["nominal_outputs"],
        sensitivity_matrix=result["sensitivity_matrix"],
        singular_values=singular_values,
        vt=result["vt"],
        rank=np.array([rank]),
        nullity=np.array([nullity]),
        threshold=np.array([threshold]),
        rank_metric=np.array([args.rank_metric], dtype=object),
        ranking_scores=result["ranking_scores"],
        nullspace_score=result["nullspace_score"],
        identifiable_score=result["identifiable_score"],
        holistic_category=holistic["recommendation_3class"].to_numpy(dtype=object),
        holistic_sens01=holistic.set_index("parameter").loc[parameter_names, "sensitivity_0to1"].to_numpy(dtype=float),
        holistic_null01=holistic.set_index("parameter").loc[parameter_names, "nullspace_0to1"].to_numpy(dtype=float),
        holistic_id01=holistic.set_index("parameter").loc[parameter_names, "identifiable_0to1"].to_numpy(dtype=float),
    )

    pd.DataFrame({"singular_value": singular_values}).to_csv(
        table_dir / f"{prefix}_singular_values.csv", index=False
    )
    pd.DataFrame(
        {
            "parameter": parameter_names,
            "ranking_score": result["ranking_scores"],
            "nullspace_score": result["nullspace_score"],
            "identifiable_score": result["identifiable_score"],
        }
    ).sort_values("ranking_score", ascending=False).to_csv(
        table_dir / f"{prefix}_ranking_and_participation.csv", index=False
    )
    holistic_path = table_dir / f"{prefix}_holistic_table.csv"
    holistic.to_csv(holistic_path, index=False)
    result["compensation_edges"].to_csv(table_dir / f"{prefix}_compensation_edges.csv", index=False)

    save_singular_values_plot(
        figure_dir / f"{prefix}_singular_values.png",
        singular_values,
        threshold,
        rank,
    )
    save_barplot(
        figure_dir / f"{prefix}_ranking_{args.rank_metric}.png",
        parameter_names,
        result["ranking_scores"],
        f"Local sensitivity ranking ({args.rank_metric})",
        "normalized score",
    )
    save_barplot(
        figure_dir / f"{prefix}_ranking_ALL_{args.rank_metric}.png",
        parameter_names,
        result["ranking_scores"],
        f"Local sensitivity-based ranking (ALL parameters)\nmetric={args.rank_metric}",
        f"Relative score (normalized) - {args.rank_metric}",
    )
    save_barplot(
        figure_dir / f"{prefix}_nullspace_participation.png",
        parameter_names,
        result["nullspace_score"],
        "Nullspace participation",
        "normalized nullspace score",
    )
    save_barplot(
        figure_dir / f"{prefix}_null_participation_ALL.png",
        parameter_names,
        result["nullspace_score"],
        "Nullspace participation (ALL parameters)\n(higher = more involved in compensations)",
        "Relative nullspace participation (normalized)",
    )
    save_nullspace_heatmap(
        figure_dir / f"{prefix}_nullspace_heatmap.png",
        nullspace,
        parameter_names,
    )
    null_report_path = figure_dir / f"{prefix}_nullspace_report.txt"
    write_nullspace_text_report(
        null_report_path,
        nullspace,
        parameter_names,
        threshold_relative=args.null_thr_rel,
        top_k=args.null_topk,
    )
    save_holistic_scatter(
        figure_dir / f"{prefix}_holistic_scatter_ns_vs_sens.png",
        holistic,
        thresholds=thresholds,
    )
    save_holistic_bar_allparams(
        figure_dir / f"{prefix}_holistic_bar_allparams.png",
        holistic,
    )
    holistic_summary_path = figure_dir / f"{prefix}_holistic_summary.txt"
    save_holistic_summary(
        holistic_summary_path,
        holistic,
        thresholds,
    )
    core_network_path = figure_dir / (
        f"{prefix}_compensation_network_core.png"
        if args.network_layout == "diagnostic"
        else f"{prefix}_compensation_network_filtered.png"
    )
    network_path = plot_compensation_network(
        result["compensation_edges"],
        nullspace,
        parameter_names,
        core_network_path,
        holistic_table=holistic,
        threshold_relative=args.null_thr_rel,
        topk=args.null_topk,
        max_edges=args.edges_top,
        include_all_nodes=False,
        layout_mode="filtered",
        include_singletons=args.network_layout != "diagnostic",
        title="Core compensation network: strongest parameter trade-offs"
        if args.network_layout == "diagnostic"
        else None,
    )
    all_network_path = None
    scenario_map_path = None
    sensitivity_rank_path = None
    scenario_classes_path = None
    if args.network_layout == "diagnostic":
        scenario_map_path = plot_parameter_scenario_map(
            holistic,
            figure_dir / f"{prefix}_parameter_scenario_map.png",
            thresholds=thresholds,
        )
        sensitivity_rank_path = plot_sensitivity_ranked_by_class(
            holistic,
            figure_dir / f"{prefix}_sensitivity_ranked_by_class.png",
        )
        scenario_classes_path = table_dir / f"{prefix}_parameter_scenario_classes.csv"
        parameter_scenario_classes(holistic, thresholds=thresholds).to_csv(
            scenario_classes_path,
            index=False,
        )
    if args.network_layout == "all-zones" or args.network_all_nodes:
        all_network_path = plot_compensation_network(
            result["compensation_edges"],
            nullspace,
            parameter_names,
            figure_dir / f"{prefix}_compensation_network_all_parameters_supplement.png",
            holistic_table=holistic,
            threshold_relative=args.null_thr_rel,
            topk=args.null_topk,
            max_edges=args.edges_top,
            include_all_nodes=True,
            layout_mode="all-zones",
        )

    print("SVD identifiability summary")
    print(f"parameters: {len(parameter_names)}")
    print(f"outputs: {', '.join(args.outputs)}")
    print(f"duration: {days:g} days, dt={args.dt:g} days")
    print(f"rank: {rank}/{len(parameter_names)}")
    print(f"nullity: {nullity}")
    print(f"threshold: {threshold:.3e}")
    print(f"holistic table: {holistic_path}")
    print(f"nullspace report: {null_report_path}")
    print(f"holistic summary: {holistic_summary_path}")
    print()
    if network_path is not None:
        print(f"Compensation network filtered: {network_path}")
    else:
        print("Compensation network filtered: skipped (empty graph or networkx unavailable)")
    if all_network_path is not None:
        print(f"Compensation network all-parameter supplement: {all_network_path}")
    if scenario_map_path is not None:
        print(f"Parameter scenario map: {scenario_map_path}")
    if sensitivity_rank_path is not None:
        print(f"Sensitivity ranked by class: {sensitivity_rank_path}")
    if scenario_classes_path is not None:
        print(f"Parameter scenario classes: {scenario_classes_path}")
    print()
    print("Top holistic recommendations")
    print(holistic.head(12).to_string(index=False))
    return holistic_path


def main() -> None:
    args = parse_args()
    parameter_names = select_parameter_names(args)

    if args.robust:
        horizons = parse_horizons(args.horizons) or [42.0, 100.0, 224.0]
        holistic_paths: list[Path] = []
        completed_horizons: list[float] = []
        for horizon in horizons:
            horizon_prefix = f"{args.prefix}_H{int(float(horizon))}"
            run_args = copy.deepcopy(args)
            run_args.robust = False
            print(f"[robust] Running horizon={horizon:g} days -> prefix={horizon_prefix}")
            try:
                holistic_paths.append(
                    run_identifiability_once(run_args, parameter_names, horizon_prefix, float(horizon))
                )
                completed_horizons.append(float(horizon))
            except Exception as exc:
                print(f"[robust][warn] Horizon {horizon:g} failed: {exc}")
        output_prefix = PROJECT_ROOT / "results" / "tables" / args.prefix
        save_robust_consensus(holistic_paths, completed_horizons, output_prefix)
        return

    run_identifiability_once(args, parameter_names, args.prefix, args.days)


if __name__ == "__main__":
    main()
