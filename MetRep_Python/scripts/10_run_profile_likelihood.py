"""Synthetic-data profile likelihood with biological admissibility penalties.

Default outputs are the measurable MetRep variables:
FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon.

The default mode is plot-first and storage-light: it saves figures plus one
small summary CSV. Detailed profile tables are optional via ``--save-tables``.
"""

from __future__ import annotations

import argparse
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

from metrep.admissibility import AdmissibilityThresholds
from metrep.initial_conditions import STATE_INDEX
from metrep.parameters import PARAMETERS, default_parameters
from metrep.profile_likelihood import (
    SyntheticData,
    classify_profile,
    make_synthetic_data,
    profile_one_parameter,
)
from metrep.scenarios import constant_non_lactating
from metrep.simulate import run_simulation

try:
    from scipy.interpolate import PchipInterpolator
except Exception:  # pragma: no cover - SciPy is a project dependency.
    PchipInterpolator = None


MEASURABLE_OUTPUTS = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "figure.titlesize": 17,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-identifiability", type=Path, help="Holistic table from script 05.")
    parser.add_argument("--params", nargs="+", help="Explicit parameters to profile.")
    parser.add_argument("--outputs", nargs="+", default=MEASURABLE_OUTPUTS)
    parser.add_argument("--days", type=float, default=90.0)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--sample-every", type=float, default=2.0)
    parser.add_argument("--noise-relative", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--profile-class",
        choices=["estimate", "anchor", "irrelevant", "all"],
        default="estimate",
        help="Parameter class selected from the local SVD holistic table.",
    )
    parser.add_argument("--max-profile-params", type=int, default=8)
    parser.add_argument("--max-nuisance-params", type=int, default=12)
    parser.add_argument("--include-anchor-limit", type=int, default=4)
    parser.add_argument("--grid-points", type=int, default=21)
    parser.add_argument("--grid-low", type=float, default=0.5)
    parser.add_argument("--grid-high", type=float, default=1.5)
    parser.add_argument(
        "--trajectory-multipliers",
        nargs="+",
        type=float,
        default=[0.5, 0.75, 1.0, 1.25, 1.5],
        help="Representative multipliers shown in perturbed trajectory plots.",
    )
    parser.add_argument(
        "--trajectory-parameter",
        default="first",
        help="Which parameter gets the expensive trajectory overlay: first, none, all, or a parameter name.",
    )
    parser.add_argument("--nuisance-low", type=float, default=0.5)
    parser.add_argument("--nuisance-high", type=float, default=1.5)
    parser.add_argument("--maxiter", type=int, default=60)
    parser.add_argument(
        "--profile-ymax",
        default="auto",
        help="Y-axis maximum for profile plots. Use a number, or 'auto' for per-plot scaling.",
    )
    parser.add_argument(
        "--profile-ymax-override",
        nargs="+",
        default=[],
        help="Per-parameter y-axis overrides, e.g. insulin_glucose_threshold=2.5 inhibin_clearance=5.",
    )
    parser.add_argument(
        "--profile-scale",
        choices=["linear", "log1p"],
        default="log1p",
        help="Y-axis scale for profile plots. log1p keeps large-loss curves visible without clipping.",
    )
    parser.add_argument("--profile-cutoff", type=float, default=1.92, help="Delta-loss cutoff shown in profile plots.")
    parser.add_argument("--admissible-only", action="store_true", help="Use only biologically admissible points in profile views/classification.")
    parser.add_argument("--rho-min", type=float, default=0.75)
    parser.add_argument("--gamma-max", type=float, default=0.30)
    parser.add_argument("--kappa-max", type=float, default=0.30)
    parser.add_argument("--max-shift-days", type=float, default=5.0)
    parser.add_argument("--penalty-weight", type=float, default=100.0)
    parser.add_argument("--prefix", default="profile_likelihood")
    parser.add_argument("--save-tables", action="store_true", help="Save detailed per-grid profile tables.")
    parser.add_argument("--save-synthetic-csv", action="store_true", help="Save synthetic observations as CSV.")
    return parser.parse_args()


def _column(table: pd.DataFrame, candidates: list[str]) -> str:
    for candidate in candidates:
        if candidate in table.columns:
            return candidate
    raise KeyError(f"None of these columns were found: {candidates}")


def parse_ymax_overrides(items: list[str]) -> dict[str, float]:
    overrides = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Expected PARAM=VALUE for --profile-ymax-override, got {item!r}")
        name, value = item.split("=", 1)
        overrides[name.strip()] = float(value)
    return overrides


def select_parameters(args: argparse.Namespace) -> tuple[list[str], list[str]]:
    """Select profile and nuisance parameters from explicit args or SVD table."""

    all_parameter_names = {parameter.name for parameter in PARAMETERS}
    if args.params:
        profile_parameters = [name for name in args.params if name in all_parameter_names]
        return profile_parameters[: args.max_profile_params], profile_parameters[: args.max_nuisance_params]

    if not args.from_identifiability:
        raise SystemExit("Use --params or --from-identifiability to select profile parameters.")

    table = pd.read_csv(args.from_identifiability)
    parameter_column = _column(table, ["parameter", "param"])
    class_column = "recommendation_3class" if "recommendation_3class" in table.columns else "recommendation"
    sensitivity_column = "sensitivity_0to1" if "sensitivity_0to1" in table.columns else "sensitivity_raw"
    ordered = table.sort_values(sensitivity_column, ascending=False)

    estimate = ordered[ordered[class_column].isin(["Estimate", "ESTIMATE_candidate", "ESTIMATE_or_constrain"])]
    anchor = ordered[
        ordered[class_column].isin(["Fix (anchor)", "FIX_anchor_compensation", "COMPENSATORY_high_impact"])
    ]
    irrelevant = ordered[
        ordered[class_column].isin(["Fix (irrelevant)", "FIX_low_sensitivity"])
    ]
    if args.profile_class == "estimate":
        selected = estimate
    elif args.profile_class == "anchor":
        selected = anchor
    elif args.profile_class == "irrelevant":
        selected = irrelevant
    else:
        selected = ordered

    profile_parameters = selected[parameter_column].tolist()
    if args.max_profile_params > 0:
        profile_parameters = profile_parameters[: args.max_profile_params]
    anchor_parameters = anchor[parameter_column].tolist()[: args.include_anchor_limit]
    nuisance = (profile_parameters + anchor_parameters)[: args.max_nuisance_params]

    profile_parameters = [name for name in profile_parameters if name in all_parameter_names]
    nuisance = [name for name in nuisance if name in all_parameter_names]
    return profile_parameters, nuisance


def should_plot_trajectory(parameter: str, profile_parameters: list[str], request: str) -> bool:
    request = request.strip()
    if request == "all":
        return True
    if request == "none":
        return False
    if request == "first":
        return bool(profile_parameters) and parameter == profile_parameters[0]
    return parameter == request


def plot_synthetic_data(data: SyntheticData, output_path: Path) -> Path:
    apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ncols = 3
    nrows = int(np.ceil(len(data.outputs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.7 * ncols, 3.0 * nrows), sharex=True, constrained_layout=True)
    axes_flat = np.asarray(axes).ravel()
    for ax, output, values in zip(axes_flat, data.outputs, data.values.T):
        ax.plot(data.times, values, marker="o", markersize=3.2, linewidth=1.4, color="tab:blue")
        ax.set_title(output)
        ax.set_ylabel(output)
        ax.grid(alpha=0.25, linewidth=0.8)
    for ax in axes_flat[len(data.outputs) :]:
        ax.axis("off")
    for ax in axes_flat[-ncols:]:
        if ax.has_data():
            ax.set_xlabel("Time (days)")
    fig.suptitle(f"Synthetic measurable outputs (relative noise = {data.noise_relative:g})", y=1.02)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def profile_view(profile: pd.DataFrame, admissible_only: bool) -> pd.DataFrame:
    """Return the profile points used for presentation views."""

    view = profile[profile["admissible"]].copy() if admissible_only else profile.copy()
    if view.empty:
        return view
    loss_column = "fit_loss" if admissible_only else "total_loss"
    view["view_delta_loss"] = view[loss_column] - float(view[loss_column].min())
    return view


def resolve_profile_ymax(
    values: np.ndarray,
    cutoff: float,
    requested: str | float,
) -> float:
    """Choose a readable profile y-axis limit."""

    if isinstance(requested, (int, float)):
        return float(requested)

    requested_text = str(requested).strip().lower()
    if requested_text not in {"auto", "adaptive"}:
        return float(requested)

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return max(5.0, 3.0 * cutoff)

    near = values[values <= max(10.0 * cutoff, cutoff + 1e-12)]
    if near.size >= 2:
        high = float(np.percentile(near, 95))
    else:
        finite_high = float(np.nanmin([np.nanmax(values), 10.0 * cutoff]))
        high = finite_high if np.isfinite(finite_high) else 3.0 * cutoff

    return max(1.35 * cutoff, 1.15 * high)


def transform_profile_y(values: np.ndarray, scale: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if scale == "log1p":
        return np.log1p(np.maximum(values, 0.0))
    return values


def plot_profile(
    profile: pd.DataFrame,
    output_path: Path,
    y_max: str | float = "auto",
    cutoff: float = 1.92,
    admissible_only: bool = False,
    profile_scale: str = "log1p",
) -> Path:
    apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    parameter = str(profile["profile_parameter"].iloc[0])
    profile = profile.sort_values("fixed_multiplier")
    view = profile_view(profile, admissible_only).sort_values("fixed_multiplier")
    fig, ax = plt.subplots(figsize=(8, 5.2), constrained_layout=True)
    if view.empty:
        ax.text(0.5, 0.5, "No biologically admissible profile points", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"Profile likelihood: {parameter}")
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return output_path

    x = view["fixed_multiplier"].to_numpy(dtype=float)
    y = view["view_delta_loss"].to_numpy(dtype=float)
    y_plot = transform_profile_y(y, profile_scale)
    cutoff_plot = float(transform_profile_y(np.array([cutoff]), profile_scale)[0])
    resolved_ymax = (
        float(y_max)
        if str(y_max).strip().lower() not in {"auto", "adaptive"}
        else max(1.15 * float(np.nanmax(y_plot)), 1.35 * cutoff_plot)
    )
    if PchipInterpolator is not None and len(view) >= 3:
        x_smooth = np.linspace(float(x.min()), float(x.max()), 240)
        y_smooth = PchipInterpolator(x, y_plot)(x_smooth)
        ax.plot(x_smooth, y_smooth, color="0.25", linewidth=2.0, label="profile curve")
    elif len(view) >= 2:
        ax.plot(x, y_plot, color="0.25", linewidth=1.8, label="profile curve")
    else:
        ax.plot([], [], color="0.25", linewidth=1.8, label="profile curve")
    ax.scatter(
        view["fixed_multiplier"],
        y_plot,
        color="tab:blue",
        label="admissible profile points" if admissible_only else "admissible",
        zorder=3,
    )
    if not admissible_only and np.any(~profile["admissible"].to_numpy(dtype=bool)):
        inadmissible_y = transform_profile_y(profile.loc[~profile["admissible"], "delta_loss"].to_numpy(dtype=float), profile_scale)
        ax.scatter(
            profile.loc[~profile["admissible"], "fixed_multiplier"],
            np.minimum(inadmissible_y, resolved_ymax),
            color="tab:red",
            label="penalized/not admissible",
            zorder=3,
            marker="x",
        )
    ax.axhline(cutoff_plot, linestyle="--", color="0.35", linewidth=1.0, label="approx. 95% cutoff")
    ax.axvline(1.0, linestyle=":", color="0.35", linewidth=1.0, label="nominal")
    ax.set_ylim(-0.05 * resolved_ymax, resolved_ymax)
    ax.set_xlabel("Fixed parameter multiplier")
    loss_label = "Delta fit loss" if admissible_only else "Delta total loss"
    ax.set_ylabel(f"log1p({loss_label})" if profile_scale == "log1p" else loss_label)
    suffix = "admissible points only" if admissible_only else f"{profile_scale} scale"
    ax.set_title(f"Profile likelihood: {parameter} ({suffix})")
    ax.grid(alpha=0.25, linewidth=0.8)
    ax.legend(loc="best", frameon=True, framealpha=0.9)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_perturbed_output_curves(
    scenario,
    reference,
    parameter: str,
    multipliers: np.ndarray,
    outputs: list[str],
    output_path: Path,
) -> Path:
    """Plot model output trajectories under direct parameter perturbations."""

    apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    base_params = default_parameters()
    cmap = plt.get_cmap("viridis")
    ncols = 3
    nrows = int(np.ceil(len(outputs) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.0 * ncols, 3.1 * nrows),
        sharex=True,
        constrained_layout=True,
    )
    axes_flat = np.asarray(axes).ravel()

    perturbed_results = []
    for multiplier in multipliers:
        params = dict(base_params)
        params[parameter] = float(base_params[parameter]) * float(multiplier)
        try:
            perturbed_results.append((float(multiplier), run_simulation(scenario, params)))
        except Exception:
            perturbed_results.append((float(multiplier), None))

    for ax, output in zip(axes_flat, outputs):
        index = STATE_INDEX[output]
        ax.plot(
            reference.t,
            reference.y[:, index],
            color="black",
            linewidth=2.4,
            label="standard",
            zorder=5,
        )
        for draw_index, (multiplier, result) in enumerate(perturbed_results):
            if result is None:
                continue
            color = cmap(draw_index / max(1, len(perturbed_results) - 1))
            ax.plot(
                result.t,
                result.y[:, index],
                color=color,
                linewidth=1.25,
                alpha=0.82,
                label=f"{multiplier:g}x" if output == outputs[0] else None,
            )
        ax.set_title(output)
        ax.set_ylabel(output)
        ax.grid(alpha=0.25, linewidth=0.8)
        ax.margins(x=0.01)

    for ax in axes_flat[len(outputs) :]:
        ax.axis("off")
    for ax in axes_flat[-ncols:]:
        if ax.has_data():
            ax.set_xlabel("Time (days)")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.035), ncol=min(6, len(labels)))
    fig.suptitle(f"Perturbed output trajectories: {parameter}", y=1.08)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_combined_profiles(
    profiles: list[pd.DataFrame],
    output_path: Path,
    y_max: str | float = "auto",
    cutoff: float = 1.92,
    admissible_only: bool = False,
    profile_scale: str = "log1p",
) -> Path:
    apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 6), constrained_layout=True)
    visible_values = []
    for profile in profiles:
        profile = profile.sort_values("fixed_multiplier")
        view = profile_view(profile, admissible_only).sort_values("fixed_multiplier")
        if view.empty:
            continue
        parameter = str(profile["profile_parameter"].iloc[0])
        x = view["fixed_multiplier"].to_numpy(dtype=float)
        y = view["view_delta_loss"].to_numpy(dtype=float)
        y_plot = transform_profile_y(y, profile_scale)
        visible_values.extend(y_plot.tolist())
        if PchipInterpolator is not None and len(view) >= 3:
            x_smooth = np.linspace(float(x.min()), float(x.max()), 240)
            y_smooth = PchipInterpolator(x, y_plot)(x_smooth)
            ax.plot(x_smooth, y_smooth, linewidth=1.9, label=parameter)
            ax.scatter(x, y_plot, s=18)
        elif len(view) >= 2:
            ax.plot(x, y_plot, marker="o", linewidth=1.7, label=parameter)
        else:
            ax.scatter(x, y_plot, s=18, label=parameter)
    cutoff_plot = float(transform_profile_y(np.array([cutoff]), profile_scale)[0])
    ax.axhline(cutoff_plot, linestyle="--", color="0.25", linewidth=1.0, label="approx. 95% cutoff")
    ax.axvline(1.0, linestyle=":", color="0.35", linewidth=1.0)
    resolved_ymax = (
        float(y_max)
        if str(y_max).strip().lower() not in {"auto", "adaptive"}
        else max(1.15 * float(np.nanmax(visible_values)), 1.35 * cutoff_plot)
    )
    ax.set_ylim(-0.05 * resolved_ymax, resolved_ymax)
    ax.set_xlabel("Fixed parameter multiplier")
    loss_label = "Delta fit loss" if admissible_only else "Delta total loss"
    ax.set_ylabel(f"log1p({loss_label})" if profile_scale == "log1p" else loss_label)
    suffix = "admissible points only" if admissible_only else f"{profile_scale} scale"
    ax.set_title(f"Profile likelihood summary ({suffix})")
    ax.grid(alpha=0.25, linewidth=0.8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def synthetic_to_dataframe(data: SyntheticData) -> pd.DataFrame:
    frame = pd.DataFrame(data.values, columns=data.outputs)
    frame.insert(0, "time_days", data.times)
    return frame


def main() -> None:
    args = parse_args()
    y_max_overrides = parse_ymax_overrides(args.profile_ymax_override)
    figure_dir = PROJECT_ROOT / "results" / "figures"
    table_dir = PROJECT_ROOT / "results" / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    profile_parameters, nuisance_parameters = select_parameters(args)
    if not profile_parameters:
        raise SystemExit("No profile parameters selected.")

    scenario = constant_non_lactating(days=args.days, step=args.dt)
    reference = run_simulation(scenario)
    synthetic = make_synthetic_data(
        reference,
        outputs=args.outputs,
        sample_every_days=args.sample_every,
        noise_relative=args.noise_relative,
        seed=args.seed,
    )
    thresholds = AdmissibilityThresholds(
        min_correlation=args.rho_min,
        max_average_difference=args.gamma_max,
        max_norm_difference=args.kappa_max,
        max_shift_days=args.max_shift_days,
        penalty_weight=args.penalty_weight,
    )
    grid = np.linspace(args.grid_low, args.grid_high, args.grid_points)

    synthetic_plot = plot_synthetic_data(synthetic, figure_dir / f"{args.prefix}_synthetic_outputs.png")
    if args.save_synthetic_csv:
        synthetic_to_dataframe(synthetic).to_csv(table_dir / f"{args.prefix}_synthetic_data.csv", index=False)

    profiles = []
    summary_rows = []
    trajectory_plots = []
    for parameter in profile_parameters:
        profile = profile_one_parameter(
            scenario=scenario,
            reference=reference,
            data=synthetic,
            profile_parameter=parameter,
            nuisance_parameters=nuisance_parameters,
            grid_multipliers=grid,
            admissibility_states=args.outputs,
            thresholds=thresholds,
            nuisance_bounds=(args.nuisance_low, args.nuisance_high),
            maxiter=args.maxiter,
        )
        profiles.append(profile)
        summary_rows.append(classify_profile(profile, threshold=args.profile_cutoff, admissible_only=args.admissible_only))
        plot_profile(
            profile,
            figure_dir / f"{args.prefix}_{parameter}_profile.png",
            y_max=y_max_overrides.get(parameter, args.profile_ymax),
            cutoff=args.profile_cutoff,
            admissible_only=args.admissible_only,
            profile_scale=args.profile_scale,
        )
        if should_plot_trajectory(parameter, profile_parameters, args.trajectory_parameter):
            trajectory_plots.append(
                plot_perturbed_output_curves(
                    scenario,
                    reference,
                    parameter,
                    np.asarray(args.trajectory_multipliers, dtype=float),
                    args.outputs,
                    figure_dir / f"{args.prefix}_{parameter}_perturbed_trajectories.png",
                )
            )
        if args.save_tables:
            profile.to_csv(table_dir / f"{args.prefix}_{parameter}_profile.csv", index=False)
        print(f"profiled {parameter}: {summary_rows[-1]['classification']}")

    summary = pd.DataFrame(summary_rows)
    summary_path = table_dir / f"{args.prefix}_summary.csv"
    summary.to_csv(summary_path, index=False)

    combined_plot = plot_combined_profiles(
        profiles,
        figure_dir / f"{args.prefix}_combined_profiles.png",
        y_max=args.profile_ymax,
        cutoff=args.profile_cutoff,
        admissible_only=args.admissible_only,
        profile_scale=args.profile_scale,
    )

    selected = pd.DataFrame(
        {
            "profile_parameter": profile_parameters,
            "used_as_nuisance": [name in nuisance_parameters for name in profile_parameters],
        }
    )
    selected["nuisance_parameters"] = ", ".join(nuisance_parameters)
    selected.to_csv(table_dir / f"{args.prefix}_selected_parameters.csv", index=False)

    print()
    print("Profile likelihood outputs")
    print(f"  synthetic plot: {synthetic_plot}")
    print(f"  combined profiles: {combined_plot}")
    print(f"  perturbed trajectory plots: {len(trajectory_plots)}")
    print(f"  compact summary: {summary_path}")


if __name__ == "__main__":
    main()
