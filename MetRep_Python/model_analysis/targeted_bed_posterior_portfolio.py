"""Targeted BED/posterior portfolio workflow linking diagnostics and broad priors.

No ODE simulations are run. Ranking stability is estimated from ODE-confirmed
admissible banks, while posterior reweighting plots use the broad +/-5% bank.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import gaussian_kde, rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "MetRep_Python"
sys.path.insert(0, str(PYTHON_ROOT))

from model_definition.parameters import PARAMETERS


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
DEFAULT_ENRICHED_FORMATTED = ROOT / "local_data/phd_bed_bank_5pct_50k_glucagon_smc_enriched"
DEFAULT_GLOBAL = ROOT / "local_outputs/global_sensitivity"
DEFAULT_UNCERTAINTY = ROOT / "local_outputs/uncertainty"
DEFAULT_PROFILE = ROOT / "results_final/tables"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run targeted BED/posterior demo from saved banks.")
    parser.add_argument("--broad-prior-dir", type=Path, required=True)
    parser.add_argument("--enriched-bank-dir", type=Path, required=True)
    parser.add_argument("--global-sensitivity-dir", type=Path, default=DEFAULT_GLOBAL)
    parser.add_argument("--uncertainty-dir", type=Path, default=DEFAULT_UNCERTAINTY)
    parser.add_argument("--profile-dir", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-observation-scenarios", type=int, default=4)
    parser.add_argument("--kde-grid-size", type=int, default=500)
    return parser.parse_args()


def mark(step: int, text: str) -> None:
    print(f"[STEP {step}/11] {text}", flush=True)


def resolve_existing(path: Path, fallback: Path | None = None) -> Path:
    if path.exists():
        return path
    if fallback is not None and fallback.exists():
        return fallback
    raise FileNotFoundError(path)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def infer_observables(outputs: pd.DataFrame) -> list[str]:
    found = sorted({m.group(1) for c in outputs.columns if (m := re.match(r"^(.+)_day_(\d+)$", c))})
    ordered = [name for name in DEFAULT_OBSERVABLE_ORDER if name in found]
    ordered.extend(name for name in found if name not in ordered)
    return ordered


def stored_days(outputs: pd.DataFrame, biomarker: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(biomarker)}_day_(\d+)$")
    return sorted(int(match.group(1)) for col in outputs.columns if (match := pattern.match(col)))


def day_columns(day: int, biomarkers: list[str]) -> list[str]:
    return [f"{name}_day_{day}" for name in biomarkers]


def load_bank(bank_dir: Path, require_admissible: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    parameters = read_csv(bank_dir / "prior_parameter_samples.csv")
    outputs = read_csv(bank_dir / "ode_output_features.csv")
    admissibility = read_csv(bank_dir / "admissibility.csv")
    metadata_path = bank_dir / "bank_metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    if not (len(parameters) == len(outputs) == len(admissibility)):
        raise ValueError(f"Bank has inconsistent row counts: {bank_dir}")
    finite = parameters.notna().all(axis=1) & outputs.notna().all(axis=1)
    if require_admissible:
        finite = finite & admissibility["admissible"].astype(bool)
    return (
        parameters.loc[finite].reset_index(drop=True),
        outputs.loc[finite].reset_index(drop=True),
        admissibility.loc[finite].reset_index(drop=True),
        metadata,
    )


def resolve_enriched_bank(path: Path) -> Path:
    if (path / "prior_parameter_samples.csv").exists():
        return path
    if DEFAULT_ENRICHED_FORMATTED.exists():
        return DEFAULT_ENRICHED_FORMATTED
    raise FileNotFoundError(
        f"{path} is not a formatted bank. Build/use {DEFAULT_ENRICHED_FORMATTED} first."
    )


def load_profile_classes(profile_dir: Path, targets: list[str]) -> pd.DataFrame:
    for filename in ["combined_parameter_summary.csv", "profile_50d_balanced_relaxed_selected_parameters.csv"]:
        path = profile_dir / filename
        if path.exists():
            table = pd.read_csv(path)
            if "parameter" in table.columns:
                cols = [c for c in ["parameter", "identifiability_class", "local_sensitivity_rank", "local_sensitivity_max_abs_auc", "local_sensitivity_main_output"] if c in table.columns]
                out = table[cols].copy()
                return out[out["parameter"].isin(targets)].drop_duplicates("parameter")
    return pd.DataFrame({"parameter": targets, "identifiability_class": ["Selected"] * len(targets)})


def load_global_links(global_dir: Path, targets: list[str]) -> pd.DataFrame:
    rows = []
    for filenames, score_col, value_col, method in [
        (("global_sensitivity_prcc_auc.csv", "global_sensitivity_98x9_prcc.csv"), "abs_prcc", "prcc", "PRCC"),
        (("global_sensitivity_spearman_auc.csv", "global_sensitivity_98x9_spearman.csv"), "abs_spearman_rho", "spearman_rho", "Spearman"),
    ]:
        path = next((global_dir / filename for filename in filenames if (global_dir / filename).exists()), None)
        if path is None:
            continue
        table = pd.read_csv(path)
        required = {"parameter", "biomarker", score_col, value_col}
        if not required.issubset(table.columns):
            continue
        table = table[table["parameter"].isin(targets)].copy()
        for parameter, group in table.groupby("parameter"):
            ranked = group.sort_values(score_col, ascending=False).head(5)
            for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
                rows.append(
                    {
                        "parameter": parameter,
                        "method": method,
                        "rank": rank,
                        "biomarker": row["biomarker"],
                        "association": float(row[value_col]),
                        "abs_association": float(row[score_col]),
                    }
                )
    return pd.DataFrame(rows, columns=["parameter", "method", "rank", "biomarker", "association", "abs_association"])


def load_uncertainty_windows(uncertainty_dir: Path, global_links: pd.DataFrame) -> pd.DataFrame:
    path = uncertainty_dir / "trajectory_quantiles.csv"
    if not path.exists() or global_links.empty:
        return pd.DataFrame()
    quantiles = pd.read_csv(path)
    biomarkers = sorted(global_links["biomarker"].unique())
    rows = []
    for biomarker in biomarkers:
        subset = quantiles[quantiles["biomarker"].eq(biomarker)].sort_values("interval_width", ascending=False).head(5)
        for rank, (_, row) in enumerate(subset.iterrows(), start=1):
            rows.append(
                {
                    "biomarker": biomarker,
                    "rank": rank,
                    "day": int(row["day"]),
                    "interval_width": float(row["interval_width"]),
                    "samples": int(row["samples"]),
                }
            )
    return pd.DataFrame(rows)


def gaussian_rank_mi(x: np.ndarray, y: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 10 or np.std(x[finite]) <= 0 or np.std(y[finite]) <= 0:
        return 0.0
    rho, _ = spearmanr(x[finite], y[finite])
    if not np.isfinite(rho):
        return 0.0
    rho = float(np.clip(rho, -0.999999, 0.999999))
    return float(-0.5 * np.log1p(-(rho * rho)))


def targeted_rankings(parameters: pd.DataFrame, outputs: pd.DataFrame, targets: list[str], biomarkers: list[str]) -> pd.DataFrame:
    rows = []
    days = sorted(set.intersection(*(set(stored_days(outputs, b)) for b in biomarkers)))
    for parameter in targets:
        x = parameters[parameter].to_numpy(float)
        for day in days:
            day_scores = []
            for biomarker in biomarkers:
                y = outputs[f"{biomarker}_day_{day}"].to_numpy(float)
                score = gaussian_rank_mi(x, y)
                day_scores.append(score)
                rows.append(
                    {
                        "parameter": parameter,
                        "day": day,
                        "biomarker": biomarker,
                        "single_biomarker_mi_proxy": score,
                    }
                )
            rows.append(
                {
                    "parameter": parameter,
                    "day": day,
                    "biomarker": "ALL_OBSERVABLES",
                    "single_biomarker_mi_proxy": float(np.sum(day_scores)),
                }
            )
    table = pd.DataFrame(rows)
    table["rank_within_parameter"] = table.groupby("parameter")["single_biomarker_mi_proxy"].rank(method="first", ascending=False).astype(int)
    return table.sort_values(["parameter", "rank_within_parameter"]).reset_index(drop=True)


def select_targets(profile: pd.DataFrame, parameters: pd.DataFrame) -> list[str]:
    available = [name for name in DEFAULT_TARGETS if name in parameters.columns]
    if len(available) == len(DEFAULT_TARGETS):
        return available
    extras = [name for name in profile["parameter"].tolist() if name in parameters.columns and name not in available]
    return (available + extras)[:9]


def load_noise(bank_dir: Path, biomarkers: list[str]) -> dict[str, float]:
    table = pd.read_csv(bank_dir / "measurement_noise.csv")
    noise = dict(zip(table["biomarker"], table["abs_error"]))
    return {name: max(float(noise[name]), 1e-12) for name in biomarkers if name in noise}


def posterior_weights(outputs: pd.DataFrame, day: int, biomarkers: list[str], noise: dict[str, float], observation_index: int) -> np.ndarray:
    y = outputs[day_columns(day, biomarkers)].to_numpy(float)
    z = y[observation_index]
    sigma = np.array([noise[name] for name in biomarkers], dtype=float)
    log_like = -0.5 * np.sum(((z[None, :] - y) / sigma[None, :]) ** 2, axis=1)
    weights = np.exp(log_like - logsumexp(log_like))
    return weights


def posterior_weights_observations(
    outputs: pd.DataFrame,
    observations: list[tuple[str, int]],
    noise: dict[str, float],
    observation_index: int,
) -> np.ndarray:
    columns = [f"{biomarker}_day_{day}" for biomarker, day in observations]
    available = [(biomarker, day, column) for (biomarker, day), column in zip(observations, columns) if column in outputs.columns and biomarker in noise]
    if not available:
        raise ValueError("No valid biomarker-day observations available for posterior weighting")
    columns = [column for _, _, column in available]
    y = outputs[columns].to_numpy(float)
    z = y[observation_index]
    sigma = np.array([noise[biomarker] for biomarker, _, _ in available], dtype=float)
    log_like = -0.5 * np.sum(((z[None, :] - y) / sigma[None, :]) ** 2, axis=1)
    return np.exp(log_like - logsumexp(log_like))


def top_guided_biomarkers(global_links: pd.DataFrame, parameter: str, limit: int = 3) -> list[dict[str, object]]:
    subset = global_links[global_links["parameter"].eq(parameter)].copy()
    if subset.empty:
        return []
    wide = subset.pivot_table(index="biomarker", columns="method", values=["association", "abs_association"], aggfunc="max")
    rows = []
    for biomarker in wide.index:
        prcc = subset[(subset["biomarker"].eq(biomarker)) & (subset["method"].eq("PRCC"))]
        spearman = subset[(subset["biomarker"].eq(biomarker)) & (subset["method"].eq("Spearman"))]
        prcc_value = float(prcc["association"].iloc[0]) if not prcc.empty else np.nan
        spearman_value = float(spearman["association"].iloc[0]) if not spearman.empty else np.nan
        score = np.nanmax([abs(prcc_value), abs(spearman_value)])
        rows.append({"biomarker": biomarker, "PRCC": prcc_value, "Spearman": spearman_value, "score": score})
    rows = sorted(rows, key=lambda item: item["score"], reverse=True)
    return rows[:limit]


def guided_candidate_day_table(
    parameter: str,
    biomarkers: list[str],
    uncertainty_windows: pd.DataFrame,
    rankings: pd.DataFrame,
    outputs: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for biomarker in biomarkers:
        wins = uncertainty_windows[uncertainty_windows["biomarker"].eq(biomarker)].head(5)
        if wins.empty:
            candidate_days = stored_days(outputs, biomarker)
            uncertainty_by_day = {day: np.nan for day in candidate_days}
        else:
            candidate_days = [int(day) for day in wins["day"].tolist()]
            uncertainty_by_day = dict(zip(wins["day"].astype(int), wins["interval_width"].astype(float)))
        subset = rankings[
            rankings["parameter"].eq(parameter)
            & rankings["biomarker"].eq(biomarker)
            & rankings["day"].isin(candidate_days)
        ].copy()
        if subset.empty:
            subset = rankings[rankings["parameter"].eq(parameter) & rankings["biomarker"].eq(biomarker)].copy()
        for _, row in subset.iterrows():
            day = int(row["day"])
            rows.append(
                {
                    "parameter": parameter,
                    "biomarker": biomarker,
                    "candidate_day": day,
                    "uncertainty_metric": float(uncertainty_by_day.get(day, np.nan)),
                    "MI_proxy_or_MI": float(row["single_biomarker_mi_proxy"]),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["parameter", "biomarker", "candidate_day", "uncertainty_metric", "MI_proxy_or_MI", "rank"])
    out = out.sort_values(["parameter", "biomarker", "MI_proxy_or_MI"], ascending=[True, True, False])
    out["rank"] = out.groupby(["parameter", "biomarker"])["MI_proxy_or_MI"].rank(method="first", ascending=False).astype(int)
    return out


def weighted_quantile(values: np.ndarray, weights: np.ndarray, qs: list[float]) -> np.ndarray:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    cdf /= cdf[-1]
    return np.interp(qs, cdf, values)


def narrowing(values: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    prior = np.full(len(values), 1.0 / len(values))
    pm = float(np.sum(prior * values))
    wm = float(np.sum(weights * values))
    pv = float(np.sum(prior * (values - pm) ** 2))
    wv = float(np.sum(weights * (values - wm) ** 2))
    p05, p95 = weighted_quantile(values, prior, [0.05, 0.95])
    q05, q95 = weighted_quantile(values, weights, [0.05, 0.95])
    return {
        "ess": float(1.0 / np.sum(weights * weights)),
        "ess_fraction": float((1.0 / np.sum(weights * weights)) / len(weights)),
        "prior_variance": pv,
        "posterior_variance": wv,
        "variance_reduction_fraction": 1.0 - wv / pv if pv > 0 else np.nan,
        "prior_90_width": float(p95 - p05),
        "posterior_90_width": float(q95 - q05),
        "interval_reduction_fraction": 1.0 - float(q95 - q05) / float(p95 - p05) if p95 > p05 else np.nan,
        "weight_warning": "uniform" if (1.0 / np.sum(weights * weights)) / len(weights) > 0.8 else ("collapsed" if (1.0 / np.sum(weights * weights)) / len(weights) < 0.02 else "informative"),
    }


def kde_curve(values: np.ndarray, weights: np.ndarray, grid: np.ndarray) -> np.ndarray:
    kde = gaussian_kde(values, weights=weights)
    y = kde(grid)
    area = np.trapezoid(y, grid)
    return y / area if area > 0 else y


def curve_rows(parameter: str, label: str, values: np.ndarray, weights: np.ndarray, grid_size: int) -> pd.DataFrame:
    span = float(values.max() - values.min())
    margin = max(0.20 * span, 1e-12)
    grid = np.linspace(values.min() - margin, values.max() + margin, grid_size)
    density = kde_curve(values, weights, grid)
    return pd.DataFrame({"parameter": parameter, "curve": label, "grid": grid, "density": density})


def plot_curve_grid(curves: pd.DataFrame, figure_dir: Path, filename: str, title: str) -> None:
    targets = list(dict.fromkeys(curves["parameter"]))
    fig, axes = plt.subplots(3, 3, figsize=(18, 13), constrained_layout=True, squeeze=False)
    for ax, parameter in zip(axes.ravel(), targets):
        subset = curves[curves["parameter"].eq(parameter)]
        for label, group in subset.groupby("curve", sort=False):
            ax.plot(group["grid"], group["density"], linewidth=1.5, label=label)
        ax.set_title(parameter, fontsize=10)
        ax.set_xlabel("Parameter value")
        ax.set_ylabel("Density")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    for ax in axes.ravel()[len(targets):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=15)
    fig.savefig(figure_dir / filename, dpi=220)
    plt.close(fig)


def plot_rankings(rankings: pd.DataFrame, selection: pd.DataFrame, figure_dir: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(18, 12), constrained_layout=True, squeeze=False)
    for ax, (parameter, group) in zip(axes.ravel(), rankings[rankings["biomarker"].eq("ALL_OBSERVABLES")].groupby("parameter")):
        ordered = group.sort_values("day")
        ax.plot(ordered["day"], ordered["single_biomarker_mi_proxy"], marker="o", linewidth=1.2)
        best_day = int(selection.loc[selection["parameter"].eq(parameter), "best_day"].iloc[0])
        ax.axvline(best_day, color="black", linestyle="--", linewidth=1)
        ax.set_title(parameter, fontsize=10)
        ax.set_xlabel("Day")
        ax.set_ylabel("Targeted MI proxy")
        ax.grid(alpha=0.25)
    for ax in axes.ravel()[selection["parameter"].nunique():]:
        ax.axis("off")
    fig.suptitle("Targeted BED ranking summary", fontsize=15)
    fig.savefig(figure_dir / "bed_targeted_ranking_summary.png", dpi=220)
    plt.close(fig)


def plot_selection_table(selection: pd.DataFrame, figure_dir: Path) -> None:
    shown = selection[["parameter", "profile_class", "top_biomarkers", "best_day", "best_species", "reason_selected"]].copy()
    fig, ax = plt.subplots(figsize=(16, 0.7 * len(shown) + 2), constrained_layout=True)
    ax.axis("off")
    table = ax.table(cellText=shown.values, colLabels=shown.columns, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    fig.savefig(figure_dir / "target_parameter_selection_table.png", dpi=220)
    plt.close(fig)


def plot_guided_biomarker_bars(guided_biomarkers: pd.DataFrame, targets: list[str], figure_dir: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 10), constrained_layout=True, squeeze=False)
    for ax, parameter in zip(axes.ravel(), targets):
        group = guided_biomarkers[guided_biomarkers["parameter"].eq(parameter)].sort_values("rank", ascending=False)
        ax.barh(group["biomarker"], group["MI_parameter_biomarker"])
        ax.set_title(parameter, fontsize=10)
        ax.set_xlabel("GSA association / MI proxy")
        ax.grid(axis="x", alpha=0.25)
    for ax in axes.ravel()[len(targets):]:
        ax.axis("off")
    fig.suptitle("Parameter-specific GSA-selected biomarkers for guided BED", fontsize=15)
    fig.savefig(figure_dir / "bed_guided_parameter_biomarker_mi_bars.png", dpi=220)
    plt.close(fig)


def plot_guided_day_curves(guided_days: pd.DataFrame, targets: list[str], figure_dir: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 10), constrained_layout=True, squeeze=False)
    x_min = int(guided_days["candidate_day"].min()) if not guided_days.empty else 54
    x_max = int(guided_days["candidate_day"].max()) if not guided_days.empty else 89
    for ax, parameter in zip(axes.ravel(), targets):
        subset = guided_days[guided_days["parameter"].eq(parameter)]
        for biomarker, group in subset.groupby("biomarker", sort=False):
            ordered = group.sort_values("candidate_day")
            ax.plot(
                ordered["candidate_day"],
                ordered["MI_proxy_or_MI"],
                marker="o",
                linestyle="None",
                markersize=5,
                label=biomarker,
            )
        ax.set_title(parameter, fontsize=10)
        ax.set_xlabel("Candidate day")
        ax.set_ylabel("MI proxy")
        ax.set_xlim(x_min - 1, x_max + 1)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=6, loc="best")
    for ax in axes.ravel()[len(targets):]:
        ax.axis("off")
    fig.suptitle("Parameter-specific candidate-day MI points over uncertainty-selected windows", fontsize=15)
    fig.savefig(figure_dir / "bed_guided_day_mi_curves.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    mark(1, "Loading broad +/-5% prior bank and enriched admissible bank")
    enriched_dir = resolve_enriched_bank(args.enriched_bank_dir)
    broad_p, broad_y, _, broad_meta = load_bank(args.broad_prior_dir, require_admissible=False)
    enriched_p, enriched_y, _, enriched_meta = load_bank(enriched_dir, require_admissible=True)
    finite_count = len(broad_p)

    mark(2, "Loading profile-likelihood parameter classes")
    profile_seed = pd.read_csv(args.profile_dir / "combined_parameter_summary.csv") if (args.profile_dir / "combined_parameter_summary.csv").exists() else pd.DataFrame({"parameter": DEFAULT_TARGETS})
    targets = select_targets(profile_seed, broad_p)
    profile = load_profile_classes(args.profile_dir, targets)

    mark(3, "Loading global sensitivity results")
    global_dir = resolve_existing(args.global_sensitivity_dir, DEFAULT_GLOBAL)
    global_links = load_global_links(global_dir, targets)

    mark(4, "Loading uncertainty propagation results")
    uncertainty_dir = resolve_existing(args.uncertainty_dir, DEFAULT_UNCERTAINTY)
    uncertainty_windows = load_uncertainty_windows(uncertainty_dir, global_links)

    mark(5, "Selecting representative target parameters")
    observables = infer_observables(enriched_y)
    noise = load_noise(args.broad_prior_dir, observables)
    selection_rows = []
    for target in targets:
        cls = profile.loc[profile["parameter"].eq(target), "identifiability_class"]
        top_biomarkers = global_links[(global_links["parameter"].eq(target)) & (global_links["method"].eq("PRCC"))].sort_values("rank")["biomarker"].head(5).tolist()
        if not top_biomarkers:
            top_biomarkers = observables
        selection_rows.append(
            {
                "parameter": target,
                "profile_class": cls.iloc[0] if not cls.empty else "Selected",
                "top_biomarkers": "|".join(top_biomarkers),
                "reason_selected": "representative profile-likelihood class; linked to global sensitivity and BED",
            }
        )
    selection = pd.DataFrame(selection_rows)

    mark(6, "Selecting candidate biomarkers/days from sensitivity and uncertainty")
    candidate_rows = []
    for _, row in selection.iterrows():
        for biomarker in str(row["top_biomarkers"]).split("|"):
            wins = uncertainty_windows[uncertainty_windows["biomarker"].eq(biomarker)].head(3)
            if wins.empty:
                for day in stored_days(enriched_y, biomarker):
                    candidate_rows.append({"parameter": row["parameter"], "biomarker": biomarker, "candidate_day": day, "source": "stored_day"})
            else:
                for _, w in wins.iterrows():
                    candidate_rows.append({"parameter": row["parameter"], "biomarker": biomarker, "candidate_day": int(w["day"]), "source": "uncertainty_window"})
    candidate_table = pd.DataFrame(candidate_rows).drop_duplicates()

    mark(7, "Running BED ranking and stability checks")
    enriched_rankings = targeted_rankings(enriched_p, enriched_y, targets, observables)
    original_keep = pd.read_csv(args.broad_prior_dir / "admissibility.csv")["admissible"].astype(bool) & broad_y.notna().all(axis=1)
    original_rankings = targeted_rankings(broad_p.loc[original_keep].reset_index(drop=True), broad_y.loc[original_keep].reset_index(drop=True), targets, observables)
    best_rows = []
    for target in targets:
        target_rank = enriched_rankings[enriched_rankings["parameter"].eq(target)]
        best_day = int(target_rank[target_rank["biomarker"].eq("ALL_OBSERVABLES")].iloc[0]["day"])
        species_order = (
            target_rank[(target_rank["day"].eq(best_day)) & (~target_rank["biomarker"].eq("ALL_OBSERVABLES"))]
            .sort_values("single_biomarker_mi_proxy", ascending=False)["biomarker"]
            .tolist()
        )
        low_day = int(target_rank[target_rank["biomarker"].eq("ALL_OBSERVABLES")].sort_values("single_biomarker_mi_proxy").iloc[0]["day"])
        best_rows.append({"parameter": target, "best_day": best_day, "low_day": low_day, "best_species": "|".join(species_order), "top_species": species_order[0]})
    best_table = pd.DataFrame(best_rows)
    selection = selection.merge(best_table, on="parameter", how="left")
    selection["candidate_bed_biomarkers_days"] = selection["best_species"] + " @ day " + selection["best_day"].astype(str)
    stability_rows = []
    for target in targets:
        old_best = original_rankings[(original_rankings["parameter"].eq(target)) & (original_rankings["biomarker"].eq("ALL_OBSERVABLES"))].iloc[0]
        new_best = enriched_rankings[(enriched_rankings["parameter"].eq(target)) & (enriched_rankings["biomarker"].eq("ALL_OBSERVABLES"))].iloc[0]
        stability_rows.append(
            {
                "parameter": target,
                "original_best_day": int(old_best["day"]),
                "enriched_best_day": int(new_best["day"]),
                "same_best_day": int(old_best["day"]) == int(new_best["day"]),
                "original_best_score": float(old_best["single_biomarker_mi_proxy"]),
                "enriched_best_score": float(new_best["single_biomarker_mi_proxy"]),
            }
        )
    stability = pd.DataFrame(stability_rows)

    mark(8, "Computing independent posterior updates from broad prior")
    observation_index = int(rng.integers(0, len(broad_y)))
    prior_weights = np.full(len(broad_p), 1.0 / len(broad_p))
    independent_rows = []
    independent_curves = []
    for _, row in selection.iterrows():
        parameter = row["parameter"]
        values = broad_p[parameter].to_numpy(float)
        independent_curves.append(curve_rows(parameter, "Broad +/-5% prior", values, prior_weights, args.kde_grid_size))
        rank = enriched_rankings[(enriched_rankings["parameter"].eq(parameter)) & (~enriched_rankings["biomarker"].eq("ALL_OBSERVABLES"))].head(args.top_observation_scenarios)
        for scenario_i, (_, r) in enumerate(rank.iterrows(), start=1):
            biomarker = r["biomarker"]
            day = int(r["day"])
            weights = posterior_weights(broad_y, day, [biomarker], noise, observation_index)
            label = f"obs {scenario_i}: {biomarker} day {day}"
            independent_curves.append(curve_rows(parameter, label, values, weights, args.kde_grid_size))
            independent_rows.append({"parameter": parameter, "scenario": label, "day": day, "biomarkers": biomarker, **narrowing(values, weights)})
    independent_curves = pd.concat(independent_curves, ignore_index=True)

    mark(9, "Computing cumulative biomarker posterior updates")
    cumulative_rows = []
    cumulative_curves = []
    for _, row in selection.iterrows():
        parameter = row["parameter"]
        values = broad_p[parameter].to_numpy(float)
        best_day = int(row["best_day"])
        order = str(row["best_species"]).split("|")
        cumulative_curves.append(curve_rows(parameter, "Broad +/-5% prior", values, prior_weights, args.kde_grid_size))
        for count in range(1, len(order) + 1):
            biomarkers = order[:count]
            weights = posterior_weights(broad_y, best_day, biomarkers, noise, observation_index)
            label = f"best {count}: {' + '.join(biomarkers)}"
            if count in {1, 2, 3, len(order)}:
                cumulative_curves.append(curve_rows(parameter, label, values, weights, args.kde_grid_size))
            cumulative_rows.append({"parameter": parameter, "day": best_day, "n_biomarkers": count, "biomarkers": "|".join(biomarkers), **narrowing(values, weights)})
    cumulative_curves = pd.concat(cumulative_curves, ignore_index=True)

    guided_biomarker_rows = []
    guided_day_frames = []
    guided_trace_rows = []
    guided_scenario_rows = []
    guided_curves = []
    for _, row in selection.iterrows():
        parameter = row["parameter"]
        profile_class = row["profile_class"]
        values = broad_p[parameter].to_numpy(float)
        guided_curves.append(curve_rows(parameter, "Broad +/-5% prior", values, prior_weights, args.kde_grid_size))
        selected = top_guided_biomarkers(global_links, parameter, limit=3)
        if not selected:
            selected = [{"biomarker": b, "PRCC": np.nan, "Spearman": np.nan, "score": np.nan} for b in observables[:3]]
        selected_biomarkers = [str(item["biomarker"]) for item in selected]
        for rank, item in enumerate(selected, start=1):
            guided_biomarker_rows.append(
                {
                    "parameter": parameter,
                    "profile_class": profile_class,
                    "biomarker": item["biomarker"],
                    "PRCC": item["PRCC"],
                    "Spearman": item["Spearman"],
                    "MI_parameter_biomarker": item["score"],
                    "rank": rank,
                }
            )
        day_candidates = guided_candidate_day_table(parameter, selected_biomarkers, uncertainty_windows, enriched_rankings, enriched_y)
        guided_day_frames.append(day_candidates)
        scenario_observations: list[tuple[str, int]] = []
        for rank, item in enumerate(selected, start=1):
            biomarker = str(item["biomarker"])
            chosen = day_candidates[day_candidates["biomarker"].eq(biomarker)].sort_values("rank").head(1)
            if chosen.empty:
                continue
            selected_day = int(chosen["candidate_day"].iloc[0])
            uncertainty_metric = float(chosen["uncertainty_metric"].iloc[0])
            day_mi = float(chosen["MI_proxy_or_MI"].iloc[0])
            scenario_observations.append((biomarker, selected_day))
            guided_trace_rows.append(
                {
                    "parameter": parameter,
                    "profile_class": profile_class,
                    "selected_biomarker": biomarker,
                    "biomarker_rank": rank,
                    "PRCC": item["PRCC"],
                    "Spearman": item["Spearman"],
                    "uncertainty_metric": uncertainty_metric,
                    "selected_day": selected_day,
                    "day_MI": day_mi,
                    "scenario_name": f"guided_{parameter}",
                    "used_in_plot": True,
                }
            )
        if scenario_observations:
            weights = posterior_weights_observations(broad_y, scenario_observations, noise, observation_index)
            label = "guided: " + " + ".join(f"{biomarker} d{day}" for biomarker, day in scenario_observations)
            guided_curves.append(curve_rows(parameter, label, values, weights, args.kde_grid_size))
            narrow = narrowing(values, weights)
            for trace_row in guided_trace_rows:
                if trace_row["parameter"] == parameter:
                    trace_row["posterior_width_reduction_percent"] = 100.0 * narrow["interval_reduction_fraction"]
                    trace_row["ESS"] = narrow["ess"]
            guided_scenario_rows.append(
                {
                    "parameter": parameter,
                    "profile_class": profile_class,
                    "scenario": f"guided_{parameter}",
                    "biomarker_day_combination": "|".join(f"{biomarker}_day_{day}" for biomarker, day in scenario_observations),
                    "reason_selected": "top parameter-specific GSA biomarkers; uncertainty-window candidate days; final days ranked by BED MI proxy",
                }
            )
    guided_biomarkers = pd.DataFrame(guided_biomarker_rows)
    guided_days = pd.concat(guided_day_frames, ignore_index=True) if guided_day_frames else pd.DataFrame()
    guided_trace = pd.DataFrame(guided_trace_rows)
    guided_scenarios = pd.DataFrame(guided_scenario_rows)
    guided_curves = pd.concat(guided_curves, ignore_index=True)

    mark(10, "Generating final plots and tables")
    selection.to_csv(args.output_dir / "target_parameter_selection.csv", index=False)
    global_links.to_csv(args.output_dir / "global_sensitivity_links_for_targets.csv", index=False)
    uncertainty_windows.to_csv(args.output_dir / "uncertainty_windows_for_targets.csv", index=False)
    enriched_rankings.to_csv(args.output_dir / "bed_rankings_targeted.csv", index=False)
    pd.DataFrame(independent_rows).to_csv(args.output_dir / "posterior_narrowing_independent_scenarios.csv", index=False)
    pd.DataFrame(cumulative_rows).to_csv(args.output_dir / "posterior_narrowing_cumulative_biomarkers.csv", index=False)
    guided_trace.to_csv(args.output_dir / "bed_guided_update_traceability.csv", index=False)
    guided_biomarkers.to_csv(args.output_dir / "bed_guided_parameter_biomarker_mi.csv", index=False)
    guided_days.to_csv(args.output_dir / "bed_guided_parameter_biomarker_day_mi.csv", index=False)
    guided_scenarios.to_csv(args.output_dir / "bed_guided_selected_observation_scenarios.csv", index=False)
    stability.to_csv(args.output_dir / "bed_ranking_stability.csv", index=False)
    final_tables = ROOT / "results_final/tables"
    final_figures = ROOT / "results_final/figures"
    final_tables.mkdir(parents=True, exist_ok=True)
    final_figures.mkdir(parents=True, exist_ok=True)
    selection.to_csv(final_tables / "bayesian_bed_representative_targets.csv", index=False)
    global_links.to_csv(final_tables / "bed_targeted_parameter_biomarker_time_links.csv", index=False)
    uncertainty_windows.to_csv(final_tables / "bed_targeted_gsa_uncertainty_observation_scenarios.csv", index=False)
    enriched_rankings.to_csv(final_tables / "bed_mi_rankings_12k_stability.csv", index=False)
    pd.DataFrame(independent_rows).to_csv(final_tables / "bed_targeted_posterior_narrowing_independent.csv", index=False)
    pd.DataFrame(cumulative_rows).to_csv(final_tables / "bed_targeted_posterior_narrowing_cumulative.csv", index=False)
    guided_trace.to_csv(final_tables / "bed_guided_update_traceability.csv", index=False)
    guided_biomarkers.to_csv(final_tables / "bed_guided_parameter_biomarker_mi.csv", index=False)
    guided_days.to_csv(final_tables / "bed_guided_parameter_biomarker_day_mi.csv", index=False)
    guided_scenarios.to_csv(final_tables / "bed_guided_selected_observation_scenarios.csv", index=False)
    guided_trace.to_csv(final_tables / "bed_guided_posterior_narrowing_parameter_specific.csv", index=False)
    stability.to_csv(final_tables / "bed_mi_stability_original_vs_12k.csv", index=False)
    pd.DataFrame(cumulative_rows + independent_rows).to_csv(final_tables / "bed_targeted_posterior_weight_diagnostics.csv", index=False)
    plot_selection_table(selection, args.figure_dir)
    plot_rankings(enriched_rankings, selection, args.figure_dir)
    plot_guided_biomarker_bars(guided_biomarkers, targets, args.figure_dir)
    plot_guided_day_curves(guided_days, targets, args.figure_dir)
    plot_curve_grid(independent_curves, args.figure_dir, "posterior_independent_observation_scenarios.png", "Independent observation-scenario posteriors from broad +/-5% prior")
    plot_curve_grid(cumulative_curves, args.figure_dir, "posterior_cumulative_biomarkers.png", "Global cumulative biomarker update from broad +/-5% prior")
    plot_curve_grid(guided_curves, args.figure_dir, "posterior_guided_parameter_specific.png", "Parameter-specific GSA + uncertainty + MI guided posterior updates")
    high_low_curves = []
    for _, row in selection.iterrows():
        parameter = row["parameter"]
        values = broad_p[parameter].to_numpy(float)
        order = str(row["best_species"]).split("|")
        high_w = posterior_weights(broad_y, int(row["best_day"]), order, noise, observation_index)
        low_w = posterior_weights(broad_y, int(row["low_day"]), order, noise, observation_index)
        high_low_curves.append(curve_rows(parameter, "Broad +/-5% prior", values, prior_weights, args.kde_grid_size))
        high_low_curves.append(curve_rows(parameter, f"high day {int(row['best_day'])}", values, high_w, args.kde_grid_size))
        high_low_curves.append(curve_rows(parameter, f"low day {int(row['low_day'])}", values, low_w, args.kde_grid_size))
    plot_curve_grid(pd.concat(high_low_curves, ignore_index=True), args.figure_dir, "high_vs_low_information_day.png", "High- versus low-information day posteriors from broad prior")
    high_low_table = []
    for _, row in selection.iterrows():
        parameter = row["parameter"]
        values = broad_p[parameter].to_numpy(float)
        order = str(row["best_species"]).split("|")
        for label, day in [("high_information_day", int(row["best_day"])), ("low_information_day", int(row["low_day"]))]:
            weights = posterior_weights(broad_y, day, order, noise, observation_index)
            high_low_table.append({"parameter": parameter, "scenario": label, "day": day, "biomarkers": "|".join(order), **narrowing(values, weights)})
    pd.DataFrame(high_low_table).to_csv(final_tables / "bed_targeted_high_vs_low_day_narrowing.csv", index=False)
    pd.DataFrame(cumulative_rows).to_csv(final_tables / "bed_targeted_gsa_uncertainty_posterior_narrowing.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    positions = np.arange(len(stability))
    ax.bar(positions, stability["same_best_day"].astype(int))
    ax.set_xticks(positions)
    ax.set_xticklabels(stability["parameter"], rotation=45, ha="right")
    ax.set_ylabel("Same best day (1=yes)")
    ax.set_title("BED ranking stability: original admissible vs enriched bank")
    fig.savefig(args.figure_dir / "bed_ranking_stability_enriched_vs_original.png", dpi=220)
    plt.close(fig)
    copy_map = [
        ("bed_targeted_ranking_summary.png", "bed_mi_ranking_stability_12k.png"),
        ("bed_ranking_stability_enriched_vs_original.png", "bed_mi_stability_original_vs_12k.png"),
        ("posterior_independent_observation_scenarios.png", "bed_targeted_independent_observation_posteriors_3x3.png"),
        ("posterior_cumulative_biomarkers.png", "bed_targeted_cumulative_biomarker_posteriors_3x3.png"),
        ("high_vs_low_information_day.png", "bed_targeted_high_vs_low_information_day_3x3.png"),
        ("posterior_guided_parameter_specific.png", "bed_targeted_gsa_uncertainty_guided_posteriors_3x3.png"),
        ("posterior_independent_observation_scenarios.png", "bed_targeted_parameter_specific_updates_3x3.png"),
        ("bed_guided_parameter_biomarker_mi_bars.png", "bed_guided_parameter_biomarker_mi_bars.png"),
        ("bed_guided_day_mi_curves.png", "bed_guided_day_mi_curves.png"),
    ]
    for src, dst in copy_map:
        shutil.copy2(args.figure_dir / src, final_figures / dst)

    mark(11, "Saving summary")
    independent = pd.DataFrame(independent_rows)
    cumulative = pd.DataFrame(cumulative_rows)
    summary = {
        "broad_prior_sample_count": int(len(pd.read_csv(args.broad_prior_dir / "prior_parameter_samples.csv"))),
        "broad_prior_finite_rows_used_for_posteriors": int(finite_count),
        "enriched_stability_bank_sample_count": int(len(enriched_p)),
        "observable_biomarkers": observables,
        "observation_index": observation_index,
        "ranking_stability_same_best_day_fraction": float(stability["same_best_day"].mean()),
        "mean_independent_variance_reduction": float(independent["variance_reduction_fraction"].mean()),
        "mean_cumulative_final_variance_reduction": float(cumulative.sort_values("n_biomarkers").groupby("parameter").tail(1)["variance_reduction_fraction"].mean()),
        "caveat": "Broad +/-5% prior used for posterior plots; enriched ODE-confirmed bank used only for BED ranking stability/robustness.",
    }
    (args.output_dir / "targeted_bed_summary.json").write_text(json.dumps(summary, indent=2))

    best_print = selection[["parameter", "profile_class", "best_day", "top_species"]].to_dict("records")
    top_biomarkers = selection[["parameter", "top_biomarkers"]].to_dict("records")
    print(f"broad prior sample count: {summary['broad_prior_sample_count']} (finite rows used: {summary['broad_prior_finite_rows_used_for_posteriors']})")
    print(f"enriched stability bank sample count: {summary['enriched_stability_bank_sample_count']}")
    print(f"selected target parameters and profile classes: {selection[['parameter','profile_class']].to_dict('records')}")
    print(f"top biomarkers selected from global sensitivity: {top_biomarkers}")
    print(f"BED best day/species per parameter: {best_print}")
    print(f"posterior narrowing summary: mean variance reduction={summary['mean_independent_variance_reduction']:.4f}")
    print(f"cumulative biomarker narrowing summary: final mean variance reduction={summary['mean_cumulative_final_variance_reduction']:.4f}")
    print(f"ranking stability conclusion: same best day fraction={summary['ranking_stability_same_best_day_fraction']:.3f}")
    print("caveat: broad +/-5% prior used for posterior plots; enriched bank used only for BED stability/ranking robustness")


if __name__ == "__main__":
    main()
