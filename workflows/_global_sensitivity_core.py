"""Global sensitivity screening from existing MetRep simulation banks.

No ODE simulations are run. The unfiltered ordinary Monte Carlo bank is used
for full-prior variance-based screening, and its biologically admissible subset
is used for Spearman and PRCC screening.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BANK = PROJECT_ROOT / "local_data/input_tables_large"
DEFAULT_OBSERVABLE_ORDER = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]
REPRESENTATIVE_CLASSES = {
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
    parser = argparse.ArgumentParser(description="Run global sensitivity screening without ODE simulations.")
    parser.add_argument("--bank-dir", type=Path, default=DEFAULT_BANK)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "local_outputs/global_sensitivity")
    parser.add_argument("--variance-bins", type=int, default=20)
    parser.add_argument("--top-per-biomarker", type=int, default=10)
    parser.add_argument(
        "--biomarkers",
        nargs="+",
        help="Observable biomarkers to analyze. Omit to use all stored observable *_day_* outputs.",
    )
    return parser.parse_args()


def stored_days(outputs: pd.DataFrame, biomarker: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(biomarker)}_day_(\d+)$")
    return sorted(int(match.group(1)) for column in outputs.columns if (match := pattern.match(column)))


def infer_observable_biomarkers(outputs: pd.DataFrame) -> list[str]:
    pattern = re.compile(r"^(.+)_day_(\d+)$")
    found = sorted({match.group(1) for column in outputs.columns if (match := pattern.match(column))})
    ordered = [name for name in DEFAULT_OBSERVABLE_ORDER if name in found]
    ordered.extend(name for name in found if name not in ordered)
    return ordered


def auc_endpoints(outputs: pd.DataFrame, biomarkers: list[str]) -> pd.DataFrame:
    auc = {}
    for biomarker in biomarkers:
        days = stored_days(outputs, biomarker)
        if len(days) < 2:
            raise ValueError(f"Need at least two stored days for {biomarker}.")
        values = outputs[[f"{biomarker}_day_{day}" for day in days]].to_numpy(dtype=float)
        auc[biomarker] = np.trapezoid(values, x=np.asarray(days, dtype=float), axis=1)
    return pd.DataFrame(auc)


def variance_screening_score(x: np.ndarray, y: np.ndarray, bins: int) -> float:
    """Estimate Var(E[Y|X])/Var(Y) using equal-count bins."""

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    total_variance = float(np.var(y))
    if len(y) < bins or total_variance <= 0:
        return np.nan
    groups = pd.qcut(rankdata(x, method="average"), q=min(bins, len(y)), labels=False, duplicates="drop")
    frame = pd.DataFrame({"group": groups, "y": y})
    grouped = frame.groupby("group", observed=True)["y"].agg(["mean", "size"])
    weighted_mean = float(np.average(grouped["mean"], weights=grouped["size"]))
    explained = float(np.average((grouped["mean"] - weighted_mean) ** 2, weights=grouped["size"]))
    return explained / total_variance


def full_prior_screening(parameters: pd.DataFrame, auc: pd.DataFrame, bins: int, biomarkers: list[str]) -> pd.DataFrame:
    rows = []
    for biomarker in biomarkers:
        y = auc[biomarker].to_numpy(dtype=float)
        for parameter in parameters.columns:
            score = variance_screening_score(parameters[parameter].to_numpy(dtype=float), y, bins)
            rows.append(
                {
                    "biomarker": biomarker,
                    "parameter": parameter,
                    "variance_screening_score": score,
                    "method": "equal-count-bin estimate of Var(E[AUC|parameter])/Var(AUC)",
                    "sampling_design": "ordinary independent uniform Monte Carlo; not Sobol/Saltelli",
                    "samples": len(parameters),
                }
            )
    return pd.DataFrame(rows)


def prcc_coefficients(parameters: pd.DataFrame, y: np.ndarray) -> np.ndarray:
    ranked_x = np.column_stack([rankdata(parameters[column]) for column in parameters.columns])
    ranked_y = rankdata(y)[:, None]
    ranked = np.column_stack([ranked_x, ranked_y]).astype(float)
    ranked -= ranked.mean(axis=0)
    std = ranked.std(axis=0, ddof=1)
    ranked /= np.where(std > 0, std, 1.0)
    precision = np.linalg.pinv(np.corrcoef(ranked, rowvar=False))
    return -precision[:-1, -1] / np.sqrt(np.maximum(precision[:-1, :-1].diagonal() * precision[-1, -1], 1e-15))


def admissible_screening(parameters: pd.DataFrame, auc: pd.DataFrame, biomarkers: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    spearman_rows = []
    prcc_rows = []
    for biomarker in biomarkers:
        y = auc[biomarker].to_numpy(dtype=float)
        prcc = prcc_coefficients(parameters, y)
        for index, parameter in enumerate(parameters.columns):
            rho, p_value = spearmanr(parameters[parameter].to_numpy(dtype=float), y)
            spearman_rows.append(
                {
                    "biomarker": biomarker,
                    "parameter": parameter,
                    "spearman_rho": rho,
                    "abs_spearman_rho": abs(rho),
                    "p_value": p_value,
                    "samples": len(parameters),
                }
            )
            prcc_rows.append(
                {
                    "biomarker": biomarker,
                    "parameter": parameter,
                    "prcc": prcc[index],
                    "abs_prcc": abs(prcc[index]),
                    "samples": len(parameters),
                }
            )
    return pd.DataFrame(spearman_rows), pd.DataFrame(prcc_rows)


def save_figures(output_dir: Path, prcc: pd.DataFrame, representatives: pd.DataFrame, top_n: int) -> None:
    top_parameters = (
        prcc.groupby("parameter", as_index=False)["abs_prcc"].max()
        .nlargest(top_n, "abs_prcc")["parameter"]
        .tolist()
    )
    heatmap = prcc[prcc["parameter"].isin(top_parameters)].pivot(index="parameter", columns="biomarker", values="prcc")
    heatmap = heatmap.loc[heatmap.abs().max(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(11, max(6, 0.32 * len(heatmap))), constrained_layout=True)
    image = ax.imshow(heatmap, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(heatmap.columns)), heatmap.columns)
    ax.set_yticks(range(len(heatmap.index)), heatmap.index, fontsize=8)
    ax.set_title("Admissible-bank PRCC sensitivity: strongest parameters by biomarker AUC")
    fig.colorbar(image, ax=ax, label="PRCC")
    fig.savefig(output_dir / "global_sensitivity_heatmap_top_parameters.png", dpi=220)
    plt.close(fig)

    rep_heatmap = representatives.pivot(index="parameter", columns="biomarker", values="prcc")
    ordered = list(REPRESENTATIVE_CLASSES)
    rep_heatmap = rep_heatmap.reindex(ordered)
    fig, ax = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    image = ax.imshow(rep_heatmap, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(rep_heatmap.columns)), rep_heatmap.columns)
    labels = [f"{name}\n{REPRESENTATIVE_CLASSES[name]}" for name in rep_heatmap.index]
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_title("Global sensitivity of representative identifiability parameters")
    fig.colorbar(image, ax=ax, label="PRCC with biomarker AUC")
    fig.savefig(output_dir / "global_sensitivity_representative_identifiability_parameters.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    parameters = pd.read_csv(args.bank_dir / "prior_parameter_samples.csv")
    outputs = pd.read_csv(args.bank_dir / "ode_output_features.csv")
    admissibility = pd.read_csv(args.bank_dir / "admissibility.csv")
    if not (len(parameters) == len(outputs) == len(admissibility)):
        raise ValueError("Stored bank tables have inconsistent row counts.")

    biomarkers = args.biomarkers or infer_observable_biomarkers(outputs)
    if not biomarkers:
        raise ValueError("No observable biomarker day columns were found.")

    finite = parameters.notna().all(axis=1) & outputs.notna().all(axis=1)
    full_parameters = parameters.loc[finite].reset_index(drop=True)
    full_auc = auc_endpoints(outputs.loc[finite].reset_index(drop=True), biomarkers)
    full = full_prior_screening(full_parameters, full_auc, args.variance_bins, biomarkers)
    full.to_csv(args.output_dir / "full_prior_variance_screening.csv", index=False)
    full.sort_values(["biomarker", "variance_screening_score"], ascending=[True, False]).groupby("biomarker").head(
        args.top_per_biomarker
    ).to_csv(args.output_dir / "full_prior_top_parameters_by_biomarker.csv", index=False)

    keep = finite & admissibility["admissible"].astype(bool)
    accepted_parameters = parameters.loc[keep].reset_index(drop=True)
    accepted_auc = auc_endpoints(outputs.loc[keep].reset_index(drop=True), biomarkers)
    spearman, prcc = admissible_screening(accepted_parameters, accepted_auc, biomarkers)
    spearman.to_csv(args.output_dir / "global_sensitivity_spearman_auc.csv", index=False)
    prcc.to_csv(args.output_dir / "global_sensitivity_prcc_auc.csv", index=False)

    top = prcc.sort_values(["biomarker", "abs_prcc"], ascending=[True, False]).groupby("biomarker").head(
        args.top_per_biomarker
    )
    top.to_csv(args.output_dir / "top_parameters_by_biomarker_auc.csv", index=False)
    representatives = prcc[prcc["parameter"].isin(REPRESENTATIVE_CLASSES)].copy()
    representatives["identifiability_class"] = representatives["parameter"].map(REPRESENTATIVE_CLASSES)
    representatives = representatives.merge(
        spearman[["biomarker", "parameter", "spearman_rho", "abs_spearman_rho"]],
        on=["biomarker", "parameter"],
    )
    representatives.to_csv(args.output_dir / "representative_identifiability_parameter_summary.csv", index=False)
    save_figures(args.output_dir, prcc, representatives, args.top_per_biomarker * 2)

    note = (
        "Local sensitivity used +1% one-at-a-time perturbations around the nominal model and AUC endpoints. "
        "Global sensitivity screening uses simulation-bank ensembles and AUC endpoints. The admissible-bank "
        "Spearman/PRCC analysis summarizes parameter-biomarker associations across biologically plausible "
        "simulations. PRCC values indicate the direction and strength of monotonic association, not formal "
        "Sobol variance decomposition. The full-prior bank is ordinary random Monte Carlo, so its variance-based result is a "
        "screening analysis, not strict Sobol indices. Perturbation scales differ because each method asks a "
        "different question. BED is handled separately.\n"
        "Glucagon was excluded from the biological admissibility filter but retained as an observable "
        "biomarker for downstream uncertainty propagation, global sensitivity, and Bayesian experimental design.\n"
    )
    (args.output_dir / "README.md").write_text(note)
    print(f"Saved global sensitivity outputs to {args.output_dir}")
    print(f"Observable biomarkers used: {', '.join(biomarkers)}")
    print(f"Full-prior rows: {len(full_parameters)}; admissible rows: {len(accepted_parameters)}; ODE simulations run: 0")


if __name__ == "__main__":
    main()
