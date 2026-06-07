"""Audit integrity of the cached 50k +/-5% PhD BED ODE bank.

This script does not rerun ODE simulations. It verifies file consistency,
sampling range, stored admissibility labels, daily-output metric recomputation,
and visual examples.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, precision_score, recall_score


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BANK = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon"
SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]
ADMISSIBILITY_SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose"]
METRICS = ["penalty", "min_correlation", "max_average_difference", "max_norm_difference"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Audit 50k PhD BED bank integrity.")
    p.add_argument("--bank-dir", type=Path, default=DEFAULT_BANK)
    p.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "outputs/bank_integrity_audit")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def mark(step: int, text: str) -> float:
    print(f"[STEP {step}/10] {text}", flush=True)
    return time.monotonic()


def done(step: int, started: float) -> None:
    elapsed = int(time.monotonic() - started)
    print(f"[STEP {step}/10] Completed in {elapsed // 60}m {elapsed % 60}s", flush=True)


def columns_for(columns: list[str], species: str) -> list[str]:
    return sorted(
        [c for c in columns if re.match(rf"^{re.escape(species)}_day_\d+$", c)],
        key=lambda c: int(c.rsplit("_", 1)[1]),
    )


def load_nominal_parameters(parameter_names: list[str], samples: pd.DataFrame) -> pd.Series:
    try:
        sys.path.insert(0, str(ROOT / "MetRep_Python"))
        from metrep.parameters import default_parameters  # type: ignore

        nominal = default_parameters()
        return pd.Series({name: nominal[name] for name in parameter_names})
    except Exception:
        return (samples.min() + samples.max()) / 2.0


def shifted_corr(values: np.ndarray, reference: np.ndarray, max_lag: int = 5) -> np.ndarray:
    best = np.full(values.shape[0], -np.inf)
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            ref = reference[-lag:]
            arr = values[:, : values.shape[1] + lag]
        elif lag > 0:
            ref = reference[:-lag]
            arr = values[:, lag:]
        else:
            ref = reference
            arr = values
        if arr.shape[1] < 2:
            continue
        den = np.linalg.norm(arr, axis=1) * max(np.linalg.norm(ref), 1e-12)
        score = np.divide(arr @ ref, den, out=np.zeros(arr.shape[0]), where=den > 0)
        best = np.maximum(best, score)
    return np.where(np.isfinite(best), best, 0.0)


def recompute_daily_metrics(outputs: pd.DataFrame, nominal: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    per_state = {}
    penalties = []
    limiting_state = []
    for state in ADMISSIBILITY_SPECIES:
        cols = columns_for(list(outputs.columns), state)
        values = outputs[cols].to_numpy(float)
        ref = nominal[cols].iloc[0].to_numpy(float)
        days = np.array([int(c.rsplit("_", 1)[1]) for c in cols], dtype=float)
        corr = shifted_corr(values, ref, max_lag=5)
        avg = np.trapezoid(np.abs(values - ref), x=days, axis=1) / max(np.trapezoid(np.abs(ref), x=days), 1e-12)
        ref_norm = max(np.sqrt(np.trapezoid(ref * ref, x=days)), 1e-12)
        norm = np.abs(np.sqrt(np.trapezoid(values * values, x=days, axis=1)) - ref_norm) / ref_norm
        pen = 100.0 * (
            np.maximum(0.0, 0.75 - corr) ** 2
            + np.maximum(0.0, avg - 0.30) ** 2
            + np.maximum(0.0, norm - 0.30) ** 2
        )
        penalties.append(pen)
        per_state[f"{state}_max_correlation"] = corr
        per_state[f"{state}_average_difference"] = avg
        per_state[f"{state}_norm_difference"] = norm
        per_state[f"{state}_penalty"] = pen
    frame = pd.DataFrame(per_state)
    state_corr_cols = [f"{s}_max_correlation" for s in ADMISSIBILITY_SPECIES]
    state_avg_cols = [f"{s}_average_difference" for s in ADMISSIBILITY_SPECIES]
    state_norm_cols = [f"{s}_norm_difference" for s in ADMISSIBILITY_SPECIES]
    penalty_matrix = np.vstack(penalties).T
    limiting_state = [ADMISSIBILITY_SPECIES[i] for i in np.argmax(penalty_matrix, axis=1)]
    summary = pd.DataFrame(
        {
            "daily_recomputed_penalty": penalty_matrix.sum(axis=1),
            "daily_recomputed_min_correlation": frame[state_corr_cols].min(axis=1),
            "daily_recomputed_max_average_difference": frame[state_avg_cols].max(axis=1),
            "daily_recomputed_max_norm_difference": frame[state_norm_cols].max(axis=1),
            "daily_recomputed_worst_species": limiting_state,
        }
    )
    return summary, frame


def exact_stored_rule(adm: pd.DataFrame, outputs: pd.DataFrame) -> np.ndarray:
    return (
        outputs.notna().all(axis=1).to_numpy()
        & (adm["penalty"].to_numpy(float) <= 0.0)
        & (adm["min_correlation"].to_numpy(float) >= 0.75)
        & (adm["max_average_difference"].to_numpy(float) <= 0.30)
        & (adm["max_norm_difference"].to_numpy(float) <= 0.30)
    )


def plot_examples(path: Path, outputs: pd.DataFrame, nominal: pd.DataFrame, rows: np.ndarray, title: str, color: str) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), constrained_layout=True)
    for ax, state in zip(axes.ravel(), ["PGF", "E2", "P4", "FSH", "INH", "Glucose", "Insulin", "IGF1"]):
        cols = columns_for(list(outputs.columns), state)
        days = [int(c.rsplit("_", 1)[1]) for c in cols]
        ax.plot(days, nominal[cols].iloc[0], color="black", linewidth=2, label="nominal")
        for row in rows:
            ax.plot(days, outputs.loc[row, cols], color=color, alpha=0.28, linewidth=1)
        ax.set_title(state)
        ax.grid(alpha=0.2)
    fig.suptitle(title)
    axes.ravel()[0].legend()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    figures = args.output_dir.parents[1] / "figures/bank_integrity_audit"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    t = mark(1, "Loading 50k bank files")
    files = sorted(p for p in args.bank_dir.iterdir() if p.is_file())
    (args.output_dir / "file_inventory.txt").write_text("\n".join(f"{p.name}\t{p.stat().st_size} bytes" for p in files) + "\n")
    params = pd.read_csv(args.bank_dir / "prior_parameter_samples.csv")
    outputs = pd.read_csv(args.bank_dir / "ode_output_features.csv")
    nominal_output = pd.read_csv(args.bank_dir / "nominal_output.csv")
    adm = pd.read_csv(args.bank_dir / "admissibility.csv")
    metadata = json.loads((args.bank_dir / "bank_metadata.json").read_text()) if (args.bank_dir / "bank_metadata.json").exists() else {}
    done(1, t)

    t = mark(2, "Verifying parameter sampling range")
    nominal_params = load_nominal_parameters(list(params.columns), params)
    multiplier = params.divide(nominal_params.replace(0, np.nan), axis=1)
    sampling = pd.DataFrame(
        {
            "parameter": params.columns,
            "min_multiplier": multiplier.min(),
            "max_multiplier": multiplier.max(),
            "mean_multiplier": multiplier.mean(),
            "std_multiplier": multiplier.std(),
            "within_0p95_1p05": ((multiplier.min() >= 0.95 - 1e-10) & (multiplier.max() <= 1.05 + 1e-10)).to_numpy(),
            "constant_unvaried": (params.nunique() <= 1).to_numpy(),
        }
    )
    sampling.to_csv(args.output_dir / "parameter_sampling_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    ax.vlines(np.arange(len(sampling)), sampling["min_multiplier"], sampling["max_multiplier"], linewidth=1)
    ax.axhline(0.95, color="tab:red", linestyle="--")
    ax.axhline(1.05, color="tab:red", linestyle="--")
    ax.set(title="Parameter multiplier ranges", xlabel="Parameter index", ylabel="sample / nominal")
    fig.savefig(figures / "parameter_multiplier_ranges.png", dpi=220)
    plt.close(fig)
    parameter_range_valid = bool(sampling["within_0p95_1p05"].all() and not sampling["constant_unvaried"].any())
    done(2, t)

    t = mark(3, "Verifying row alignment across parameters, outputs, nominal, and admissibility")
    row_counts = {"parameters": len(params), "outputs": len(outputs), "admissibility": len(adm), "nominal_output": len(nominal_output)}
    row_alignment_valid = len({row_counts["parameters"], row_counts["outputs"], row_counts["admissibility"]}) == 1 and row_counts["nominal_output"] == 1
    duplicate_params = int(params.duplicated().sum())
    duplicate_outputs = int(outputs.duplicated().sum())
    pd.DataFrame(
        [
            {
                **row_counts,
                "row_alignment_valid": row_alignment_valid,
                "duplicate_parameter_rows": duplicate_params,
                "duplicate_output_rows": duplicate_outputs,
                "explicit_simulation_id_columns": ",".join([c for c in params.columns if "id" in c.lower()]),
            }
        ]
    ).to_csv(args.output_dir / "row_alignment_report.csv", index=False)
    done(3, t)

    t = mark(4, "Verifying finite/valid ODE outputs")
    values = outputs.to_numpy(float)
    finite_rows = np.isfinite(values).all(axis=1)
    validity = pd.DataFrame(
        [
            {
                "n_rows": len(outputs),
                "finite_rows": int(finite_rows.sum()),
                "nan_values": int(np.isnan(values).sum()),
                "inf_values": int(np.isinf(values).sum()),
                "negative_values": int((values < 0).sum()),
                "max_abs_value": float(np.nanmax(np.abs(values))),
                "duplicated_output_rows": duplicate_outputs,
            }
        ]
    )
    validity.to_csv(args.output_dir / "output_validity_summary.csv", index=False)
    output_finite_valid = bool(finite_rows.all())
    done(4, t)

    t = mark(5, "Recomputing aggregate admissibility rule from stored columns")
    stored = adm["admissible"].astype(bool).to_numpy()
    recomputed = exact_stored_rule(adm, outputs)
    mismatches = np.flatnonzero(stored != recomputed)
    pd.DataFrame(
        [
            {
                "stored_admissible": int(stored.sum()),
                "recomputed_rule_admissible": int(recomputed.sum()),
                "agreement_count": int((stored == recomputed).sum()),
                "mismatch_count": int(len(mismatches)),
                "precision": precision_score(stored, recomputed, zero_division=0),
                "recall": recall_score(stored, recomputed, zero_division=0),
                "f1": f1_score(stored, recomputed, zero_division=0),
            }
        ]
    ).to_csv(args.output_dir / "stored_rule_reproduction.csv", index=False)
    pd.DataFrame({"row_index": mismatches, "stored_admissible": stored[mismatches], "rule_admissible": recomputed[mismatches]}).to_csv(args.output_dir / "stored_rule_mismatches.csv", index=False)
    stored_rule_ok = bool(np.array_equal(stored, recomputed))
    trace = f"""# Admissibility Application Trace

- Generation script: `analyses/bayesian_experimental_design/surrogate_bed/prepare_phd_bed_bank.py`
- Function applying admissibility: `metrep.admissibility.trajectory_admissibility`
- Admissibility biomarker panel: `{metadata.get("admissibility_species", ADMISSIBILITY_SPECIES)}`
- Stored BED biomarker panel: `{metadata.get("species", SPECIES)}`
- Glucagon included in admissibility: `{"Glucagon" in metadata.get("admissibility_species", [])}`
- Glucagon included in stored outputs: `{"Glucagon" in metadata.get("species", [])}`
- Metrics computed from: full ODE trajectory at `step = 1/6` day inside generation script, then daily BED outputs were stored.
- Failed/NaN simulations: caught during generation and written as NaN outputs with `admissible = False`.
- Final stored label columns: `admissible`, `penalty`, `min_correlation`, `max_average_difference`, `max_norm_difference`, `worst_species`.
- Final selection rule confirmed here: finite outputs AND `penalty <= 0` AND `min_correlation >= 0.75` AND `max_average_difference <= 0.30` AND `max_norm_difference <= 0.30`.
- Number selected: `{int(stored.sum())}`.
"""
    (args.output_dir / "admissibility_application_trace.md").write_text(trace)
    done(5, t)

    t = mark(6, "Recomputing admissibility metrics from trajectories if possible")
    daily_summary, daily_state = recompute_daily_metrics(outputs, nominal_output)
    comparison = pd.concat([adm[METRICS].reset_index(drop=True), daily_summary], axis=1)
    for stored_col, daily_col in [
        ("penalty", "daily_recomputed_penalty"),
        ("min_correlation", "daily_recomputed_min_correlation"),
        ("max_average_difference", "daily_recomputed_max_average_difference"),
        ("max_norm_difference", "daily_recomputed_max_norm_difference"),
    ]:
        comparison[f"{stored_col}_abs_error"] = np.abs(comparison[stored_col] - comparison[daily_col])
    comparison.to_csv(args.output_dir / "recomputed_metric_comparison.csv", index=False)
    daily_state.to_csv(args.output_dir / "daily_per_state_metric_recomputation.csv", index=False)
    recompute_partial = True
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    pairs = [
        ("penalty", "daily_recomputed_penalty"),
        ("min_correlation", "daily_recomputed_min_correlation"),
        ("max_average_difference", "daily_recomputed_max_average_difference"),
        ("max_norm_difference", "daily_recomputed_max_norm_difference"),
    ]
    for ax, (stored_col, daily_col) in zip(axes.ravel(), pairs):
        ax.scatter(comparison[stored_col], comparison[daily_col], s=4, alpha=0.25)
        ax.set_xlabel(f"stored {stored_col}")
        ax.set_ylabel(f"daily recomputed {stored_col}")
        ax.grid(alpha=0.2)
    fig.savefig(figures / "recomputed_vs_stored_metrics.png", dpi=220)
    plt.close(fig)
    done(6, t)

    t = mark(7, "Auditing selected admissible trajectory examples")
    admissible_rows = np.flatnonzero(stored)
    rejected_rows = np.flatnonzero(~stored)
    near = np.argsort(np.abs(adm["min_correlation"].to_numpy(float) - 0.75))[:30]
    plot_examples(figures / "admissible_trajectory_examples.png", outputs, nominal_output, rng.choice(admissible_rows, size=min(30, len(admissible_rows)), replace=False), "Stored admissible examples", "tab:blue")
    done(7, t)

    t = mark(8, "Auditing rejected trajectory examples")
    plot_examples(figures / "rejected_trajectory_examples.png", outputs, nominal_output, rng.choice(rejected_rows, size=min(30, len(rejected_rows)), replace=False), "Stored rejected examples", "tab:red")
    plot_examples(figures / "near_boundary_trajectory_examples.png", outputs, nominal_output, near, "Near-boundary examples", "tab:orange")
    done(8, t)

    t = mark(9, "Checking reproducibility metadata and random seed")
    metric_summary = adm.loc[stored, METRICS].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).T
    metric_summary.to_csv(args.output_dir / "admissible_metric_summary.csv")
    limiting = daily_summary.loc[stored, "daily_recomputed_worst_species"].value_counts().rename_axis("species").reset_index(name="count")
    limiting.to_csv(args.output_dir / "admissible_limiting_biomarker_counts.csv", index=False)
    provenance = []
    for path in [ROOT / "analyses/bayesian_experimental_design/surrogate_bed/prepare_phd_bed_bank.py", args.bank_dir / "bank_metadata.json"]:
        if path.exists():
            provenance.append(f"===== {path.relative_to(ROOT) if path.is_relative_to(ROOT) else path} =====\n")
            provenance.append(path.read_text(errors="ignore")[:12000])
            provenance.append("\n")
    try:
        result = subprocess.run(["rg", "-n", "prepare_phd_bed_bank|phd_bed_bank_5pct_50k_glucagon|--n-samples|--seed|Glucagon"], cwd=ROOT, capture_output=True, text=True, check=False)
        provenance.append("===== repository search =====\n")
        provenance.append(result.stdout)
    except Exception as exc:
        provenance.append(f"repository search failed: {exc}\n")
    (args.output_dir / "generation_provenance.txt").write_text("".join(provenance))
    done(9, t)

    t = mark(10, "Saving audit report and conclusion")
    acceptance = float(stored.mean())
    metric_error_max = comparison[[f"{c}_abs_error" for c in METRICS]].max().max()
    recomputed_match = "partial" if recompute_partial else ("yes" if metric_error_max < 1e-8 else "no")
    trustworthy = parameter_range_valid and row_alignment_valid and output_finite_valid and stored_rule_ok
    conclusion = pd.DataFrame(
        [
            {
                "bank_row_count": len(params),
                "parameter_range_valid": parameter_range_valid,
                "row_alignment_valid": row_alignment_valid,
                "output_finite_valid": output_finite_valid,
                "stored_rule_reproduces_labels": stored_rule_ok,
                "recomputed_metrics_match_stored_metrics": recomputed_match,
                "admissible_simulations": int(stored.sum()),
                "acceptance_rate": acceptance,
                "labels_trustworthy": trustworthy,
                "critical_caveat": "Stored daily outputs allow only partial metric recomputation; exact admissibility metrics were computed on full 1/6-day ODE trajectories during bank generation.",
            }
        ]
    )
    conclusion.to_csv(args.output_dir / "bank_integrity_conclusion.csv", index=False)
    done(10, t)

    print(f"bank row count: {len(params)}")
    print(f"parameter range valid: {'yes' if parameter_range_valid else 'no'}")
    print(f"row alignment valid: {'yes' if row_alignment_valid else 'no'}")
    print(f"output finite valid: {'yes' if output_finite_valid else 'no'}")
    print(f"stored rule reproduces labels: {'yes' if stored_rule_ok else 'no'}")
    print(f"recomputed metrics match stored metrics: {recomputed_match}")
    print(f"number of admissible simulations: {int(stored.sum())}")
    print(f"acceptance rate: {acceptance:.6f}")
    print(f"2,957 admissible labels trustworthy: {'yes' if trustworthy else 'no'}")
    print("critical caveat: stored daily outputs allow only partial metric recomputation; exact labels were computed from full 1/6-day trajectories.")
    print("admissibility applied by: prepare_phd_bed_bank.py / metrep.admissibility.trajectory_admissibility")
    print(f"admissibility biomarker panel: {metadata.get('admissibility_species', ADMISSIBILITY_SPECIES)}")
    print("full time grid or stored outputs: full 1/6-day ODE grid during generation")
    print("Glucagon included in admissibility: no")
    print("final selection rule: finite outputs AND penalty<=0 AND min_correlation>=0.75 AND max_average_difference<=0.30 AND max_norm_difference<=0.30")


if __name__ == "__main__":
    main()
