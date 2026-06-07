"""Adaptive SMC enrichment of ODE-confirmed admissible parameter sets.

The machine-learning model is used only to rank SMC proposals before ODE
simulation. Exact admissibility is always recomputed from the full 1/6-day ODE
trajectory, using the same aggregate rule as the 50k PhD BED bank.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "MetRep_Python"
sys.path.insert(0, str(PYTHON_ROOT))

from model_definition.admissibility import AdmissibilityThresholds, trajectory_admissibility
from model_definition.initial_conditions import STATE_INDEX
from model_definition.parameters import PARAMETERS, default_parameters
from model_definition.scenarios import constant_non_lactating
from model_definition.simulate import run_simulation


DEFAULT_BANK = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "outputs/smc_ml_enrichment"
SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose", "Glucagon"]
ADMISSIBILITY_SPECIES = ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose"]
DAYS = list(range(54, 90))
PRIOR_HALF_RANGE = 0.05
AGG_METRICS = ["penalty", "min_correlation", "max_average_difference", "max_norm_difference"]
STATE_METRICS = [
    "PGF_max_correlation",
    "E2_max_correlation",
    "PGF_average_difference",
    "E2_average_difference",
    "PGF_norm_difference",
    "E2_norm_difference",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Adaptive SMC + ML pre-screening for admissible MetRep parameters.")
    parser.add_argument("--bank-dir", type=Path, default=DEFAULT_BANK)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-rounds", type=int, default=5)
    parser.add_argument("--proposals-per-round", type=int, default=20_000)
    parser.add_argument("--ode-audit-per-round", type=int, default=2_000)
    parser.add_argument("--initial-noise", type=float, default=0.01)
    parser.add_argument("--min-noise", type=float, default=0.0025)
    parser.add_argument("--max-noise", type=float, default=0.03)
    parser.add_argument("--increase-threshold", type=float, default=0.25)
    parser.add_argument("--decrease-threshold", type=float, default=0.10)
    parser.add_argument("--increase-factor", type=float, default=1.20)
    parser.add_argument("--decrease-factor", type=float, default=0.70)
    parser.add_argument("--trees", type=int, default=500)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--resume-existing", action="store_true", help="Continue from an existing SMC output directory.")
    parser.add_argument("--start-round", type=int, default=None, help="First round number for continuation; inferred from existing summary if omitted.")
    parser.add_argument("--min-rms-distance", type=float, default=0.0, help="Minimum RMS distance from accepted particles in normalized prior coordinates.")
    parser.add_argument("--cluster-resampling", action="store_true", help="Sample seed particles approximately evenly across source rounds.")
    parser.add_argument("--stop-below-acceptance", type=float, default=0.10)
    parser.add_argument("--stop-on-low-acceptance", action="store_true")
    parser.add_argument("--method", default="BDF")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def mark(step: int, text: str) -> float:
    print(f"[STEP {step}/10] {text}", flush=True)
    return time.monotonic()


def done(step: int, started: float) -> None:
    elapsed = int(time.monotonic() - started)
    print(f"[STEP {step}/10] Completed in {elapsed // 60}m {elapsed % 60}s", flush=True)


def exact_rule(admissibility: pd.DataFrame, outputs: pd.DataFrame) -> np.ndarray:
    finite = outputs.notna().all(axis=1).to_numpy()
    return (
        finite
        & (admissibility["penalty"].to_numpy(float) <= 0.0)
        & (admissibility["min_correlation"].to_numpy(float) >= 0.75)
        & (admissibility["max_average_difference"].to_numpy(float) <= 0.30)
        & (admissibility["max_norm_difference"].to_numpy(float) <= 0.30)
    )


def normalize_parameters(frame: pd.DataFrame, columns: list[str], nominal_vector: pd.Series) -> np.ndarray:
    values = frame[columns].to_numpy(float)
    nominal_values = nominal_vector.reindex(columns).to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = (values / nominal_values - 1.0) / PRIOR_HALF_RANGE
    return np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)


def rms_distance_to_pool(candidates: pd.DataFrame, pool: pd.DataFrame, columns: list[str], nominal_vector: pd.Series) -> np.ndarray:
    if candidates.empty or pool.empty:
        return np.full(len(candidates), np.inf)
    x = normalize_parameters(candidates, columns, nominal_vector)
    ref = normalize_parameters(pool, columns, nominal_vector)
    nn = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(ref)
    dist, _ = nn.kneighbors(x)
    return dist[:, 0] / np.sqrt(len(columns))


def sample_particle_base(
    particle_pool: pd.DataFrame,
    parameter_columns: list[str],
    count: int,
    rng: np.random.Generator,
    cluster_resampling: bool,
) -> np.ndarray:
    if not cluster_resampling or "source_round" not in particle_pool.columns:
        return particle_pool.sample(
            n=count,
            replace=True,
            random_state=int(rng.integers(0, 2**31 - 1)),
        )[parameter_columns].to_numpy(float)
    groups = list(particle_pool.groupby("source_round"))
    per_group = int(np.ceil(count / max(len(groups), 1)))
    parts = []
    for _, group in groups:
        parts.append(
            group.sample(
                n=per_group,
                replace=True,
                random_state=int(rng.integers(0, 2**31 - 1)),
            )[parameter_columns]
        )
    sampled = pd.concat(parts, ignore_index=True).sample(
        n=count,
        replace=False,
        random_state=int(rng.integers(0, 2**31 - 1)),
    )
    return sampled.to_numpy(float)


def columns_for(columns: list[str], state: str) -> list[str]:
    prefix = f"{state}_day_"
    selected = [col for col in columns if col.startswith(prefix)]
    return sorted(selected, key=lambda value: int(value.rsplit("_", 1)[1]))


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
        numerator = arr @ ref
        denominator = np.linalg.norm(arr, axis=1) * max(np.linalg.norm(ref), 1e-12)
        score = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
        best = np.maximum(best, score)
    return np.where(np.isfinite(best), best, 0.0)


def per_state_targets(outputs: pd.DataFrame, nominal: pd.DataFrame) -> pd.DataFrame:
    data: dict[str, np.ndarray] = {}
    for state in ["PGF", "E2"]:
        cols = columns_for(list(outputs.columns), state)
        values = outputs[cols].to_numpy(float)
        reference = nominal[cols].iloc[0].to_numpy(float)
        days = np.array([int(col.rsplit("_", 1)[1]) for col in cols], dtype=float)
        data[f"{state}_max_correlation"] = shifted_corr(values, reference)
        ref_auc = max(np.trapezoid(np.abs(reference), x=days), 1e-12)
        data[f"{state}_average_difference"] = np.trapezoid(np.abs(values - reference), x=days, axis=1) / ref_auc
        ref_norm = max(np.sqrt(np.trapezoid(reference * reference, x=days)), 1e-12)
        data[f"{state}_norm_difference"] = np.abs(np.sqrt(np.trapezoid(values * values, x=days, axis=1)) - ref_norm) / ref_norm
    return pd.DataFrame(data)


def boundary_score(metrics: pd.DataFrame) -> np.ndarray:
    return np.maximum.reduce(
        [
            metrics["penalty"].to_numpy(float),
            0.75 - metrics["min_correlation"].to_numpy(float),
            metrics["max_average_difference"].to_numpy(float) - 0.30,
            metrics["max_norm_difference"].to_numpy(float) - 0.30,
        ]
    )


def pgf_e2_score(metrics: pd.DataFrame) -> np.ndarray:
    margins = np.vstack(
        [
            metrics["PGF_max_correlation"].to_numpy(float) - 0.75,
            metrics["E2_max_correlation"].to_numpy(float) - 0.75,
            0.30 - metrics["PGF_average_difference"].to_numpy(float),
            0.30 - metrics["E2_average_difference"].to_numpy(float),
            0.30 - metrics["PGF_norm_difference"].to_numpy(float),
            0.30 - metrics["E2_norm_difference"].to_numpy(float),
        ]
    )
    return -np.min(margins, axis=0)


def score_frame(agg: pd.DataFrame, state: pd.DataFrame, scaler: StandardScaler | None = None) -> pd.DataFrame:
    raw = pd.DataFrame(
        {
            "aggregate_margin_score": boundary_score(agg),
            "pgf_e2_margin_score": pgf_e2_score(state),
        }
    )
    if scaler is None:
        scaled = StandardScaler().fit_transform(raw)
    else:
        scaled = scaler.transform(raw)
    raw["combined_score"] = scaled[:, 0] + scaled[:, 1]
    raw["conservative_combined_score"] = np.maximum(scaled[:, 0], scaled[:, 1])
    raw["bottleneck_score"] = np.maximum(raw["aggregate_margin_score"], raw["pgf_e2_margin_score"])
    return raw


def extract_days(simulation) -> dict[str, float]:
    return {
        f"{state}_day_{day}": float(np.interp(day, simulation.t, simulation.y[:, STATE_INDEX[state]]))
        for day in DAYS
        for state in SPECIES
    }


def audit_parameters(
    candidates: pd.DataFrame,
    round_index: int,
    output_dir: Path,
    scenario,
    reference,
    nominal: dict[str, float],
    parameter_columns: list[str],
    args: argparse.Namespace,
) -> pd.DataFrame:
    path = output_dir / f"round_{round_index:02d}_ode_audit_results.csv"
    results = pd.read_csv(path).to_dict("records") if path.exists() else []
    thresholds = AdmissibilityThresholds()
    start = len(results)
    if start:
        print(f"  -> resuming round {round_index} audit at {start}/{len(candidates)}", flush=True)
    for row_index in range(start, len(candidates)):
        row = candidates.iloc[row_index]
        params = dict(nominal)
        params.update(row[parameter_columns].to_dict())
        result = {
            "round": round_index,
            "audit_row": row_index,
            "ranking_score": float(row["ranking_score"]),
            "predicted_probability_proxy": float(row["predicted_probability_proxy"]),
        }
        for col in parameter_columns:
            result[col] = float(row[col])
        try:
            simulation = run_simulation(
                scenario,
                parameters=params,
                method=args.method,
                rtol=args.rtol,
                atol=args.atol,
            )
            output_features = extract_days(simulation)
            _, summary = trajectory_admissibility(
                reference,
                simulation,
                states=ADMISSIBILITY_SPECIES,
                thresholds=thresholds,
            )
            result.update(output_features)
            result.update({f"ode_{key}": value for key, value in summary.items()})
            finite = bool(np.isfinite(list(output_features.values())).all())
            exact = (
                finite
                and float(summary["penalty"]) <= 0.0
                and float(summary["min_correlation"]) >= 0.75
                and float(summary["max_average_difference"]) <= 0.30
                and float(summary["max_norm_difference"]) <= 0.30
            )
            result["ode_confirmed_admissible"] = bool(exact)
            result["error"] = ""
        except Exception as exc:
            result.update(
                {
                    "ode_admissible": False,
                    "ode_penalty": np.inf,
                    "ode_min_correlation": np.nan,
                    "ode_max_average_difference": np.nan,
                    "ode_max_norm_difference": np.nan,
                    "ode_worst_species": "simulation_failed",
                    "ode_confirmed_admissible": False,
                    "error": str(exc),
                }
            )
        results.append(result)
        if (row_index + 1) % args.checkpoint_every == 0 or row_index + 1 == len(candidates):
            pd.DataFrame(results).to_csv(path, index=False)
            print(f"  -> round {round_index}: ODE audit {row_index + 1}/{len(candidates)}", flush=True)
    return pd.DataFrame(results)


def plot_outputs(output_dir: Path, figures: Path) -> None:
    round_path = output_dir / "smc_round_summary.csv"
    score_path = output_dir / "smc_score_samples.csv"
    if not round_path.exists():
        return
    summary = pd.read_csv(round_path)
    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.plot(summary["round"], summary["ode_acceptance_rate"], marker="o")
    ax.set(xlabel="SMC round", ylabel="ODE acceptance rate", title="ODE-confirmed acceptance by SMC round")
    ax.grid(alpha=0.25)
    fig.savefig(figures / "acceptance_by_round.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.plot(summary["round"], summary["perturbation_scale_start"], marker="o", label="start")
    ax.plot(summary["round"], summary["perturbation_scale_next"], marker="o", label="next")
    ax.set(xlabel="SMC round", ylabel="Multiplicative noise SD", title="Adaptive perturbation scale")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(figures / "perturbation_scale_by_round.png", dpi=220)
    plt.close(fig)

    if score_path.exists():
        scores = pd.read_csv(score_path)
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        for round_index, group in scores.groupby("round"):
            ax.hist(group["ranking_score"], bins=50, alpha=0.35, density=True, label=f"round {round_index}")
        ax.set(xlabel="ML ranking score; lower is better", ylabel="Density", title="Proposal score distribution by round")
        ax.legend(fontsize=8)
        fig.savefig(figures / "score_distribution_by_round.png", dpi=220)
        plt.close(fig)

    if "median_selected_nearest_particle_rms_distance" in summary.columns:
        fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
        ax.plot(summary["round"], summary["median_selected_nearest_particle_rms_distance"], marker="o")
        ax.set(
            xlabel="SMC round",
            ylabel="Median RMS distance to accepted pool",
            title="Accepted-pool novelty of ODE-audited proposals",
        )
        ax.grid(alpha=0.25)
        fig.savefig(figures / "diversity_by_round.png", dpi=220)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir.parent.parent / "figures/smc_ml_enrichment"
    figures.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    t = mark(1, "Loading 50k bank")
    parameters = pd.read_csv(args.bank_dir / "prior_parameter_samples.csv")
    outputs = pd.read_csv(args.bank_dir / "ode_output_features.csv")
    admissibility = pd.read_csv(args.bank_dir / "admissibility.csv")
    nominal_output = pd.read_csv(args.bank_dir / "nominal_output.csv")
    metadata = json.loads((args.bank_dir / "bank_metadata.json").read_text()) if (args.bank_dir / "bank_metadata.json").exists() else {}
    labels = exact_rule(admissibility, outputs)
    stored = admissibility["admissible"].astype(bool).to_numpy()
    if not np.array_equal(labels, stored):
        raise RuntimeError("Stored admissible labels do not match the exact aggregate rule.")
    parameter_columns = list(parameters.columns)
    initial_admissible_count = int(labels.sum())
    done(1, t)

    t = mark(2, "Initializing particles from ODE-confirmed admissible set")
    nominal = default_parameters()
    nominal_vector = pd.Series({param.name: nominal[param.name] for param in PARAMETERS})
    nominal_vector = nominal_vector.reindex(parameter_columns)
    lower = (0.95 * nominal_vector).to_numpy(float)
    upper = (1.05 * nominal_vector).to_numpy(float)
    lower = np.minimum(lower, upper)
    upper = np.maximum(lower, upper)
    particle_path = args.output_dir / "smc_particle_pool.csv"
    if args.resume_existing and particle_path.exists():
        particle_pool = pd.read_csv(particle_path)
        missing = sorted(set(parameter_columns).difference(particle_pool.columns))
        if missing:
            raise RuntimeError(f"Existing particle pool is missing parameter columns: {missing[:5]}")
    else:
        particle_pool = parameters.loc[labels].copy().reset_index(drop=True)
        particle_pool.insert(0, "source_round", 0)
        particle_pool.insert(1, "source", "initial_50k_ode_confirmed")
    particle_pool.to_csv(args.output_dir / "smc_particle_pool.csv", index=False)
    done(2, t)

    t = mark(3, "Training ML score model from existing bank")
    clean = np.isfinite(admissibility[AGG_METRICS].to_numpy(float)).all(axis=1) & outputs.notna().all(axis=1).to_numpy()
    x_train = parameters.loc[clean].reset_index(drop=True)
    y_agg = admissibility.loc[clean, AGG_METRICS].reset_index(drop=True)
    y_state = per_state_targets(outputs.loc[clean].reset_index(drop=True), nominal_output)
    agg_model = ExtraTreesRegressor(
        n_estimators=args.trees,
        min_samples_leaf=2,
        n_jobs=args.n_jobs,
        random_state=args.seed,
    ).fit(x_train, y_agg)
    state_model = ExtraTreesRegressor(
        n_estimators=args.trees,
        min_samples_leaf=2,
        n_jobs=args.n_jobs,
        random_state=args.seed + 1,
    ).fit(x_train, y_state[STATE_METRICS])
    agg_fit = pd.DataFrame(agg_model.predict(x_train), columns=AGG_METRICS)
    state_fit = pd.DataFrame(state_model.predict(x_train), columns=STATE_METRICS)
    scaler = StandardScaler().fit(score_frame(agg_fit, state_fit)[["aggregate_margin_score", "pgf_e2_margin_score"]])
    done(3, t)

    t = mark(4, "Starting SMC generation loop")
    scenario = constant_non_lactating(days=100.0, step=1.0 / 6.0)
    reference = run_simulation(
        scenario,
        parameters=nominal,
        method=args.method,
        rtol=args.rtol,
        atol=args.atol,
    )
    existing_round_path = args.output_dir / "smc_round_summary.csv"
    existing_rounds = pd.read_csv(existing_round_path) if args.resume_existing and existing_round_path.exists() else pd.DataFrame()
    if args.resume_existing and not existing_rounds.empty and args.initial_noise == 0.01:
        scale = float(existing_rounds["perturbation_scale_next"].iloc[-1])
    else:
        scale = float(args.initial_noise)
    round_rows: list[dict[str, object]] = existing_rounds.to_dict("records") if not existing_rounds.empty else []
    added_round_rows: list[dict[str, object]] = []
    confirmed_path = args.output_dir / "smc_confirmed_admissible_parameters.csv"
    rejected_path = args.output_dir / "smc_rejected_audit_parameters.csv"
    score_path = args.output_dir / "smc_score_samples.csv"
    all_confirmed: list[pd.DataFrame] = [pd.read_csv(confirmed_path)] if args.resume_existing and confirmed_path.exists() else []
    all_rejected: list[pd.DataFrame] = [pd.read_csv(rejected_path)] if args.resume_existing and rejected_path.exists() else []
    score_samples: list[pd.DataFrame] = [pd.read_csv(score_path)] if args.resume_existing and score_path.exists() else []
    if args.start_round is not None:
        first_round = args.start_round
    elif args.resume_existing and not existing_rounds.empty:
        first_round = int(existing_rounds["round"].max()) + 1
    else:
        first_round = 1
    last_round = first_round + args.n_rounds - 1
    done(4, t)

    for round_index in range(first_round, last_round + 1):
        t = mark(5, "Proposing local perturbations around admissible particles")
        base = sample_particle_base(
            particle_pool,
            parameter_columns,
            args.proposals_per_round,
            rng,
            args.cluster_resampling,
        )
        proposal = base * (1.0 + rng.normal(0.0, scale, size=base.shape))
        proposal = np.clip(proposal, lower, upper)
        proposals = pd.DataFrame(proposal, columns=parameter_columns)
        proposals = proposals.drop_duplicates().reset_index(drop=True)
        done(5, t)

        t = mark(6, "ML pre-screening proposals")
        agg_pred = pd.DataFrame(agg_model.predict(proposals), columns=AGG_METRICS)
        state_pred = pd.DataFrame(state_model.predict(proposals), columns=STATE_METRICS)
        ranked_scores = score_frame(agg_pred, state_pred, scaler)
        ranked = proposals.copy()
        ranked["ranking_score"] = ranked_scores["conservative_combined_score"].to_numpy()
        ranked["predicted_probability_proxy"] = 1.0 / (1.0 + np.exp(np.clip(ranked["ranking_score"], -60, 60)))
        for col in AGG_METRICS:
            ranked[f"predicted_{col}"] = agg_pred[col].to_numpy()
        for col in STATE_METRICS:
            ranked[f"predicted_{col}"] = state_pred[col].to_numpy()
        ranked["nearest_particle_rms_distance"] = rms_distance_to_pool(ranked, particle_pool, parameter_columns, nominal_vector)
        if args.min_rms_distance > 0:
            filtered = ranked.loc[ranked["nearest_particle_rms_distance"] >= args.min_rms_distance].copy()
            if len(filtered) >= min(args.ode_audit_per_round, len(ranked)):
                ranked = filtered
        sample_scores = ranked[["ranking_score", "predicted_probability_proxy"]].sample(
            n=min(5000, len(ranked)),
            random_state=args.seed + round_index,
        )
        sample_scores.insert(0, "round", round_index)
        score_samples.append(sample_scores)
        selected = ranked.nsmallest(min(args.ode_audit_per_round, len(ranked)), "ranking_score").reset_index(drop=True)
        selected.to_csv(args.output_dir / f"round_{round_index:02d}_ml_selected_proposals.csv", index=False)
        done(6, t)

        t = mark(7, "Running ODE audit on selected proposals")
        audited = audit_parameters(
            selected,
            round_index,
            args.output_dir,
            scenario,
            reference,
            nominal,
            parameter_columns,
            args,
        )
        done(7, t)

        t = mark(8, "Applying exact admissibility rule")
        confirmed_mask = audited["ode_confirmed_admissible"].astype(bool).to_numpy()
        confirmed = audited.loc[confirmed_mask].copy()
        rejected = audited.loc[~confirmed_mask].copy()
        if not confirmed.empty:
            all_confirmed.append(confirmed)
            new_particles = confirmed[parameter_columns].copy().reset_index(drop=True)
            new_particles.insert(0, "source_round", round_index)
            new_particles.insert(1, "source", "smc_ode_confirmed")
            particle_pool = pd.concat([particle_pool, new_particles], ignore_index=True)
            particle_pool.to_csv(args.output_dir / "smc_particle_pool.csv", index=False)
        if not rejected.empty:
            all_rejected.append(rejected)
        done(8, t)

        t = mark(9, "Updating particle pool and perturbation scale")
        accepted = int(confirmed_mask.sum())
        audited_count = int(len(audited))
        acceptance_rate = accepted / max(audited_count, 1)
        median_distance = float(selected["nearest_particle_rms_distance"].median()) if "nearest_particle_rms_distance" in selected else np.nan
        near_duplicate_selected_rate = (
            float((selected["nearest_particle_rms_distance"] <= args.min_rms_distance).mean())
            if args.min_rms_distance > 0 and "nearest_particle_rms_distance" in selected
            else 0.0
        )
        next_scale = scale
        if acceptance_rate >= args.increase_threshold:
            next_scale = min(args.max_noise, scale * args.increase_factor)
        elif acceptance_rate <= args.decrease_threshold:
            next_scale = max(args.min_noise, scale * args.decrease_factor)
        round_rows.append(
            {
                "round": round_index,
                "proposals_generated": int(len(proposals)),
                "ode_audited": audited_count,
                "ode_confirmed_admissible": accepted,
                "ode_rejected": int(audited_count - accepted),
                "ode_acceptance_rate": acceptance_rate,
                "particle_pool_size": int(len(particle_pool)),
                "new_particle_pool_size": int(len(particle_pool) - initial_admissible_count),
                "perturbation_scale_start": scale,
                "perturbation_scale_next": next_scale,
                "ranking_score_min": float(selected["ranking_score"].min()) if not selected.empty else np.nan,
                "ranking_score_median": float(selected["ranking_score"].median()) if not selected.empty else np.nan,
                "median_selected_nearest_particle_rms_distance": median_distance,
                "near_duplicate_selected_rate": near_duplicate_selected_rate,
            }
        )
        added_round_rows.append(round_rows[-1])
        pd.DataFrame(round_rows).to_csv(args.output_dir / "smc_round_summary.csv", index=False)
        pd.DataFrame(added_round_rows).to_csv(args.output_dir / "smc_round_summary_continued.csv", index=False)
        pd.DataFrame(round_rows)[["round", "ode_acceptance_rate", "ode_audited", "ode_confirmed_admissible"]].to_csv(
            args.output_dir / "smc_acceptance_by_round.csv",
            index=False,
        )
        if score_samples:
            pd.concat(score_samples, ignore_index=True).to_csv(args.output_dir / "smc_score_samples.csv", index=False)
        if all_confirmed:
            pd.concat(all_confirmed, ignore_index=True).to_csv(args.output_dir / "smc_confirmed_admissible_parameters.csv", index=False)
        if all_rejected:
            pd.concat(all_rejected, ignore_index=True).to_csv(args.output_dir / "smc_rejected_audit_parameters.csv", index=False)
        scale = next_scale
        done(9, t)
        if args.stop_on_low_acceptance and acceptance_rate < args.stop_below_acceptance:
            print(f"  -> stopping: acceptance {acceptance_rate:.4f} below {args.stop_below_acceptance:.4f}", flush=True)
            break

    t = mark(10, "Saving enriched ODE-confirmed bank")
    if not all_confirmed and not (args.output_dir / "smc_confirmed_admissible_parameters.csv").exists():
        pd.DataFrame(columns=["round", "audit_row", *parameter_columns]).to_csv(
            args.output_dir / "smc_confirmed_admissible_parameters.csv",
            index=False,
        )
    if not all_rejected and not (args.output_dir / "smc_rejected_audit_parameters.csv").exists():
        pd.DataFrame(columns=["round", "audit_row", *parameter_columns]).to_csv(
            args.output_dir / "smc_rejected_audit_parameters.csv",
            index=False,
        )
    if not round_rows:
        pd.DataFrame(
            columns=[
                "round",
                "proposals_generated",
                "ode_audited",
                "ode_confirmed_admissible",
                "ode_rejected",
                "ode_acceptance_rate",
                "particle_pool_size",
                "new_particle_pool_size",
                "perturbation_scale_start",
                "perturbation_scale_next",
                "ranking_score_min",
                "ranking_score_median",
                "median_selected_nearest_particle_rms_distance",
                "near_duplicate_selected_rate",
            ]
        ).to_csv(args.output_dir / "smc_round_summary.csv", index=False)
        pd.DataFrame(columns=["round", "ode_acceptance_rate", "ode_audited", "ode_confirmed_admissible"]).to_csv(
            args.output_dir / "smc_acceptance_by_round.csv",
            index=False,
        )
    plot_outputs(args.output_dir, figures)
    summary = pd.DataFrame(round_rows)
    new_count = int(summary["ode_confirmed_admissible"].sum()) if not summary.empty else 0
    added_count = int(pd.DataFrame(added_round_rows)["ode_confirmed_admissible"].sum()) if added_round_rows else 0
    best_acceptance = float(summary["ode_acceptance_rate"].max()) if not summary.empty else 0.0
    final_summary = {
        "initial_admissible_count": initial_admissible_count,
        "new_ode_confirmed_admissible_count": new_count,
        "new_ode_confirmed_admissible_count_this_run": added_count,
        "total_enriched_admissible_bank_size": initial_admissible_count + new_count,
        "best_round_acceptance": best_acceptance,
        "final_perturbation_scale": scale,
        "bank_dir": str(args.bank_dir),
        "metadata": metadata,
        "caveat": "Only ODE-confirmed samples are scientific truth; ML scores only rank proposals before ODE audit.",
    }
    (args.output_dir / "smc_final_summary.json").write_text(json.dumps(final_summary, indent=2, default=str))
    done(10, t)

    rates = pd.DataFrame(added_round_rows)["ode_acceptance_rate"].round(4).tolist() if added_round_rows else []
    print(f"initial admissible count: {initial_admissible_count}")
    print(f"new ODE-confirmed admissible count: {added_count}")
    print(f"ODE acceptance rate per round: {rates}")
    print(f"total enriched admissible bank size: {initial_admissible_count + new_count}")
    print(f"best round acceptance: {best_acceptance:.4f}")
    print(f"final perturbation scale: {scale:.5f}")
    print("caveat: only ODE-confirmed samples are scientific truth")


if __name__ == "__main__":
    main()
