"""Surrogate-assisted Bayesian experimental design for precomputed MetRep runs.

This script does not run the ODE model. It expects tables produced from ODE
simulations, trains a validated surrogate, and then uses the surrogate outputs
to estimate posterior narrowing and mutual information for candidate
measurements.

Required inputs:
    --parameters-csv
        Prior samples. Rows are samples; columns are uncertain parameters.
    --outputs-csv
        ODE output features for the same rows. Columns are measurable outputs,
        usually species-at-time features such as FSH_day_65 or P4_day_81.
    --target-column
        Parameter column used as the BED target W.

Optional inputs:
    --candidate-map-csv
        Table with columns candidate, columns. The columns field is a
        pipe-separated list of output feature names for that candidate design.
        If omitted, every output column is treated as one scalar candidate.
    --nominal-output-csv
        One-row table containing the nominal/reference model prediction for the
        output columns. If omitted, the mean ODE output is used only as a
        fallback for synthetic posterior plots.
    --admissibility-csv
        Table with one Boolean column named admissible. Rows match samples.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import digamma
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KDTree, KernelDensity, NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import QuantileTransformer, StandardScaler


@dataclass(frozen=True)
class CandidateDesign:
    name: str
    columns: list[str]


@dataclass
class SurrogateModel:
    x_columns: list[str]
    y_columns: list[str]
    x_pipeline: Pipeline
    y_scaler: StandardScaler
    y_pca: PCA

    def predict(self, parameters: pd.DataFrame) -> pd.DataFrame:
        scores = self.x_pipeline.predict(parameters[self.x_columns])
        y_scaled = self.y_pca.inverse_transform(scores)
        y = self.y_scaler.inverse_transform(y_scaled)
        return pd.DataFrame(y, columns=self.y_columns, index=parameters.index)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train and validate a surrogate-assisted BED workflow from precomputed ODE samples."
    )
    parser.add_argument("--parameters-csv", type=Path, required=True)
    parser.add_argument("--outputs-csv", type=Path, required=True)
    parser.add_argument("--target-column", required=True)
    parser.add_argument("--candidate-map-csv", type=Path)
    parser.add_argument("--nominal-output-csv", type=Path)
    parser.add_argument("--admissibility-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("surrogate_bed_outputs"))
    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument("--validation-size", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--regressor", choices=["extra_trees", "random_forest"], default="extra_trees")
    parser.add_argument("--n-estimators", type=int, default=800)
    parser.add_argument("--min-samples-leaf", type=int, default=2)
    parser.add_argument("--pca-variance", type=float, default=0.995)
    parser.add_argument("--relative-noise", type=float, default=0.05)
    parser.add_argument("--noise-floor", type=float, default=1e-8)
    parser.add_argument(
        "--mi-noisy-measurements",
        action="store_true",
        help="Estimate MI after adding Gaussian measurement noise to predicted outputs, closer to the MATLAB BED workflow.",
    )
    parser.add_argument("--synthetic-observation-seed", type=int, default=7)
    parser.add_argument("--posterior-candidate", help="Candidate name used for posterior plots. Defaults to top MI.")
    parser.add_argument("--posterior-grid-size", type=int, default=400)
    parser.add_argument("--kde-bandwidth", type=float, default=0.05)
    parser.add_argument("--mi-neighbors", type=int, default=5)
    parser.add_argument("--mi-repeats", type=int, default=10)
    parser.add_argument(
        "--convergence-sizes",
        nargs="+",
        type=int,
        default=[250, 500, 1000, 2000, 5000, 10000],
    )
    parser.add_argument(
        "--learning-curve-sizes",
        nargs="+",
        type=float,
        default=[0.20, 0.40, 0.60, 0.80, 1.00],
        help="Fractions of the training set used for surrogate learning-curve checks.",
    )
    parser.add_argument("--max-candidates", type=int, help="Optional limit for quick diagnostic runs.")
    parser.add_argument("--make-plots", action="store_true")
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Input table not found: {path}. Replace the example path with a real CSV file "
            "containing precomputed ODE samples."
        )
    table = pd.read_csv(path)
    if table.empty:
        raise ValueError(f"Empty table: {path}")
    return table


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    parameters = read_table(args.parameters_csv)
    outputs = read_table(args.outputs_csv)
    if len(parameters) != len(outputs):
        raise ValueError("Parameter and output tables must have the same number of rows.")
    if args.target_column not in parameters.columns:
        raise ValueError(f"Target column must be in parameters table: {args.target_column}")

    admissible = np.ones(len(parameters), dtype=bool)
    if args.admissibility_csv:
        admissibility = read_table(args.admissibility_csv)
        if "admissible" not in admissibility.columns:
            raise ValueError("Admissibility table must contain an 'admissible' column.")
        if len(admissibility) != len(parameters):
            raise ValueError("Admissibility table must match the sample count.")
        admissible = admissibility["admissible"].astype(bool).to_numpy()

    parameters = parameters.loc[admissible].reset_index(drop=True)
    outputs = outputs.loc[admissible].reset_index(drop=True)
    if len(parameters) < 100:
        raise ValueError("Fewer than 100 admissible samples remain. Add more ODE samples before BED.")
    return parameters, outputs, admissible


def load_candidate_designs(args: argparse.Namespace, output_columns: list[str]) -> list[CandidateDesign]:
    if args.candidate_map_csv is None:
        designs = [CandidateDesign(name=column, columns=[column]) for column in output_columns]
    else:
        table = read_table(args.candidate_map_csv)
        required = {"candidate", "columns"}
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"Candidate map is missing columns: {sorted(missing)}")
        designs = []
        for _, row in table.iterrows():
            columns = [item.strip() for item in str(row["columns"]).split("|") if item.strip()]
            unknown = sorted(set(columns).difference(output_columns))
            if unknown:
                raise ValueError(f"Candidate {row['candidate']} uses unknown output columns: {unknown}")
            designs.append(CandidateDesign(name=str(row["candidate"]), columns=columns))

    if args.max_candidates:
        designs = designs[: args.max_candidates]
    if not designs:
        raise ValueError("No candidate designs were defined.")
    return designs


def make_regressor(args: argparse.Namespace):
    common = {
        "n_estimators": args.n_estimators,
        "random_state": args.seed,
        "min_samples_leaf": args.min_samples_leaf,
        "n_jobs": -1,
    }
    if args.regressor == "random_forest":
        return RandomForestRegressor(**common)
    return ExtraTreesRegressor(**common)


def train_surrogate(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    args: argparse.Namespace,
) -> SurrogateModel:
    y_scaler = StandardScaler()
    y_train_scaled = y_scaler.fit_transform(y_train)
    y_pca = PCA(n_components=args.pca_variance, svd_solver="full")
    y_scores = y_pca.fit_transform(y_train_scaled)

    x_pipeline = Pipeline(
        steps=[
            ("xscale", StandardScaler()),
            ("regressor", make_regressor(args)),
        ]
    )
    x_pipeline.fit(x_train, y_scores)
    return SurrogateModel(
        x_columns=list(x_train.columns),
        y_columns=list(y_train.columns),
        x_pipeline=x_pipeline,
        y_scaler=y_scaler,
        y_pca=y_pca,
    )


def split_data(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x_work, x_test, y_work, y_test = train_test_split(
        parameters,
        outputs,
        test_size=args.test_size,
        random_state=args.seed,
    )
    validation_fraction_of_work = args.validation_size / (1.0 - args.test_size)
    x_train, x_valid, y_train, y_valid = train_test_split(
        x_work,
        y_work,
        test_size=validation_fraction_of_work,
        random_state=args.seed + 1,
    )
    return x_train, x_valid, x_test, y_train, y_valid, y_test


def evaluate_predictions(observed: pd.DataFrame, predicted: pd.DataFrame, split: str) -> pd.DataFrame:
    rows = []
    for column in observed.columns:
        obs = observed[column].to_numpy(dtype=float)
        pred = predicted[column].to_numpy(dtype=float)
        rmse = float(np.sqrt(mean_squared_error(obs, pred)))
        scale = float(np.nanmax(obs) - np.nanmin(obs))
        rows.append(
            {
                "split": split,
                "output": column,
                "mae": float(mean_absolute_error(obs, pred)),
                "rmse": rmse,
                "nrmse_range": rmse / scale if scale > 0 else np.nan,
                "r2": float(r2_score(obs, pred)) if len(np.unique(obs)) > 1 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def learning_curve(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_valid: pd.DataFrame,
    y_valid: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    rng = np.random.default_rng(args.seed)
    rows = []
    n_train = len(x_train)
    for fraction in args.learning_curve_sizes:
        fraction = float(fraction)
        n_subset = max(50, min(n_train, int(round(fraction * n_train))))
        indices = rng.choice(n_train, size=n_subset, replace=False)
        model = train_surrogate(x_train.iloc[indices], y_train.iloc[indices], args)
        pred = model.predict(x_valid)
        metrics = evaluate_predictions(y_valid, pred, split="validation")
        rows.append(
            {
                "train_fraction": fraction,
                "n_train": n_subset,
                "median_r2": float(metrics["r2"].median()),
                "median_nrmse_range": float(metrics["nrmse_range"].median()),
                "max_nrmse_range": float(metrics["nrmse_range"].max()),
            }
        )
    return pd.DataFrame(rows)


def load_nominal_observation(args: argparse.Namespace, output_columns: list[str], fallback: pd.DataFrame) -> np.ndarray:
    if args.nominal_output_csv:
        nominal = read_table(args.nominal_output_csv)
        missing = sorted(set(output_columns).difference(nominal.columns))
        if missing:
            raise ValueError(f"Nominal output table is missing output columns: {missing}")
        return nominal.loc[nominal.index[0], output_columns].to_numpy(dtype=float)
    return fallback[output_columns].mean(axis=0).to_numpy(dtype=float)


def synthetic_observation(nominal: np.ndarray, relative_noise: float, noise_floor: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    sigma = np.maximum(relative_noise * np.abs(nominal), noise_floor)
    observation = nominal + sigma * rng.normal(size=len(nominal))
    return observation, sigma


def posterior_weights(predicted: pd.DataFrame, observation: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    residual = (predicted.to_numpy(dtype=float) - observation[None, :]) / sigma[None, :]
    log_likelihood = -0.5 * np.sum(residual * residual, axis=1)
    shifted = log_likelihood - float(np.max(log_likelihood))
    weights = np.exp(shifted)
    total = float(np.sum(weights))
    if total == 0 or not np.isfinite(total):
        raise ValueError("Posterior weights collapsed. Increase noise level or check candidate outputs.")
    return weights / total


def effective_sample_size(weights: np.ndarray) -> float:
    return float(1.0 / np.sum(np.square(weights)))


def weighted_kde(values: np.ndarray, weights: np.ndarray, bandwidth: float, grid_size: int) -> pd.DataFrame:
    values = np.asarray(values, dtype=float)
    lo = float(np.min(values))
    hi = float(np.max(values))
    pad = 0.05 * (hi - lo) if hi > lo else 1.0
    grid = np.linspace(lo - pad, hi + pad, grid_size)

    prior_kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    prior_kde.fit(values[:, None])
    posterior_kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    posterior_kde.fit(values[:, None], sample_weight=weights)

    prior_density = np.exp(prior_kde.score_samples(grid[:, None]))
    posterior_density = np.exp(posterior_kde.score_samples(grid[:, None]))
    prior_area = np.trapezoid(prior_density, grid)
    posterior_area = np.trapezoid(posterior_density, grid)
    if prior_area > 0:
        prior_density = prior_density / prior_area
    if posterior_area > 0:
        posterior_density = posterior_density / posterior_area

    return pd.DataFrame(
        {
            "target_grid": grid,
            "prior_density": prior_density,
            "posterior_density": posterior_density,
        }
    )


def rank_gaussianize(values: np.ndarray, seed: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    transformer = QuantileTransformer(
        n_quantiles=min(1000, len(values)),
        output_distribution="normal",
        random_state=seed,
        subsample=None,
    )
    return transformer.fit_transform(values)


def ksg_mutual_information(target: np.ndarray, measurement: np.ndarray, k: int, seed: int) -> float:
    """Kraskov-style k-nearest-neighbor MI estimate for scalar W and vector Z.

    The estimator uses the Chebyshev metric. Rank-Gaussianization is applied
    before neighbor counting to reduce scale effects between species.
    """

    target = np.asarray(target, dtype=float).reshape(-1, 1)
    measurement = np.asarray(measurement, dtype=float)
    if measurement.ndim == 1:
        measurement = measurement.reshape(-1, 1)
    n = len(target)
    if n <= k + 2:
        return np.nan

    w = rank_gaussianize(target, seed=seed)
    z = rank_gaussianize(measurement, seed=seed + 1)
    wz = np.column_stack([w, z])

    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="chebyshev")
    nbrs.fit(wz)
    distances, _ = nbrs.kneighbors(wz)
    eps = np.nextafter(distances[:, k], 0.0)

    tree_w = KDTree(w, metric="chebyshev")
    tree_z = KDTree(z, metric="chebyshev")
    n_w = np.array([len(item) - 1 for item in tree_w.query_radius(w, r=eps)], dtype=float)
    n_z = np.array([len(item) - 1 for item in tree_z.query_radius(z, r=eps)], dtype=float)

    mi = digamma(k) + digamma(n) - np.mean(digamma(n_w + 1.0) + digamma(n_z + 1.0))
    return float(max(mi, 0.0))


def estimate_mi_for_designs(
    target: np.ndarray,
    predicted_outputs: pd.DataFrame,
    designs: list[CandidateDesign],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(args.seed)
    max_n = len(target)
    summary_rows = []
    convergence_rows = []

    for design in designs:
        full_measurement = predicted_outputs[design.columns].to_numpy(dtype=float)
        if args.mi_noisy_measurements:
            sigma = np.maximum(args.relative_noise * np.abs(np.mean(full_measurement, axis=0)), args.noise_floor)
            full_measurement = full_measurement + sigma[None, :] * rng.normal(size=full_measurement.shape)
        full_mi = ksg_mutual_information(target, full_measurement, args.mi_neighbors, args.seed)
        summary_rows.append(
            {
                "candidate": design.name,
                "columns": "|".join(design.columns),
                "n_columns": len(design.columns),
                "n_samples": max_n,
                "mi_nats": full_mi,
            }
        )

        for requested_n in args.convergence_sizes:
            n = min(int(requested_n), max_n)
            if n <= args.mi_neighbors + 2:
                continue
            mi_values = []
            for repeat in range(args.mi_repeats):
                indices = rng.choice(max_n, size=n, replace=False)
                mi_values.append(
                    ksg_mutual_information(
                        target[indices],
                        full_measurement[indices],
                        args.mi_neighbors,
                        args.seed + 1000 * repeat + n,
                    )
                )
            convergence_rows.append(
                {
                    "candidate": design.name,
                    "n_samples": n,
                    "mi_mean_nats": float(np.nanmean(mi_values)),
                    "mi_sd_nats": float(np.nanstd(mi_values, ddof=1)) if len(mi_values) > 1 else 0.0,
                    "n_repeats": len(mi_values),
                }
            )

    mi_summary = pd.DataFrame(summary_rows).sort_values("mi_nats", ascending=False).reset_index(drop=True)
    mi_convergence = pd.DataFrame(convergence_rows)
    return mi_summary, mi_convergence


def posterior_for_candidate(
    candidate: CandidateDesign,
    target_values: np.ndarray,
    predicted_outputs: pd.DataFrame,
    nominal_outputs: np.ndarray,
    output_columns: list[str],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    column_index = [output_columns.index(column) for column in candidate.columns]
    nominal = nominal_outputs[column_index]
    observation, sigma = synthetic_observation(
        nominal,
        relative_noise=args.relative_noise,
        noise_floor=args.noise_floor,
        seed=args.synthetic_observation_seed,
    )
    predicted_candidate = predicted_outputs[candidate.columns]
    weights = posterior_weights(predicted_candidate, observation, sigma)
    kde = weighted_kde(target_values, weights, args.kde_bandwidth, args.posterior_grid_size)
    diagnostic = pd.DataFrame(
        {
            "candidate": [candidate.name],
            "effective_sample_size": [effective_sample_size(weights)],
            "n_samples": [len(weights)],
            "ess_fraction": [effective_sample_size(weights) / len(weights)],
            "relative_noise": [args.relative_noise],
        }
    )
    return kde, diagnostic


def save_plots(
    output_dir: Path,
    validation_metrics: pd.DataFrame,
    learning: pd.DataFrame,
    mi_summary: pd.DataFrame,
    mi_convergence: pd.DataFrame,
    posterior: pd.DataFrame,
) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    top_validation = validation_metrics[validation_metrics["split"] == "test"].copy()
    top_validation = top_validation.sort_values("nrmse_range", ascending=False).head(25)
    plt.figure(figsize=(9, 7))
    plt.barh(top_validation["output"], top_validation["nrmse_range"])
    plt.xlabel("Test NRMSE divided by output range")
    plt.ylabel("Output")
    plt.title("Surrogate validation error")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(figure_dir / "surrogate_validation_error.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(learning["n_train"], learning["median_nrmse_range"], marker="o", label="median")
    plt.plot(learning["n_train"], learning["max_nrmse_range"], marker="o", label="max")
    plt.xlabel("Training samples")
    plt.ylabel("Validation NRMSE divided by output range")
    plt.title("Surrogate learning curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_dir / "surrogate_learning_curve.png", dpi=200)
    plt.close()

    top_mi = mi_summary.head(20).copy()
    plt.figure(figsize=(9, 7))
    plt.barh(top_mi["candidate"], top_mi["mi_nats"])
    plt.xlabel("Mutual information estimate (nats)")
    plt.ylabel("Candidate design")
    plt.title("BED candidate ranking")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(figure_dir / "mi_candidate_ranking.png", dpi=200)
    plt.close()

    top_candidates = mi_summary.head(5)["candidate"].tolist()
    plt.figure(figsize=(8, 5))
    for candidate in top_candidates:
        subset = mi_convergence[mi_convergence["candidate"] == candidate].sort_values("n_samples")
        if subset.empty:
            continue
        plt.errorbar(
            subset["n_samples"],
            subset["mi_mean_nats"],
            yerr=subset["mi_sd_nats"],
            marker="o",
            capsize=3,
            label=candidate,
        )
    plt.xlabel("Monte Carlo samples")
    plt.ylabel("Mutual information estimate (nats)")
    plt.title("MI convergence for top candidates")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "mi_convergence.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(posterior["target_grid"], posterior["prior_density"], "--", label="Prior")
    plt.plot(posterior["target_grid"], posterior["posterior_density"], label="Posterior")
    plt.xlabel("Target parameter")
    plt.ylabel("Density")
    plt.title("Posterior narrowing for selected candidate")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_dir / "posterior_comparison.png", dpi=200)
    plt.close()


def write_run_metadata(
    output_dir: Path,
    args: argparse.Namespace,
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    designs: list[CandidateDesign],
    surrogate: SurrogateModel,
) -> None:
    metadata = {
        "n_samples_after_admissibility": len(parameters),
        "n_parameters": parameters.shape[1],
        "n_outputs": outputs.shape[1],
        "n_candidate_designs": len(designs),
        "target_column": args.target_column,
        "regressor": args.regressor,
        "n_estimators": args.n_estimators,
        "min_samples_leaf": args.min_samples_leaf,
        "pca_variance_requested": args.pca_variance,
        "pca_components_used": int(surrogate.y_pca.n_components_),
        "pca_variance_explained": float(np.sum(surrogate.y_pca.explained_variance_ratio_)),
        "relative_noise": args.relative_noise,
        "mi_neighbors": args.mi_neighbors,
        "mi_repeats": args.mi_repeats,
        "seed": args.seed,
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    parameters, outputs, _ = load_inputs(args)
    designs = load_candidate_designs(args, list(outputs.columns))
    x_train, x_valid, x_test, y_train, y_valid, y_test = split_data(parameters, outputs, args)

    learning = learning_curve(x_train, y_train, x_valid, y_valid, args)
    learning.to_csv(args.output_dir / "surrogate_learning_curve.csv", index=False)

    surrogate = train_surrogate(x_train, y_train, args)
    valid_pred = surrogate.predict(x_valid)
    test_pred = surrogate.predict(x_test)
    validation_metrics = pd.concat(
        [
            evaluate_predictions(y_valid, valid_pred, split="validation"),
            evaluate_predictions(y_test, test_pred, split="test"),
        ],
        ignore_index=True,
    )
    validation_metrics.to_csv(args.output_dir / "surrogate_validation_metrics.csv", index=False)

    final_surrogate = train_surrogate(parameters, outputs, args)
    predicted_all = final_surrogate.predict(parameters)
    predicted_all.to_csv(args.output_dir / "surrogate_predicted_outputs.csv", index=False)

    target_values = parameters[args.target_column].to_numpy(dtype=float)
    mi_summary, mi_convergence = estimate_mi_for_designs(target_values, predicted_all, designs, args)
    mi_summary.to_csv(args.output_dir / "mi_candidate_ranking.csv", index=False)
    mi_convergence.to_csv(args.output_dir / "mi_convergence.csv", index=False)

    selected_name = args.posterior_candidate or str(mi_summary.loc[0, "candidate"])
    selected = next((candidate for candidate in designs if candidate.name == selected_name), None)
    if selected is None:
        raise ValueError(f"Posterior candidate was not found: {selected_name}")

    nominal_outputs = load_nominal_observation(args, list(outputs.columns), outputs)
    posterior, posterior_diagnostic = posterior_for_candidate(
        selected,
        target_values,
        predicted_all,
        nominal_outputs,
        list(outputs.columns),
        args,
    )
    posterior.to_csv(args.output_dir / "posterior_prior_comparison.csv", index=False)
    posterior_diagnostic.to_csv(args.output_dir / "posterior_diagnostics.csv", index=False)

    write_run_metadata(args.output_dir, args, parameters, outputs, designs, final_surrogate)
    if args.make_plots:
        save_plots(args.output_dir, validation_metrics, learning, mi_summary, mi_convergence, posterior)

    print(f"Saved surrogate BED outputs to {args.output_dir}")
    print("Main files:")
    print("  surrogate_validation_metrics.csv")
    print("  surrogate_learning_curve.csv")
    print("  mi_candidate_ranking.csv")
    print("  mi_convergence.csv")
    print("  posterior_prior_comparison.csv")
    print("  posterior_diagnostics.csv")
    print("Review validation, MI convergence, and posterior ESS before comparing with thesis figures.")


if __name__ == "__main__":
    main()
