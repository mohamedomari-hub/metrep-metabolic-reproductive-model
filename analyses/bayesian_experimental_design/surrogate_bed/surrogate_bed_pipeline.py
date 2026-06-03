"""Surrogate-assisted BED scaffold.

This script is intentionally a scaffold: it trains and validates a surrogate
from precomputed ODE samples, then computes posterior weights and a Monte Carlo
mutual-information estimate. It does not run the MetRep ODE model.

Expected inputs:
    --parameters-csv: rows are prior samples, columns are uncertain parameters
    --outputs-csv: rows match parameter samples, columns are output features
    --target-column: parameter/output column used as the BED target W

Optional:
    --admissibility-csv: rows match samples and contain an `admissible` column
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KernelDensity
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parameters-csv", type=Path, required=True)
    parser.add_argument("--outputs-csv", type=Path, required=True)
    parser.add_argument("--target-column", required=True)
    parser.add_argument("--admissibility-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("surrogate_bed_outputs"))
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pca-components", type=float, default=0.99)
    parser.add_argument("--n-estimators", type=int, default=600)
    parser.add_argument("--relative-noise", type=float, default=0.05)
    parser.add_argument("--kde-bandwidth", type=float, default=0.05)
    parser.add_argument("--convergence-sizes", nargs="+", type=int, default=[250, 500, 1000, 2000, 5000])
    return parser.parse_args()


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    parameters = pd.read_csv(args.parameters_csv)
    outputs = pd.read_csv(args.outputs_csv)
    if len(parameters) != len(outputs):
        raise ValueError("Parameter and output tables must have the same number of rows.")

    admissible = np.ones(len(parameters), dtype=bool)
    if args.admissibility_csv:
        admissibility = pd.read_csv(args.admissibility_csv)
        if "admissible" not in admissibility.columns:
            raise ValueError("Admissibility table must contain an 'admissible' column.")
        if len(admissibility) != len(parameters):
            raise ValueError("Admissibility table must match the sample count.")
        admissible = admissibility["admissible"].astype(bool).to_numpy()

    return parameters.loc[admissible].reset_index(drop=True), outputs.loc[admissible].reset_index(drop=True), admissible


def train_surrogate(
    parameters: pd.DataFrame,
    outputs: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[Pipeline, pd.DataFrame]:
    x_train, x_test, y_train, y_test = train_test_split(
        parameters,
        outputs,
        test_size=args.test_size,
        random_state=args.seed,
    )

    model = Pipeline(
        steps=[
            ("xscale", StandardScaler()),
            ("regressor", ExtraTreesRegressor(
                n_estimators=args.n_estimators,
                random_state=args.seed,
                min_samples_leaf=2,
                n_jobs=-1,
            )),
        ]
    )

    # PCA is applied to outputs manually because scikit-learn's Pipeline is
    # input-oriented. The surrogate predicts PCA scores, then reconstructs.
    y_scaler = StandardScaler()
    y_train_scaled = y_scaler.fit_transform(y_train)
    pca = PCA(n_components=args.pca_components, svd_solver="full", random_state=args.seed)
    y_train_scores = pca.fit_transform(y_train_scaled)

    model.fit(x_train, y_train_scores)
    pred_scores = model.predict(x_test)
    pred_outputs = y_scaler.inverse_transform(pca.inverse_transform(pred_scores))
    pred_outputs = pd.DataFrame(pred_outputs, columns=outputs.columns, index=y_test.index)

    metrics = []
    for column in outputs.columns:
        metrics.append(
            {
                "output": column,
                "mae": mean_absolute_error(y_test[column], pred_outputs[column]),
                "r2": r2_score(y_test[column], pred_outputs[column]),
            }
        )

    model.output_scaler_ = y_scaler
    model.output_pca_ = pca
    return model, pd.DataFrame(metrics)


def predict_outputs(model: Pipeline, parameters: pd.DataFrame) -> pd.DataFrame:
    scores = model.predict(parameters)
    outputs = model.output_scaler_.inverse_transform(model.output_pca_.inverse_transform(scores))
    return pd.DataFrame(outputs)


def posterior_weights(predicted: pd.DataFrame, observation: np.ndarray, relative_noise: float) -> np.ndarray:
    sigma = np.maximum(relative_noise * np.abs(observation), 1e-8)
    residual = (predicted.to_numpy(dtype=float) - observation[None, :]) / sigma[None, :]
    log_likelihood = -0.5 * np.sum(residual * residual, axis=1)
    shifted = log_likelihood - float(np.max(log_likelihood))
    weights = np.exp(shifted)
    return weights / float(np.sum(weights))


def weighted_kde(values: np.ndarray, weights: np.ndarray, bandwidth: float, grid: np.ndarray) -> np.ndarray:
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    kde.fit(values[:, None], sample_weight=weights)
    density = np.exp(kde.score_samples(grid[:, None]))
    area = np.trapz(density, grid)
    return density / area if area > 0 else density


def mutual_information_knn(target: np.ndarray, measurement: np.ndarray, bandwidth: float) -> float:
    """Simple KDE plug-in MI estimator for a scalar target and vector measurement."""

    target = np.asarray(target, dtype=float).reshape(-1, 1)
    measurement = np.asarray(measurement, dtype=float)
    joint = np.column_stack([target, measurement])

    kde_w = KernelDensity(bandwidth=bandwidth).fit(target)
    kde_z = KernelDensity(bandwidth=bandwidth).fit(measurement)
    kde_wz = KernelDensity(bandwidth=bandwidth).fit(joint)

    log_w = kde_w.score_samples(target)
    log_z = kde_z.score_samples(measurement)
    log_wz = kde_wz.score_samples(joint)
    return float(np.mean(log_wz - log_w - log_z))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    parameters, outputs, _ = load_inputs(args)
    if args.target_column not in parameters.columns and args.target_column not in outputs.columns:
        raise ValueError(f"Target column not found: {args.target_column}")

    model, metrics = train_surrogate(parameters, outputs, args)
    metrics.to_csv(args.output_dir / "surrogate_validation_metrics.csv", index=False)

    predicted = predict_outputs(model, parameters)
    predicted.columns = outputs.columns
    nominal_observation = outputs.mean(axis=0).to_numpy(dtype=float)
    weights = posterior_weights(predicted, nominal_observation, args.relative_noise)

    target_values = (
        parameters[args.target_column].to_numpy(dtype=float)
        if args.target_column in parameters.columns
        else outputs[args.target_column].to_numpy(dtype=float)
    )
    grid = np.linspace(float(np.min(target_values)), float(np.max(target_values)), 300)
    posterior = weighted_kde(target_values, weights, args.kde_bandwidth, grid)
    pd.DataFrame({"target_grid": grid, "posterior_density": posterior}).to_csv(
        args.output_dir / "posterior_weighted_kde.csv",
        index=False,
    )

    convergence = []
    rng = np.random.default_rng(args.seed)
    max_n = len(target_values)
    for requested in args.convergence_sizes:
        n = min(int(requested), max_n)
        if n < 20:
            continue
        index = rng.choice(max_n, size=n, replace=False)
        mi = mutual_information_knn(
            target_values[index],
            predicted.iloc[index].to_numpy(dtype=float),
            args.kde_bandwidth,
        )
        convergence.append({"n_samples": n, "mi_estimate": mi})
    pd.DataFrame(convergence).to_csv(args.output_dir / "mi_convergence.csv", index=False)

    print(f"Saved surrogate BED scaffold outputs to {args.output_dir}")
    print("Review validation metrics and MI convergence before reporting any surrogate BED result.")


if __name__ == "__main__":
    main()
