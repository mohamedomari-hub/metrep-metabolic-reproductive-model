"""Diagnose diversity of SMC-enriched ODE-confirmed admissible parameters."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "MetRep_Python"
sys.path.insert(0, str(PYTHON_ROOT))

from model_definition.parameters import PARAMETERS, default_parameters


DEFAULT_BANK = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon"
DEFAULT_SMC = Path(__file__).resolve().parent / "outputs/smc_ml_enrichment"
PRIOR_HALF_RANGE = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Diagnose SMC admissible-bank diversity.")
    parser.add_argument("--bank-dir", type=Path, default=DEFAULT_BANK)
    parser.add_argument("--smc-dir", type=Path, default=DEFAULT_SMC)
    parser.add_argument("--near-duplicate-rms-threshold", type=float, default=0.01)
    parser.add_argument("--n-clusters", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def parameter_columns(frame: pd.DataFrame) -> list[str]:
    names = [p.name for p in PARAMETERS]
    return [name for name in names if name in frame.columns]


def normalize_parameters(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    nominal = default_parameters()
    nominal_values = np.array([nominal[col] for col in columns], dtype=float)
    values = frame[columns].to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = (values / nominal_values - 1.0) / PRIOR_HALF_RANGE
    return np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)


def rms_nearest(query: np.ndarray, reference: np.ndarray) -> np.ndarray:
    if len(query) == 0 or len(reference) == 0:
        return np.full(len(query), np.nan)
    nn = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(reference)
    dist, _ = nn.kneighbors(query)
    return dist[:, 0] / np.sqrt(query.shape[1])


def sequential_previous_distances(values: np.ndarray) -> np.ndarray:
    out = np.full(len(values), np.nan)
    previous: list[np.ndarray] = []
    for i, row in enumerate(values):
        if previous:
            ref = np.vstack(previous)
            out[i] = rms_nearest(row.reshape(1, -1), ref)[0]
        previous.append(row)
    return out


def exact_duplicate_summary(original: pd.DataFrame, new: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    original_keys = set(pd.util.hash_pandas_object(original[columns], index=False).astype(str))
    new_hash = pd.util.hash_pandas_object(new[columns], index=False).astype(str)
    duplicate_with_original = new_hash.isin(original_keys).to_numpy()
    duplicate_within_new = new_hash.duplicated(keep=False).to_numpy()
    return pd.DataFrame(
        [
            {
                "new_samples": int(len(new)),
                "exact_duplicates_with_original": int(duplicate_with_original.sum()),
                "exact_duplicate_rate_with_original": float(duplicate_with_original.mean()) if len(new) else 0.0,
                "exact_duplicates_within_new": int(duplicate_within_new.sum()),
                "exact_duplicate_rate_within_new": float(duplicate_within_new.mean()) if len(new) else 0.0,
            }
        ]
    )


def coverage_summary(original: pd.DataFrame, enriched: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows = []
    for name in columns:
        for label, frame in [("original", original), ("enriched", enriched)]:
            values = frame[name].to_numpy(float)
            rows.append(
                {
                    "parameter": name,
                    "set": label,
                    "min": float(np.min(values)),
                    "q01": float(np.quantile(values, 0.01)),
                    "q05": float(np.quantile(values, 0.05)),
                    "median": float(np.quantile(values, 0.50)),
                    "q95": float(np.quantile(values, 0.95)),
                    "q99": float(np.quantile(values, 0.99)),
                    "max": float(np.max(values)),
                    "range": float(np.max(values) - np.min(values)),
                }
            )
    return pd.DataFrame(rows)


def plot_figures(figures: Path, original_z: np.ndarray, new_z: np.ndarray, dist: pd.DataFrame, coverage: pd.DataFrame) -> None:
    combined = np.vstack([original_z, new_z]) if len(new_z) else original_z
    labels = np.array(["original"] * len(original_z) + ["new_smc"] * len(new_z))
    pca = PCA(n_components=2, random_state=42).fit_transform(combined)
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for label, color, alpha in [("original", "tab:blue", 0.25), ("new_smc", "tab:orange", 0.65)]:
        mask = labels == label
        ax.scatter(pca[mask, 0], pca[mask, 1], s=10, alpha=alpha, label=label, color=color)
    ax.set(xlabel="PC1", ylabel="PC2", title="PCA projection of normalized parameter multipliers")
    ax.legend()
    fig.savefig(figures / "pca_projection_original_vs_enriched.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.hist(dist["nearest_original_rms_distance"].dropna(), bins=60, alpha=0.75)
    ax.set(xlabel="RMS distance to nearest original admissible", ylabel="New SMC samples", title="New-sample distance from original bank")
    fig.savefig(figures / "nearest_original_distance_distribution.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.hist(dist["nearest_previous_smc_rms_distance"].dropna(), bins=60, alpha=0.75, color="tab:green")
    ax.set(xlabel="RMS distance to nearest previous SMC sample", ylabel="New SMC samples", title="Sequential SMC novelty")
    fig.savefig(figures / "nearest_previous_smc_distance_distribution.png", dpi=220)
    plt.close(fig)

    pivot = coverage.pivot(index="parameter", columns="set", values="range").dropna()
    ratio = (pivot["enriched"] / pivot["original"].replace(0, np.nan)).sort_values()
    fig, ax = plt.subplots(figsize=(8, max(5, 0.16 * len(ratio))), constrained_layout=True)
    ax.barh(ratio.index, ratio.values)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1)
    ax.set(xlabel="Enriched/original parameter range", title="Parameter coverage change after SMC enrichment")
    fig.savefig(figures / "parameter_coverage_original_vs_enriched.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = args.smc_dir / "diversity_diagnostics"
    figures = args.smc_dir.parent.parent / "figures/smc_ml_enrichment/diversity_diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    print("[STEP 1/7] Loading original and SMC admissible banks", flush=True)
    parameters = pd.read_csv(args.bank_dir / "prior_parameter_samples.csv")
    admissibility = pd.read_csv(args.bank_dir / "admissibility.csv")
    original = parameters.loc[admissibility["admissible"].astype(bool)].reset_index(drop=True)
    particle_pool = pd.read_csv(args.smc_dir / "smc_particle_pool.csv")
    confirmed_path = args.smc_dir / "smc_confirmed_admissible_parameters.csv"
    confirmed = pd.read_csv(confirmed_path) if confirmed_path.exists() else pd.DataFrame()
    columns = parameter_columns(parameters)
    new = particle_pool.loc[particle_pool["source"].eq("smc_ode_confirmed"), columns].reset_index(drop=True)
    enriched = pd.concat([original[columns], new[columns]], ignore_index=True)

    print("[STEP 2/7] Checking duplicates and near-duplicates", flush=True)
    duplicate = exact_duplicate_summary(original, new, columns)
    original_z = normalize_parameters(original, columns)
    new_z = normalize_parameters(new, columns)
    enriched_z = normalize_parameters(enriched, columns)
    nearest_original = rms_nearest(new_z, original_z)
    nearest_enriched = rms_nearest(new_z, enriched_z)
    previous = sequential_previous_distances(new_z)
    distance = pd.DataFrame(
        {
            "new_index": np.arange(len(new)),
            "nearest_original_rms_distance": nearest_original,
            "nearest_any_enriched_rms_distance": nearest_enriched,
            "nearest_previous_smc_rms_distance": previous,
            "near_duplicate_original": nearest_original <= args.near_duplicate_rms_threshold,
            "near_duplicate_previous_smc": previous <= args.near_duplicate_rms_threshold,
        }
    )
    if not confirmed.empty and "round" in confirmed:
        distance["round"] = confirmed.loc[confirmed["ode_confirmed_admissible"].astype(bool), "round"].to_numpy()[: len(distance)]

    print("[STEP 3/7] Comparing parameter coverage", flush=True)
    coverage = coverage_summary(original[columns], enriched[columns], columns)

    print("[STEP 4/7] Computing PCA and cluster diversity", flush=True)
    n_clusters = min(args.n_clusters, max(2, len(enriched_z)))
    kmeans = MiniBatchKMeans(n_clusters=n_clusters, random_state=args.seed, n_init=10, batch_size=2048)
    cluster = kmeans.fit_predict(enriched_z)
    original_cluster = cluster[: len(original_z)]
    new_cluster = cluster[len(original_z) :]
    cluster_counts = pd.DataFrame(
        {
            "cluster": np.arange(n_clusters),
            "original_count": np.bincount(original_cluster, minlength=n_clusters),
            "new_smc_count": np.bincount(new_cluster, minlength=n_clusters) if len(new_cluster) else np.zeros(n_clusters, dtype=int),
        }
    )
    cluster_counts["enriched_count"] = cluster_counts["original_count"] + cluster_counts["new_smc_count"]
    occupied_original = int((cluster_counts["original_count"] > 0).sum())
    occupied_enriched = int((cluster_counts["enriched_count"] > 0).sum())

    print("[STEP 5/7] Summarizing diversity by round", flush=True)
    if "round" in distance:
        by_round = (
            distance.groupby("round")
            .agg(
                new_samples=("new_index", "size"),
                median_nearest_original_rms_distance=("nearest_original_rms_distance", "median"),
                near_duplicate_original_rate=("near_duplicate_original", "mean"),
                median_nearest_previous_smc_rms_distance=("nearest_previous_smc_rms_distance", "median"),
                near_duplicate_previous_smc_rate=("near_duplicate_previous_smc", "mean"),
            )
            .reset_index()
        )
    else:
        by_round = pd.DataFrame()

    print("[STEP 6/7] Saving tables and figures", flush=True)
    duplicate.to_csv(output_dir / "duplicate_summary.csv", index=False)
    distance.to_csv(output_dir / "nearest_neighbor_distances.csv", index=False)
    coverage.to_csv(output_dir / "parameter_coverage_original_vs_enriched.csv", index=False)
    cluster_counts.to_csv(output_dir / "cluster_counts.csv", index=False)
    by_round.to_csv(output_dir / "diversity_by_round.csv", index=False)
    plot_figures(figures, original_z, new_z, distance, coverage)

    print("[STEP 7/7] Saving concise recommendation", flush=True)
    near_original_rate = float(distance["near_duplicate_original"].mean()) if len(distance) else 0.0
    median_original_distance = float(np.nanmedian(nearest_original)) if len(nearest_original) else np.nan
    median_previous_distance = float(np.nanmedian(previous)) if len(previous) else np.nan
    diversity_acceptable = bool(near_original_rate < 0.20 and occupied_enriched >= occupied_original)
    recommendation = pd.DataFrame(
        [
            {
                "original_admissible_count": int(len(original)),
                "new_smc_admissible_count": int(len(new)),
                "enriched_admissible_count": int(len(enriched)),
                "exact_duplicate_rate_with_original": float(duplicate.iloc[0]["exact_duplicate_rate_with_original"]),
                "near_duplicate_original_rate": near_original_rate,
                "median_nearest_original_rms_distance": median_original_distance,
                "median_nearest_previous_smc_rms_distance": median_previous_distance,
                "occupied_clusters_original": occupied_original,
                "occupied_clusters_enriched": occupied_enriched,
                "diversity_acceptable": diversity_acceptable,
                "recommendation": (
                    "continue SMC with duplicate/near-duplicate filtering"
                    if diversity_acceptable
                    else "pause continuation and increase diversity controls"
                ),
            }
        ]
    )
    recommendation.to_csv(output_dir / "diversity_recommendation.csv", index=False)
    print(f"original admissible count: {len(original)}")
    print(f"current enriched count: {len(enriched)}")
    print(f"duplicate/near-duplicate rate: {duplicate.iloc[0]['exact_duplicate_rate_with_original']:.4f} / {near_original_rate:.4f}")
    print(f"diversity acceptable: {'yes' if diversity_acceptable else 'no'}")
    print(f"median nearest original RMS distance: {median_original_distance:.5f}")
    print(f"occupied clusters original/enriched: {occupied_original}/{occupied_enriched}")


if __name__ == "__main__":
    main()
