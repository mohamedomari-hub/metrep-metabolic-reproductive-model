"""Merge existing local sensitivity, profile, and global sensitivity results."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = Path(__file__).resolve().parent
REPRESENTATIVES = [
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


def strongest_by_parameter(table: pd.DataFrame, magnitude: str) -> pd.DataFrame:
    return table.loc[table.groupby("parameter")[magnitude].idxmax()].reset_index(drop=True)


def display_class(classification: str) -> str:
    if classification == "practically identifiable":
        return "Practically identifiable"
    if classification == "boundary-limited":
        return "Boundary-limited"
    return "Weak/flat/non-identifiable"


def interpretation(row: pd.Series, local_high_threshold: float) -> str:
    local_high = row["local_sensitivity_max_abs_auc"] >= local_high_threshold
    global_high = max(row["prcc_max_abs"], row["spearman_max_abs"]) >= 0.30
    global_low = max(row["prcc_max_abs"], row["spearman_max_abs"]) <= 0.10
    cls = row["identifiability_class"]
    if cls == "Practically identifiable" and local_high and global_high:
        return "robust influential and identifiable"
    if cls == "Boundary-limited" and global_high:
        return "influential but boundary-limited"
    if cls == "Weak/flat/non-identifiable" and global_high:
        return "global effect but weakly identifiable"
    if local_high and global_low:
        return "nominal-local effect only"
    if cls == "Weak/flat/non-identifiable" and global_low:
        return "weak global signal / weak identifiability"
    return "mixed or moderate evidence"


def build_summary() -> pd.DataFrame:
    local = pd.read_csv(ROOT / "results_final/tables/metrep_non_lactating_standard_local_sensitivity_allparams.csv")
    profiles = pd.read_csv(ROOT / "results_final/tables/profile_50d_balanced_relaxed_summary.csv")
    prcc = pd.read_csv(ROOT / "analyses/global_sensitivity/outputs/global_sensitivity_prcc_auc.csv")
    spearman = pd.read_csv(ROOT / "analyses/global_sensitivity/outputs/global_sensitivity_spearman_auc.csv")

    local["local_sensitivity_max_abs_auc"] = local["relative_sensitivity"].abs()
    local_max = strongest_by_parameter(local, "local_sensitivity_max_abs_auc")
    local_max["local_sensitivity_rank"] = (
        local_max["local_sensitivity_max_abs_auc"].rank(method="min", ascending=False).astype(int)
    )
    local_max = local_max.rename(columns={"output": "local_sensitivity_main_output"})

    prcc_max = strongest_by_parameter(prcc, "abs_prcc").rename(
        columns={"abs_prcc": "prcc_max_abs", "biomarker": "prcc_main_biomarker"}
    )
    prcc_max["prcc_direction"] = np.where(prcc_max["prcc"] >= 0, "positive", "negative")
    spearman_max = strongest_by_parameter(spearman, "abs_spearman_rho").rename(
        columns={"abs_spearman_rho": "spearman_max_abs", "biomarker": "spearman_main_biomarker"}
    )

    summary = pd.DataFrame({"parameter": REPRESENTATIVES})
    summary = summary.merge(
        profiles[["parameter", "classification"]], on="parameter", how="left"
    ).merge(
        local_max[
            [
                "parameter",
                "local_sensitivity_rank",
                "local_sensitivity_max_abs_auc",
                "local_sensitivity_main_output",
            ]
        ],
        on="parameter",
        how="left",
    ).merge(
        prcc_max[["parameter", "prcc_max_abs", "prcc_main_biomarker", "prcc_direction"]],
        on="parameter",
        how="left",
    ).merge(
        spearman_max[["parameter", "spearman_max_abs", "spearman_main_biomarker"]],
        on="parameter",
        how="left",
    )
    summary["identifiability_class"] = summary["classification"].map(display_class)
    local_high_threshold = float(local_max["local_sensitivity_max_abs_auc"].quantile(0.75))
    summary["interpretation_tag"] = summary.apply(
        interpretation, axis=1, local_high_threshold=local_high_threshold
    )
    return summary[
        [
            "parameter",
            "identifiability_class",
            "local_sensitivity_rank",
            "local_sensitivity_max_abs_auc",
            "local_sensitivity_main_output",
            "prcc_max_abs",
            "prcc_main_biomarker",
            "prcc_direction",
            "spearman_max_abs",
            "spearman_main_biomarker",
            "interpretation_tag",
        ]
    ]


def save_markdown(summary: pd.DataFrame) -> None:
    display = summary.copy()
    for column in ["local_sensitivity_max_abs_auc", "prcc_max_abs", "spearman_max_abs"]:
        display[column] = display[column].map(lambda value: f"{value:.3f}")
    note = (
        "# Combined Parameter Summary\n\n"
        "Local sensitivity measures nominal one-at-a-time AUC response. PRCC/Spearman summarize "
        "parameter-biomarker AUC associations across biologically admissible simulations. "
        "Identifiability class comes from representative profile likelihood behavior. These metrics "
        "are complementary and should not be interpreted as identical sensitivity measures.\n\n"
        "Interpretation tags use a high-global-association threshold of 0.30, a low-global-association "
        "threshold of 0.10, and the upper quartile of all-parameter local sensitivity as the high-local threshold.\n\n"
    )
    (OUTPUT_DIR / "combined_parameter_summary.md").write_text(note + display.to_markdown(index=False) + "\n")


def save_figure(summary: pd.DataFrame) -> None:
    metric_columns = ["local_sensitivity_max_abs_auc", "prcc_max_abs", "spearman_max_abs"]
    values = summary[metric_columns].to_numpy(dtype=float)
    normalized = values / np.maximum(values.max(axis=0, keepdims=True), 1e-12)
    row_labels = [
        f"{row.parameter}\n{row.identifiability_class}"
        for row in summary[["parameter", "identifiability_class"]].itertuples(index=False)
    ]
    column_labels = ["Local sensitivity\nmax |AUC response|", "PRCC\nmax |association|", "Spearman\nmax |association|"]
    fig, ax = plt.subplots(figsize=(10.5, 7.5), constrained_layout=True)
    image = ax.imshow(normalized, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(3), column_labels)
    ax.set_yticks(range(len(row_labels)), row_labels, fontsize=8)
    ax.set_title("Combined parameter diagnostics: local sensitivity, identifiability, and global association")
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            color = "white" if normalized[row, column] > 0.58 else "black"
            ax.text(column, row, f"{values[row, column]:.3f}", ha="center", va="center", color=color, fontsize=8)
    fig.colorbar(image, ax=ax, label="Column-normalized strength")
    fig.savefig(OUTPUT_DIR / "combined_parameter_summary_heatmap.png", dpi=220)
    plt.close(fig)


def main() -> None:
    summary = build_summary()
    summary.to_csv(OUTPUT_DIR / "combined_parameter_summary.csv", index=False)
    save_markdown(summary)
    save_figure(summary)
    print(f"Saved combined model diagnostics to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
