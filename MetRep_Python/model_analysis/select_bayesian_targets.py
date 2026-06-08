"""Select fixed 3x3 Bayesian inference/BED target parameters."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


TARGET_CLASSES = {
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
    parser = argparse.ArgumentParser(description="Select Bayesian inference target parameters.")
    parser.add_argument("--profile-dir", type=Path, required=True)
    parser.add_argument("--global-sensitivity-dir", type=Path, required=True)
    parser.add_argument("--uncertainty-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = args.profile_dir / "combined_parameter_summary.csv"
    combined = pd.read_csv(combined_path) if combined_path.exists() else pd.DataFrame({"parameter": list(TARGET_CLASSES)})
    prcc_path = args.global_sensitivity_dir / "global_sensitivity_98x9_prcc.csv"
    if not prcc_path.exists():
        prcc_path = args.global_sensitivity_dir / "global_sensitivity_prcc_auc.csv"
    prcc = pd.read_csv(prcc_path)
    uncertainty = pd.read_csv(args.uncertainty_dir / "trajectory_quantiles.csv")

    rows, link_rows, window_rows = [], [], []
    for parameter, cls in TARGET_CLASSES.items():
        base = combined[combined["parameter"].eq(parameter)]
        local_rank = base["local_sensitivity_rank"].iloc[0] if "local_sensitivity_rank" in base and not base.empty else pd.NA
        local_score = base["local_sensitivity_max_abs_auc"].iloc[0] if "local_sensitivity_max_abs_auc" in base and not base.empty else pd.NA
        top = prcc[prcc["parameter"].eq(parameter)].sort_values("abs_prcc", ascending=False).head(5)
        biomarkers = top["biomarker"].tolist()
        rows.append(
            {
                "parameter": parameter,
                "profile_likelihood_class": cls,
                "local_sensitivity_rank": local_rank,
                "local_sensitivity_max_abs_auc": local_score,
                "top_global_sensitivity_biomarkers": "|".join(biomarkers),
                "reason_selected": "fixed representative 3x3 profile-likelihood example",
            }
        )
        for rank, (_, item) in enumerate(top.iterrows(), start=1):
            link_rows.append({"parameter": parameter, "rank": rank, **item.to_dict()})
            windows = uncertainty[uncertainty["biomarker"].eq(item["biomarker"])].sort_values("interval_width", ascending=False).head(3)
            for wrank, (_, window) in enumerate(windows.iterrows(), start=1):
                window_rows.append(
                    {
                        "parameter": parameter,
                        "biomarker": item["biomarker"],
                        "biomarker_rank": rank,
                        "window_rank": wrank,
                        "day": int(window["day"]),
                        "interval_width": float(window["interval_width"]),
                    }
                )
    targets = pd.DataFrame(rows)
    links = pd.DataFrame(link_rows)
    windows = pd.DataFrame(window_rows)
    final_tables = Path("results_final/tables")
    final_tables.mkdir(parents=True, exist_ok=True)
    targets.to_csv(args.output_dir / "bayesian_target_parameters.csv", index=False)
    targets.to_csv(final_tables / "bayesian_bed_representative_targets.csv", index=False)
    links.to_csv(final_tables / "bed_targeted_parameter_biomarker_time_links.csv", index=False)
    windows.to_csv(final_tables / "bed_targeted_gsa_uncertainty_observation_scenarios.csv", index=False)
    print("9 representative profile-likelihood parameters used")
    print(f"top biomarkers per representative parameter: {targets[['parameter','top_global_sensitivity_biomarkers']].to_dict('records')}")
    print("uncertainty-informed time windows selected")


if __name__ == "__main__":
    main()
