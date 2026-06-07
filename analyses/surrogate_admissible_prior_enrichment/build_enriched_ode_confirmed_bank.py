"""Build a downstream bank from original plus SMC ODE-confirmed admissible rows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "MetRep_Python"
sys.path.insert(0, str(PYTHON_ROOT))

from model_definition.parameters import PARAMETERS


DEFAULT_ORIGINAL = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon"
DEFAULT_SMC = Path(__file__).resolve().parent / "outputs/smc_ml_enrichment"
DEFAULT_OUTPUT = ROOT / "analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon_smc_enriched"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create enriched ODE-confirmed admissible bank for downstream analyses.")
    parser.add_argument("--original-bank-dir", type=Path, default=DEFAULT_ORIGINAL)
    parser.add_argument("--smc-dir", type=Path, default=DEFAULT_SMC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def observable_columns(outputs: pd.DataFrame) -> list[str]:
    pattern = re.compile(r"^(.+)_day_(\d+)$")
    return [column for column in outputs.columns if pattern.match(column)]


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    parameter_columns = [parameter.name for parameter in PARAMETERS]
    original_parameters = pd.read_csv(args.original_bank_dir / "prior_parameter_samples.csv")
    original_outputs = pd.read_csv(args.original_bank_dir / "ode_output_features.csv")
    original_admissibility = pd.read_csv(args.original_bank_dir / "admissibility.csv")
    original_keep = original_admissibility["admissible"].astype(bool) & original_outputs.notna().all(axis=1)

    confirmed_path = args.smc_dir / "smc_confirmed_admissible_parameters.csv"
    confirmed = pd.read_csv(confirmed_path)
    confirmed = confirmed.loc[confirmed["ode_confirmed_admissible"].astype(bool)].reset_index(drop=True)

    output_columns = observable_columns(original_outputs)
    missing_parameters = sorted(set(parameter_columns).difference(confirmed.columns))
    missing_outputs = sorted(set(output_columns).difference(confirmed.columns))
    if missing_parameters:
        raise ValueError(f"SMC confirmed file is missing parameter columns: {missing_parameters[:5]}")
    if missing_outputs:
        raise ValueError(f"SMC confirmed file is missing observable output columns: {missing_outputs[:5]}")

    enriched_parameters = pd.concat(
        [
            original_parameters.loc[original_keep, parameter_columns],
            confirmed[parameter_columns],
        ],
        ignore_index=True,
    )
    enriched_outputs = pd.concat(
        [
            original_outputs.loc[original_keep, output_columns],
            confirmed[output_columns],
        ],
        ignore_index=True,
    )
    smc_admissibility = pd.DataFrame(
        {
            "admissible": True,
            "penalty": confirmed["ode_penalty"].to_numpy(float),
            "min_correlation": confirmed["ode_min_correlation"].to_numpy(float),
            "max_average_difference": confirmed["ode_max_average_difference"].to_numpy(float),
            "max_norm_difference": confirmed["ode_max_norm_difference"].to_numpy(float),
            "worst_species": confirmed["ode_worst_species"].astype(str).to_numpy(),
        }
    )
    enriched_admissibility = pd.concat(
        [
            original_admissibility.loc[
                original_keep,
                ["admissible", "penalty", "min_correlation", "max_average_difference", "max_norm_difference", "worst_species"],
            ],
            smc_admissibility,
        ],
        ignore_index=True,
    )
    provenance = pd.DataFrame(
        {
            "row_index": np.arange(len(enriched_parameters)),
            "source": ["original_50k_ode_confirmed"] * int(original_keep.sum()) + ["smc_ode_confirmed"] * len(confirmed),
            "source_round": [0] * int(original_keep.sum()) + confirmed["round"].astype(int).tolist(),
        }
    )

    duplicate_mask = enriched_parameters.duplicated(keep="first")
    if duplicate_mask.any():
        keep = ~duplicate_mask
        enriched_parameters = enriched_parameters.loc[keep].reset_index(drop=True)
        enriched_outputs = enriched_outputs.loc[keep].reset_index(drop=True)
        enriched_admissibility = enriched_admissibility.loc[keep].reset_index(drop=True)
        provenance = provenance.loc[keep].reset_index(drop=True)
        provenance["row_index"] = np.arange(len(provenance))

    enriched_parameters.to_csv(args.output_dir / "prior_parameter_samples.csv", index=False)
    enriched_outputs.to_csv(args.output_dir / "ode_output_features.csv", index=False)
    enriched_admissibility.to_csv(args.output_dir / "admissibility.csv", index=False)
    provenance.to_csv(args.output_dir / "sample_provenance.csv", index=False)

    for filename in ["nominal_output.csv", "measurement_noise.csv"]:
        shutil.copy2(args.original_bank_dir / filename, args.output_dir / filename)

    original_metadata_path = args.original_bank_dir / "bank_metadata.json"
    original_metadata = json.loads(original_metadata_path.read_text()) if original_metadata_path.exists() else {}
    metadata = {
        **original_metadata,
        "workflow": "SMC-enriched ODE-confirmed admissible bank for downstream analyses",
        "original_admissible_rows": int(original_keep.sum()),
        "smc_ode_confirmed_rows": int(len(confirmed)),
        "duplicate_parameter_rows_removed": int(duplicate_mask.sum()),
        "admissible_rows": int(len(enriched_parameters)),
        "all_rows_are_ode_confirmed_admissible": True,
        "admissibility_species": original_metadata.get(
            "admissibility_species",
            ["FSH", "PGF", "P4", "E2", "INH", "IGF1", "Insulin", "Glucose"],
        ),
        "observable_species": sorted({column.rsplit("_day_", 1)[0] for column in output_columns}),
        "glucagon_note": (
            "Glucagon was excluded from the biological admissibility filter but retained as an "
            "observable biomarker for downstream uncertainty propagation, global sensitivity, and BED."
        ),
    }
    (args.output_dir / "bank_metadata.json").write_text(json.dumps(metadata, indent=2))
    pd.DataFrame(
        [
            {
                "original_admissible_rows": int(original_keep.sum()),
                "smc_ode_confirmed_rows": int(len(confirmed)),
                "duplicate_parameter_rows_removed": int(duplicate_mask.sum()),
                "enriched_admissible_rows": int(len(enriched_parameters)),
                "observable_biomarkers": "|".join(metadata["observable_species"]),
            }
        ]
    ).to_csv(args.output_dir / "enriched_bank_summary.csv", index=False)

    print(f"Saved enriched ODE-confirmed bank to {args.output_dir}")
    print(f"Original admissible rows: {int(original_keep.sum())}")
    print(f"SMC ODE-confirmed rows: {len(confirmed)}")
    print(f"Duplicate parameter rows removed: {int(duplicate_mask.sum())}")
    print(f"Enriched admissible rows: {len(enriched_parameters)}")


if __name__ == "__main__":
    main()
