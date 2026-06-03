# Surrogate BED Pilot Run

This folder contains a pilot surrogate-assisted BED run generated from the
translated Python ODE model.

## Input Generation

Input tables were generated with:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/prepare_surrogate_inputs.py \
  --n-samples 500 \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/input_tables
```

The local prior sampled all 98 analyzed parameters uniformly between 0.995 and
1.005 times their nominal value.

Generated input table sizes:

| Table | Size |
|---|---:|
| `prior_parameter_samples.csv` | 500 x 98 |
| `ode_output_features.csv` | 500 x 324 |
| `candidate_map.csv` | 360 x 2 |
| `nominal_output.csv` | 1 x 324 |
| `admissibility.csv` | 500 x 6 |

Biological admissibility retained 389 of 500 samples.

## Surrogate Run

The surrogate pipeline was run with:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/surrogate_bed_pipeline.py \
  --parameters-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables/prior_parameter_samples.csv \
  --outputs-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables/ode_output_features.csv \
  --target-column insulin_glucose_threshold \
  --candidate-map-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables/candidate_map.csv \
  --nominal-output-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables/nominal_output.csv \
  --admissibility-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables/admissibility.csv \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/run_outputs \
  --n-estimators 600 \
  --convergence-sizes 100 200 300 389 \
  --mi-repeats 5 \
  --make-plots
```

## Diagnostic Summary

The run used 389 admissible samples, 98 parameters, 324 output features, and
360 candidate measurement designs.

Median held-out surrogate accuracy:

| Split | Median R2 | Median normalized RMSE |
|---|---:|---:|
| validation | 0.631 | 0.149 |
| test | 0.572 | 0.168 |

Top candidate for the target `insulin_glucose_threshold`:

```text
all_species_day_88
```

Estimated mutual information for that candidate:

| Samples | Mean MI |
|---:|---:|
| 100 | 0.514 |
| 200 | 0.688 |
| 300 | 0.795 |
| 389 | 0.860 |

Posterior effective sample size for the top candidate:

```text
37.45 / 389 = 0.096
```

## Interpretation

This run is useful for comparing the Python surrogate workflow against the
thesis-style BED figures. It should still be treated as a pilot result rather
than a final public BED result because the MI estimate continues to increase
with sample size. A stronger reportable run should use more ODE samples and
should show that the top-candidate ranking is stable as the Monte Carlo sample
size increases.
