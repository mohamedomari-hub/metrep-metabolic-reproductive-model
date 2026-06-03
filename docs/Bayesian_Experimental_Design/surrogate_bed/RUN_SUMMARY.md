# Surrogate BED Pilot Run

This folder contains a pilot surrogate-assisted BED run generated from the
translated Python ODE model.

## Input Generation

Input tables were generated with:

```bash
python docs/Bayesian_Experimental_Design/surrogate_bed/prepare_surrogate_inputs.py \
  --n-samples 500 \
  --output-dir docs/Bayesian_Experimental_Design/surrogate_bed/input_tables
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
python docs/Bayesian_Experimental_Design/surrogate_bed/surrogate_bed_pipeline.py \
  --parameters-csv docs/Bayesian_Experimental_Design/surrogate_bed/input_tables/prior_parameter_samples.csv \
  --outputs-csv docs/Bayesian_Experimental_Design/surrogate_bed/input_tables/ode_output_features.csv \
  --target-column insulin_glucose_threshold \
  --candidate-map-csv docs/Bayesian_Experimental_Design/surrogate_bed/input_tables/candidate_map.csv \
  --nominal-output-csv docs/Bayesian_Experimental_Design/surrogate_bed/input_tables/nominal_output.csv \
  --admissibility-csv docs/Bayesian_Experimental_Design/surrogate_bed/input_tables/admissibility.csv \
  --output-dir docs/Bayesian_Experimental_Design/surrogate_bed/run_outputs \
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

## Thesis-Style Surrogate Figures

The thesis-style summary figure was generated with:

```bash
python docs/Bayesian_Experimental_Design/surrogate_bed/plot_thesis_style_surrogate.py \
  --benchmark-samples 12 \
  --surrogate-benchmark-repeats 200
```

Generated files:

| File | Content |
|---|---|
| `run_outputs/figures/surrogate_bed_thesis_style_summary.png` | Nominal follicle/P4 trajectory, all-species MI by day, per-species MI at the best day, and posterior curves for selected measurement designs |
| `run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png` | Local ODE-versus-surrogate prediction timing comparison |
| `run_outputs/thesis_style_posterior_curves.csv` | Prior and posterior density curves used in the summary plot |
| `run_outputs/thesis_style_posterior_diagnostics.csv` | Effective sample size for each plotted posterior curve |
| `run_outputs/thesis_style_speed_benchmark.csv` | Timing values for ODE simulation, surrogate prediction, and surrogate training |

The thesis-style plot selected day 88 as the best all-species candidate and
day 56 as the lowest all-species candidate for the target
`insulin_glucose_threshold`.

Local timing benchmark:

| Method | Seconds per sample |
|---|---:|
| ODE | 0.327 |
| Surrogate prediction | 0.00152 |

This corresponds to an observed surrogate prediction speedup of about `215x`
per parameter sample on this machine. Surrogate training took `0.313` seconds
for the 389 admissible pilot samples; training is a one-time cost, while
prediction is the repeated operation used during candidate exploration.

## Interpretation

This run is useful for comparing the Python surrogate workflow against the
thesis-style BED figures. It should still be treated as a pilot result rather
than a final public BED result because the MI estimate continues to increase
with sample size. A stronger reportable run should use more ODE samples and
should show that the top-candidate ranking is stable as the Monte Carlo sample
size increases.
