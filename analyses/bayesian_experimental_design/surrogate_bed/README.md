# Surrogate-Assisted BED Workflow

This folder contains a Python workflow for comparing surrogate-assisted
Bayesian experimental design (BED) results with the original thesis results.
The script does not run the MetRep ODE model. It uses ODE simulations that were
already generated and stored as tables.

The key rule is:

The ODE model remains the reference. The surrogate is only an accelerator.

## Step 1: Prepare Input Tables

Generate the required CSV files from the translated Python ODE model:

```bash
cd /Users/omari/Documents/Academic/BovSys/MetRep_Model

python analyses/bayesian_experimental_design/surrogate_bed/prepare_surrogate_inputs.py \
  --n-samples 500 \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/input_tables
```

This creates:

```text
analyses/bayesian_experimental_design/surrogate_bed/input_tables/prior_parameter_samples.csv
analyses/bayesian_experimental_design/surrogate_bed/input_tables/ode_output_features.csv
analyses/bayesian_experimental_design/surrogate_bed/input_tables/candidate_map.csv
analyses/bayesian_experimental_design/surrogate_bed/input_tables/nominal_output.csv
analyses/bayesian_experimental_design/surrogate_bed/input_tables/admissibility.csv
```

The default parameter prior is local: each analyzed parameter is sampled
uniformly between 0.995 and 1.005 times its nominal value. This narrow local
prior is intentional because the biological admissibility filter can reject
very wide all-parameter perturbations. For a quick test, use fewer samples such
as `--n-samples 120`. For a more stable BED comparison, increase the sample
count after confirming the workflow runs correctly.

## Step 2: Run The Surrogate BED Pipeline

After the input tables exist, run:

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

Only `--parameters-csv`, `--outputs-csv`, and `--target-column` are strictly
required. The other inputs make the comparison more scientifically defensible.

## Input Tables

### Parameter Table

Rows are Monte Carlo prior samples. Columns are uncertain parameters.

Example:

| insulin_glucose_threshold | insulin_clearance | p4_clearance |
|---:|---:|---:|
| 0.91 | 1.12 | 0.83 |
| 1.04 | 0.96 | 1.18 |

The target parameter named in `--target-column` must be one of these columns.

### Output Table

Rows must match the parameter table. Columns are simulated ODE outputs used as
candidate measurements. A useful naming convention is species plus sampling
day.

Example:

| FSH_day_65 | P4_day_65 | E2_day_65 | Glucose_day_65 |
|---:|---:|---:|---:|
| 1.45 | 0.32 | 0.12 | 0.47 |
| 1.10 | 0.51 | 0.08 | 0.49 |

### Candidate Map

The candidate map defines what each possible design measures. It must contain
two columns:

| candidate | columns |
|---|---|
| FSH_day_65 | FSH_day_65 |
| endocrine_panel_day_65 | FSH_day_65\|P4_day_65\|E2_day_65\|INH_day_65 |

If this file is omitted, every output column is treated as one scalar
candidate design.

### Nominal Output Table

This optional one-row table contains the nominal/reference prediction. It is
used to generate a fixed synthetic observation for posterior plots.

### Admissibility Table

This optional table contains one Boolean column:

| admissible |
|---|
| true |
| false |

Rows marked `false` are removed before training, validation, posterior
calculation, and mutual-information calculation.

## Method

### Surrogate Model

Let the uncertain parameter vector be:

$$
\theta_i \in \mathbb{R}^{p}
$$

and let the ODE output feature vector be:

$$
y_i = y(\theta_i) \in \mathbb{R}^{m}
$$

The surrogate learns an approximation:

$$
\widehat{y}(\theta_i) \approx y(\theta_i)
$$

The output matrix is first standardized and compressed with principal component
analysis (PCA):

$$
\widetilde{Y} = \mathrm{scale}(Y)
$$

$$
\widetilde{Y} \approx A C^{T}
$$

where `A` contains PCA scores and `C` contains PCA loading vectors. The tree
ensemble is trained to predict the PCA scores from the parameter vector:

$$
\widehat{A} = f(\theta)
$$

The predicted output is reconstructed by inverse PCA and inverse scaling:

$$
\widehat{Y} =
\mathrm{scale}^{-1}
\left(
  \widehat{A} C^{T}
\right)
$$

This is more stable than fitting one independent model per output because the
multi-species trajectory is treated as one correlated output object.

### Training, Validation, And Test Split

The samples are split into three groups:

- training samples: fit the surrogate
- validation samples: tune/check the workflow and learning curve
- test samples: final held-out accuracy check

For each output feature, the script reports:

$$
\mathrm{RMSE}_{j}
=
\sqrt{
  \frac{1}{n}
  \sum_{i=1}^{n}
  \left(
    y_{ij} - \widehat{y}_{ij}
  \right)^2
}
$$

and:

$$
\mathrm{NRMSE}_{j}
=
\frac{\mathrm{RMSE}_{j}}
{\max_i(y_{ij}) - \min_i(y_{ij})}
$$

The script also reports `R2`, mean absolute error, and a learning curve. The
learning curve repeats training with increasing fractions of the training data.
If the validation error is still decreasing strongly at the largest training
size, more ODE simulations are needed before reporting surrogate BED results.

### Synthetic Observation And Posterior

For posterior plots, the script uses one fixed synthetic observation vector:

$$
z_{\mathrm{obs}}
=
\mu + \sigma \odot \epsilon
$$

where:

$$
\epsilon \sim \mathcal{N}(0,I)
$$

and:

$$
\sigma_j =
\max
\left(
  r |\mu_j|,
  \sigma_{\min}
\right)
$$

Here, `r` is the relative noise level and `sigma_min` is a small noise floor.
The default relative noise is 0.05.

For each prior sample, the Gaussian log-likelihood is:

$$
\ell_i
=
-
\frac{1}{2}
\sum_j
\left(
  \frac{
    z_{\mathrm{obs},j} - \widehat{y}_{ij}
  }
  {\sigma_j}
\right)^2
$$

The likelihood weights are stabilized and normalized:

$$
\widetilde{w}_i =
\exp
\left(
  \ell_i - \max_k \ell_k
\right)
$$

$$
w_i =
\frac{\widetilde{w}_i}
{\sum_k \widetilde{w}_k}
$$

The posterior density of the target parameter is estimated with a weighted
kernel density estimate:

$$
p(\theta_p \mid z_{\mathrm{obs}})
\approx
\sum_i
w_i
K_h
\left(
  \theta_p - \theta_{i,p}
\right)
$$

The effective sample size is also reported:

$$
\mathrm{ESS}
=
\frac{1}
{\sum_i w_i^2}
$$

A very small ESS means that the posterior is dominated by too few samples. In
that case, increase the Monte Carlo sample size or use a less restrictive
noise model before interpreting the posterior.

### Mutual Information

For each candidate design, the script estimates how informative the simulated
measurement is about the target parameter.

Let:

$$
W_i = \theta_{i,p}
$$

be the target parameter value for sample `i`, and let:

$$
Z_i =
\widehat{y}_{i,D}
$$

be the surrogate-predicted output vector for candidate design `D`.

The quantity of interest is:

$$
I(W;Z)
=
\int
\int
p(w,z)
\log
\left(
  \frac{p(w,z)}
  {p(w)p(z)}
\right)
dw dz
$$

The script estimates this using a k-nearest-neighbor mutual-information
estimator with the Chebyshev distance. Before neighbor counting, `W` and `Z`
are rank-Gaussianized so that species with large numerical units do not
dominate the distance calculation.

For reporting, the important point is practical:

- high `I(W;Z)` means the candidate measurement is expected to reduce
  uncertainty about the target parameter;
- low `I(W;Z)` means the candidate measurement is weakly informative for that
  target;
- stable rankings across increasing sample sizes are more important than a
  single absolute MI number.

### Convergence

The script estimates MI repeatedly at increasing Monte Carlo sample sizes:

```text
250, 500, 1000, 2000, 5000, 10000
```

For each size, it draws repeated subsamples and reports the mean and standard
deviation of the MI estimate. A candidate design should only be interpreted if
its MI curve is reasonably stable and its ranking does not change strongly
with more samples.

## Output Files

The script writes:

| File | Meaning |
|---|---|
| `surrogate_validation_metrics.csv` | Held-out validation/test accuracy for every output feature |
| `surrogate_learning_curve.csv` | Whether more ODE training samples are still needed |
| `surrogate_predicted_outputs.csv` | Surrogate predictions for all admissible samples |
| `mi_candidate_ranking.csv` | Candidate designs ranked by estimated mutual information |
| `mi_convergence.csv` | Repeated MI estimates across Monte Carlo sample sizes |
| `posterior_prior_comparison.csv` | Prior and posterior KDE for the selected candidate |
| `posterior_diagnostics.csv` | Effective sample size and posterior diagnostic values |
| `run_metadata.json` | Settings used for the run |

With `--make-plots`, the script also writes:

| Figure | Meaning |
|---|---|
| `surrogate_validation_error.png` | Worst held-out prediction errors |
| `surrogate_learning_curve.png` | Whether surrogate accuracy improves with more training data |
| `mi_candidate_ranking.png` | Most informative candidate measurements |
| `mi_convergence.png` | Stability of MI estimates with increasing Monte Carlo size |
| `posterior_comparison.png` | Prior versus posterior density for the selected candidate |

## Thesis-Style Summary Plot

After the surrogate BED pipeline has been run, create a thesis-style summary
figure and an ODE-versus-surrogate speed comparison with:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/plot_thesis_style_surrogate.py \
  --benchmark-samples 12 \
  --surrogate-benchmark-repeats 200
```

This writes:

| Figure or table | Meaning |
|---|---|
| `run_outputs/figures/surrogate_bed_thesis_style_summary.png` | Nominal follicle/P4 trajectory, all-species MI by day, per-species MI, and posterior comparison curves |
| `run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png` | Measured ODE solve time versus surrogate prediction time |
| `run_outputs/thesis_style_speed_benchmark.csv` | Timing values behind the speed plot |
| `run_outputs/thesis_style_posterior_curves.csv` | Density curves behind the posterior panel |
| `run_outputs/thesis_style_posterior_diagnostics.csv` | Effective sample size for each posterior curve |

## Reporting Rule

Do not include surrogate BED figures in `results_final/` until these checks are
acceptable:

- held-out surrogate prediction errors are small for the measured species and
  days used in BED;
- learning-curve error has approximately plateaued;
- MI rankings are stable across increasing Monte Carlo sample sizes;
- posterior effective sample size is not too small;
- biological admissibility filtering is applied and documented.

## Faithful PhD Parameter BED

The final scientific workflow uses the original independent uniform
parameter prior of nominal `+/-5%`. The earlier `input_tables_large/` bank
uses `+/-0.5%` and is a pilot only; the faithful analyzer rejects it.

The workflow calculates:

1. global single-day information `I(Theta; Y_day)` using the full parameter
   vector and all eight MATLAB biomarkers;
2. individual biomarker information at the best global day;
3. cumulative biomarker acquisition at the global best day, ordered by
   descending full-vector `I(Theta; Y_i)`;
4. normalized parameter posterior curves for the representative
   identifiability classes, with a highest-versus-lowest-information-day
   comparison;
5. parameter-specific `I(theta_i; Y_day)` as secondary interpretation.

The posterior calculation follows the thesis normalization:

```text
p(theta_i | z*) = p(theta_i, z*) / integral p(theta_i, z*) d theta_i
```

`posterior_diagnostics.csv` reports the numerical density integral, effective
sample size, and posterior-to-prior variance ratio for every cumulative
measurement set.

`posterior_learning_by_identifiability_class.csv` summarizes the final
all-biomarker posterior at the globally highest-information day. It provides
the direct BED-to-identifiability comparison without assuming that every
class must show its expected narrowing pattern.

Generate the faithful bank once. The generator checkpoints every 100
simulations and resumes from the existing checkpoint after interruption:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/prepare_phd_bed_bank.py \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct \
  --n-samples 10000
```

Then run the faithful MATLAB-style density-ratio BED:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/matlab_faithful_parameter_bed.py \
  --parameters-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct/prior_parameter_samples.csv \
  --outputs-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct/ode_output_features.csv \
  --nominal-output-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct/nominal_output.csv \
  --measurement-noise-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct/measurement_noise.csv \
  --admissibility-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct/admissibility.csv \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_results
```

By default, the analyzer uses every biologically admissible sample as both
the evaluation ensemble and Gaussian-mixture reference ensemble, matching
the original MATLAB workflow. Runtime-limiting subsampling options exist only
for diagnostics and must not be used for the final PhD reproduction.

### Earlier 7,874-Sample Bank Comparison

The earlier `+/-0.5%` bank can be analyzed with the same BED outputs and
figures for direct comparison. It remains scientifically separate from the
faithful `+/-5%` PhD result:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/matlab_faithful_parameter_bed.py \
  --parameters-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables_large/prior_parameter_samples.csv \
  --outputs-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables_large/ode_output_features.csv \
  --nominal-output-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables_large/nominal_output.csv \
  --measurement-noise-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct/measurement_noise.csv \
  --admissibility-csv analyses/bayesian_experimental_design/surrogate_bed/input_tables_large/admissibility.csv \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/bed_results_first_bank_05pct \
  --required-prior-half-range 0.005 \
  --prior-range-tolerance 0.001 \
  --analysis-label "Earlier +/-0.5% bank comparison"
```

The figures state the admissible-bank, MI-evaluation, and mixture-reference
sample counts. Their global BED barplots also include a KDE-smoothed visual
guide; the bars remain the density-ratio MI estimates used for ranking.

### Extended 50,000-Simulation PhD Bank With Glucagon

This extended bank preserves the original `+/-5%` prior and original
eight-biomarker biological admissibility filter, while storing Glucagon as a
ninth BED biomarker. It is kept separate from the original eight-biomarker
PhD reproduction:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/prepare_phd_bed_bank.py \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --n-samples 50000 \
  --checkpoint-every 100 \
  --species FSH PGF P4 E2 INH IGF1 Insulin Glucose Glucagon
```

After generation completes, run the extended MATLAB-style BED:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/matlab_faithful_parameter_bed.py \
  --parameters-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon/prior_parameter_samples.csv \
  --outputs-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon/ode_output_features.csv \
  --nominal-output-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon/nominal_output.csv \
  --measurement-noise-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon/measurement_noise.csv \
  --admissibility-csv analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon/admissibility.csv \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_results_5pct_50k_glucagon \
  --species FSH PGF P4 E2 INH IGF1 Insulin Glucose Glucagon \
  --analysis-label "Extended +/-5% PhD bank with Glucagon"
```
