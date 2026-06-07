# MetRep Metabolic-Reproductive Model

This repository contains a curated Python implementation of a mechanistic
metabolic-reproductive ODE model originally developed during PhD work, together
with model diagnostics for sensitivity, identifiability, uncertainty
propagation, and Bayesian experimental design. The original MATLAB model is
preserved as a scientific reference.

## Scientific Motivation

MetRep links metabolic regulation with reproductive endocrine dynamics. The
model is used to study which mechanisms control observable biomarkers, which
parameters can be estimated, how robust the calibrated model is, and which
future measurements would be most informative.

## Model And Biomarkers

The mechanistic ODE model couples reproductive hormone regulation, ovarian
dynamics, and glucose-insulin-IGF metabolism. GitHub-facing analyses emphasize
observable biomarkers:

`FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon`.

## Analysis Workflow

```text
Baseline ODE simulation
-> Local sensitivity
-> Biological admissibility filtering
-> Global sensitivity on ODE-confirmed admissible ensemble
-> Identifiability analysis
-> Profile likelihood
-> Uncertainty propagation
-> Bayesian inference and Bayesian experimental design
   -> posterior reweighting, PhD style
   -> reduced ODE-archive posterior inference
   -> archive-based sequential ABC filtering
   -> SMC+ML admissible-bank enrichment for BED stability
```

## Key Results

### Baseline Dynamics

![Baseline selected states](results_final/figures/baseline_selected_states.png)

The baseline simulation reproduces coupled metabolic and reproductive
endocrine dynamics under the non-lactating baseline scenario.

### Local Sensitivity

Local sensitivity uses a `+1%` one-at-a-time parameter perturbation and AUC
endpoints across all 98 parameters. It identifies mechanisms that strongly
affect biomarker exposure near the nominal calibrated model.

### Identifiability

![Representative profile likelihood classes](results_final/figures/profile_likelihood_representative_3x3.png)

SVD screening and profile likelihood separate parameters into practically
identifiable, boundary-limited, and weak/non-identifiable classes.

### Combined Parameter Diagnostics

![Combined parameter diagnostics](results_final/figures/combined_parameter_diagnostics.png)

Representative parameters are compared across local sensitivity, admissible-bank
global association, and identifiability class. Practically identifiable
parameters generally show stronger and more consistent diagnostic signal,
while weak/non-identifiable mechanisms show limited or inconsistent signal.

### Global Sensitivity And Admissible-Bank Association

![Representative global sensitivity](results_final/figures/global_sensitivity_representative_identifiability_parameters.png)

PRCC and Spearman associations summarize the full 98 x 9
parameter-biomarker AUC relationships across the enriched 12,721-row
ODE-confirmed admissible ensemble. They indicate monotonic association, not
strict Sobol variance decomposition. The full-prior bank uses ordinary Monte
Carlo sampling, so its variance-based results are labeled screening rather than
Sobol indices.

### Uncertainty Propagation

![Uncertainty propagation](results_final/figures/uncertainty_readme_summary.png)

Uncertainty propagation uses ODE-confirmed biologically admissible ensembles.
Shading shows the 5th-95th percentile range, the solid blue line shows the
ensemble median, and the dashed black line shows the nominal trajectory. Narrow
admissible banks should be read as local robustness around the calibrated model;
the expanded `+/-5%` enriched bank is used for the final downstream analyses.

### Bayesian Experimental Design

BED scripts and outputs are available under
`analyses/bayesian_experimental_design/`. The final BED layer uses the broad
`+/-5%` prior for posterior plots and the 12,721-row ODE-confirmed enriched
bank for MI stability and ranking robustness. It includes independent
observation scenarios, a global cumulative best-1 through best-9 biomarker
update, high- versus low-information day comparisons, and parameter-specific
GSA + uncertainty + MI guided posterior updates for the fixed 3x3
representative parameter set.

The final Bayesian/BED interpretation is consistent with the profile
likelihood diagnosis: practically identifiable parameters show stronger
posterior narrowing, boundary-limited parameters show partial or one-sided
learning, and weak/flat parameters often remain broad even under guided
observations. BED identifies informative measurements; it does not by itself
rescue structurally weak parameter directions.

### Bayesian Inference

Reduced posterior and archive-based sequential ABC workflows are available under
`analyses/bayesian_inference/`. The reduced posterior script uses likelihood
weights on real broad-prior ODE archive rows and does not use surrogate
predictions. Archive-based sequential ABC filtering uses real ODE archive rows
and does not treat surrogate-predicted candidates as truth.

## Repository Structure

```text
MetRep_Matlab/                  Original MATLAB reference implementation
MetRep_Python/                  Reproducible Python model and core scripts
analyses/                       Sensitivity, identifiability, uncertainty, and BED workflows
results_final/                  Curated publication/GitHub-facing figures and tables
docs/                           Concise model, methodology, results, and reproducibility notes
```

## Reproducibility

Install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Core model and diagnostics:

```bash
python MetRep_Python/scripts/02_run_baseline.py
python MetRep_Python/scripts/04_run_sensitivity.py
python MetRep_Python/scripts/05_run_identifiability.py
python MetRep_Python/scripts/10_run_profile_likelihood.py
python analyses/global_sensitivity/run_global_sensitivity_enriched_98x9.py
python analyses/uncertainty/run_uncertainty_propagation.py
python analyses/model_diagnostics/build_combined_parameter_summary.py
```

Final Bayesian/BED commands are listed in:

- `analyses/bayesian_experimental_design/surrogate_bed/README.md`
- `analyses/bayesian_inference/README.md`

See [docs/reproducibility.md](docs/reproducibility.md) for workflow details.

## Important Methodological Note

Different analyses use different perturbation scales because they answer
different questions:

- Local sensitivity: `+1%` one-at-a-time perturbation near the nominal model.
- SVD identifiability: small numerical derivative step.
- Profile likelihood: wider parameter profiling range.
- Global association and uncertainty: Monte Carlo simulation banks.
- BED: prior-based information calculation.

## Documentation

- [Model overview](docs/model_overview.md)
- [Methodology](docs/methodology.md)
- [Results summary](docs/results_summary.md)
- [Reproducibility](docs/reproducibility.md)
