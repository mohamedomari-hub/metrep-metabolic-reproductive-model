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
-> Identifiability analysis
-> Global sensitivity / admissible-bank association
-> Uncertainty propagation
-> Bayesian experimental design
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

PRCC and Spearman associations summarize parameter-biomarker AUC relationships
across biologically admissible simulations. They indicate monotonic association,
not strict Sobol variance decomposition. The full-prior bank uses ordinary
Monte Carlo sampling, so its variance-based results are labeled screening
rather than Sobol indices.

### Uncertainty Propagation

![Uncertainty propagation](results_final/figures/uncertainty_readme_summary.png)

Uncertainty propagation uses the biologically admissible `+/-0.5%` simulation
bank. Shading shows the 5th-95th percentile range, the solid blue line shows the
ensemble median, and the dashed black line shows the nominal trajectory.
Because the ensemble is narrow and filtered, these bands represent local
robustness around the calibrated model rather than full population variability.

### Bayesian Experimental Design

BED scripts and outputs are available under
`analyses/bayesian_experimental_design/`. Final curated BED figures will be
added after the expanded prior-based analysis is regenerated.

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
python analyses/global_sensitivity/run_global_sensitivity.py
python analyses/uncertainty/run_uncertainty_propagation.py
python analyses/model_diagnostics/build_combined_parameter_summary.py
```

See [docs/reproducibility.md](docs/reproducibility.md) for workflow details.

## Important Methodological Note

Different analyses use different perturbation scales because they answer
different questions:

- Local sensitivity: `+1%` one-at-a-time perturbation near the nominal model.
- SVD identifiability: small numerical derivative step.
- Profile likelihood: wider parameter profiling range.
- Global association and uncertainty: Monte Carlo simulation banks.
- BED: prior-based information calculation.

Different perturbation scales are used because each analysis answers a
different question: local sensitivity uses `+1%` one-at-a-time perturbations,
SVD identifiability uses small finite differences, profile likelihood explores
a wider parameter range, global sensitivity uses simulation-bank associations,
uncertainty propagation uses biologically admissible ensembles, and BED uses
prior-based information calculations.

## Documentation

- [Model overview](docs/model_overview.md)
- [Methodology](docs/methodology.md)
- [Results summary](docs/results_summary.md)
- [Reproducibility](docs/reproducibility.md)
