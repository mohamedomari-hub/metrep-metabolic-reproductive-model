# BovSys / MetRep Model

Mechanistic metabolic-reproductive modeling, identifiability analysis,
Bayesian experimental design, and dexamethasone perturbation validation.

This repository is a reproducible companion to PhD work on the BovSys/MetRep
dairy cow model. It preserves the original MATLAB reference implementation,
provides an open Python translation for users without MATLAB, and organizes the
analysis logic from model construction through experimental-design improvement
and Dexa validation.

## Scientific Logic

The project follows this modeling arc:

```text
1. Build the MetRep mechanistic model
   -> couple reproductive endocrine dynamics with glucose-insulin metabolism

2. Translate the model from MATLAB to Python
   -> make the model accessible and reproducible without a MATLAB license

3. Run classical model analysis
   -> local sensitivity analysis
   -> SVD/local identifiability screening
   -> profile-likelihood confirmation

4. Identify information gaps
   -> which parameters are sensitive, estimable, weak, compensatory, or fixed

5. Use Bayesian experimental design
   -> ask which sampling days and measured species would improve information
      for weak or non-identifiable model directions

6. Extend and validate with Dexa perturbation
   -> use the dexamethasone scenario as an external pharmacological challenge
      to test whether the model reproduces expected metabolic responses
```

In this structure, Bayesian experimental design is not a separate add-on. It is
the answer to the identifiability problem: after finding weakly informed
parameters, BED asks how future experiments should be designed to make the model
more informative. The Dexa module then acts as a perturbation-based validation
endpoint.

## What Is Included

- `MetRep_Matlab/` contains the original MATLAB reference implementation.
- `MetRep_Python/model_definition/` contains the translated 22-state Python
  core model equations, parameters, scenarios, simulation, analysis, plotting
  functions, and the optional 25-state Dexa extension.
- `MetRep_Python/model_running/` contains runnable scripts for validation,
  baseline simulation, sensitivity, identifiability, and profile likelihood.
- `Project_Documentation/sensitivity_identifiability_bayesian_design.md`
  summarizes sensitivity, SVD identifiability, profile likelihood, and the
  link to BED.
- `Project_Documentation/Bayesian_Experimental_Design/` documents BED, including a
  GitHub-facing v3 baseline MATLAB port and historical PhD provenance files.
- `Project_Documentation/dexa_python_implementation.md` documents the Python
  Dexa implementation and how it is kept separate from baseline analyses.
- `results_final/` contains curated tables and figures for public reporting.
- `Project_Documentation/` contains simulation instructions, method notes, BED
  materials, and result interpretation pages.

## Repository Status

The Python implementation covers the 22-state metabolic-reproductive core
model and an optional Dexa perturbation workflow with states 23-25. Ordinary
scenario simulations, sensitivity analysis, SVD identifiability, and profile
likelihood use the non-Dexa core parameter set by default. Dexa is activated
only through the dedicated Dexa runner.

## Core Workflow

The recommended Python workflow is:

```text
01_validate_against_matlab.py
  Confirms Python parameters, state order, and initial conditions match MATLAB.

02_run_baseline.py
  Runs the baseline Python MetRep simulation and saves core model outputs.

04_run_sensitivity.py
  Measures local parameter influence on selected outputs.

05_run_identifiability.py
  Uses SVD of the sensitivity matrix to classify Estimate / Fix parameters.

10_run_profile_likelihood.py
  Confirms practical identifiability for selected parameters.
```

The analyses are intended as a staged workflow:

```text
sensitivity analysis
-> local SVD identifiability screen
-> profile likelihood confirmation
```

Sensitivity analysis identifies high-impact parameters. SVD then evaluates
whether high-impact parameters are separable or compensatory. Profile
likelihood provides a nonlinear practical-identifiability confirmation for the
selected parameters.

The mathematical definitions for sensitivity, SVD/nullspace identifiability,
profile likelihood, and Bayesian experimental design are documented in:

- `Project_Documentation/sensitivity_identifiability_bayesian_design.md`
- `Project_Documentation/results_sensitivity_identifiability_bayesian_design.md`

## Quick Start

Create an environment and install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Validate the Python translation metadata:

```bash
python MetRep_Python/model_running/01_validate_against_matlab.py
```

Run the baseline non-Dexa simulation:

```bash
python MetRep_Python/model_running/02_run_baseline.py
```

List all built-in Python model scenarios:

```bash
python MetRep_Python/model_running/08_run_model_scenarios.py --list
```

Run all built-in non-Dexa Python scenarios:

```bash
python MetRep_Python/model_running/08_run_model_scenarios.py --scenario all
```

Run the optional Python Dexa perturbation:

```bash
python MetRep_Python/model_running/07_run_dexa_scenarios.py \
  --scenario baseline_non_lactating \
  --days 3 \
  --dose-day 0 \
  --figure-dir results_final/figures \
  --table-dir results_final/tables \
  --prefix dexa_non_lactating_standard_3d
```

Run the standard scenario simulation plus sensitivity, SVD identifiability, and
uncertainty analyses:

```bash
python MetRep_Python/model_running/09_run_standard_analysis.py
```

Run the MATLAB Dexa reference simulation:

```text
Open MATLAB from the repository root and run:
BovSys_run_dexa_v3
```

Run the main identifiability screen:

```bash
python MetRep_Python/model_running/05_run_identifiability.py \
  --outputs FSH PGF P4 E2 INH IGF1 Insulin Glucose Glucagon \
  --days 50 \
  --dt 2 \
  --prefix structid_50d_measurable
```

Run the profile likelihood confirmation:

```bash
python MetRep_Python/model_running/10_run_profile_likelihood.py \
  --from-identifiability results_final/tables/structid_50d_measurable_holistic_table.csv \
  --profile-class estimate \
  --max-profile-params 9 \
  --max-nuisance-params 8 \
  --outputs FSH PGF P4 E2 INH IGF1 Insulin Glucose Glucagon \
  --days 50 \
  --grid-low 0.80 \
  --grid-high 1.20 \
  --grid-points 9 \
  --maxiter 20 \
  --admissible-only \
  --rho-min 0.65 \
  --gamma-max 0.45 \
  --kappa-max 0.45 \
  --max-shift-days 7 \
  --profile-scale log1p \
  --trajectory-parameter none \
  --prefix profile_50d_balanced_relaxed
```

## Current Identifiability Result

Using measurable outputs
`FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon`, the current
profile-likelihood confirmation gives:

- 60 profiled parameters
- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable
- 2 flat/non-identifiable

In plain language, 51 parameters had profile curves with clear enough minima
inside the tested range. Six had best fits at the edge of the tested range, one
was only weakly bounded, and two stayed too flat to support reliable estimation
from the current output panel.

See `Project_Documentation/sensitivity_identifiability_bayesian_design.md`,
`Project_Documentation/results_sensitivity_identifiability_bayesian_design.md`,
`Project_Documentation/plot_interpretation_guide.md`, and `results_final/` for
the curated summary.

## Bayesian Experimental Design

The BED implementation is kept as MATLAB reference code in
`Project_Documentation/Bayesian_Experimental_Design/matlab_original/`. For GitHub, the
single recommended entry point is `BED_1M_ALL.m`, which calls
`BovSys_run_v3_baseline.m` and uses the published v3 model equations with Dexa
PK/PD switched off.

The full BED result is reported in the PhD thesis. This repository keeps the
MATLAB BED code and methodology notes, but does not currently include BED/RF
surrogate figures in `results_final` because that workflow still needs a
documented public validation path. A compact Python BED reproduction can be
added later using the translated Python model.

For all simulation and analysis commands, see
`Project_Documentation/simulation_and_analysis_reproducibility.md`.

## Dexa Perturbation Validation

The original MATLAB model includes a dexamethasone PK/PD extension. The Python
translation implements the same optional three-state PK/PD structure and can
compare Dexa trajectories against the matching no-Dexa baseline. The Dexa
simulation result is reported in the Dexa paper; this repository keeps the
MATLAB reference code and a Python runner for reproducible perturbation tests.

Dexa PK/PD constants are not included in the sensitivity, SVD identifiability,
or profile-likelihood parameter list. Those analyses remain focused on the
98-parameter non-Dexa core model.

## License And Citation

Add a license after confirming what can be distributed for code, data, and
published-paper material.
