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
- `MetRep_Python/model_definition/` contains the translated 22-state non-Dexa
  Python model equations, parameters, scenarios, simulation, analysis, and
  plotting functions.
- `MetRep_Python/model_running/` contains runnable scripts for validation,
  baseline simulation, sensitivity, identifiability, and profile likelihood.
- `Project_Documentation/sensitivity_identifiability_bayesian_design.md`
  summarizes sensitivity, SVD identifiability, profile likelihood, and the
  link to BED.
- `Project_Documentation/Bayesian_Experimental_Design/` documents BED, including a
  GitHub-facing v3 baseline MATLAB port and historical PhD provenance files.
- `Project_Documentation/dexa_extension_plan.md` documents the Dexa extension
  status and plan.
- `results_final/` contains curated tables and figures for public reporting.
- `Project_Documentation/` contains simulation instructions, method notes, BED
  materials, and result interpretation pages.

## Repository Status

The Python implementation currently covers the 22-state metabolic-reproductive
core model. The original MATLAB reference includes the Dexa PK/PD extension
with states 23-25. The Python Dexa extension is planned and documented, but not
yet presented as completed. The Dexa workflow is included in the scientific
story as the validation/extension endpoint.

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

Run the standard scenario simulation plus sensitivity, SVD identifiability, and
uncertainty analyses:

```bash
python MetRep_Python/model_running/09_run_standard_analysis.py
```

Run the Dexa reference simulation:

```text
Open MATLAB from the repository root and run:
BovSys_run_dexa_v3
```

The Python Dexa script currently documents the planned extension and is kept as
a placeholder:

```bash
python MetRep_Python/model_running/07_run_dexa_scenarios.py
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

The original MATLAB model includes a dexamethasone PK/PD extension. In the
project logic, this is the final validation/extension step. The Dexa simulation
result is reported in the Dexa paper; this repository keeps the MATLAB
reference code and documents the planned Python Dexa implementation.

## License And Citation

Add a license after confirming what can be distributed for code, data, and
published-paper material.
