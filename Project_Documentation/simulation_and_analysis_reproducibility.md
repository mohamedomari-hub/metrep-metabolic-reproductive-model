# Simulation And Analysis Reproducibility

This page lists the main commands needed to reproduce the public Python
simulation workflows and the MATLAB Dexa reference workflow.

## Environment

Run commands from the repository root.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Translation Check

```bash
python MetRep_Python/model_running/01_validate_against_matlab.py
```

This checks the Python parameter count, state order, and initial conditions
against the MATLAB-derived reference metadata.

## Baseline Python Simulation

```bash
python MetRep_Python/model_running/02_run_baseline.py
```

This runs the standard 22-state non-Dexa Python model and writes simulation
outputs and figures under `MetRep_Python/results/`.

## Built-In Scenario Simulations

List available scenarios:

```bash
python MetRep_Python/model_running/08_run_model_scenarios.py --list
```

Run all built-in non-Dexa scenarios:

```bash
python MetRep_Python/model_running/08_run_model_scenarios.py --scenario all
```

Run one named scenario:

```bash
python MetRep_Python/model_running/08_run_model_scenarios.py --scenario baseline_non_lactating
```

Run a custom feeding schedule from CSV:

```bash
python MetRep_Python/model_running/08_run_model_scenarios.py \
  --forcing-csv MetRep_Python/data/my_scenario.csv \
  --mode lactating \
  --name custom_lactating_scenario
```

The custom CSV must contain:

```text
time_days,DMI,Milk
```

## Standard Analysis Workflow

Run scenario simulations plus all-parameter local sensitivity, SVD
identifiability, and uncertainty analysis:

```bash
python MetRep_Python/model_running/09_run_standard_analysis.py
```

## Sensitivity

```bash
python MetRep_Python/model_running/04_run_sensitivity.py
```

## SVD Identifiability

```bash
python MetRep_Python/model_running/05_run_identifiability.py \
  --outputs FSH PGF P4 E2 INH IGF1 Insulin Glucose Glucagon \
  --days 50 \
  --dt 2 \
  --prefix structid_50d_measurable
```

## Profile Likelihood

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

## Dexa Reference Simulation

The Dexa perturbation is currently represented by the original MATLAB
reference implementation. From MATLAB, set the repository root as the working
directory and run:

```matlab
BovSys_run_dexa_v3
```

The Python Dexa workflow is intentionally a placeholder until the 22-state
Python core is fully validated against MATLAB outputs:

```bash
python MetRep_Python/model_running/07_run_dexa_scenarios.py
```

## Results Policy

Curated public outputs are stored in `results_final/`. Development outputs from
rerunning scripts are written under `MetRep_Python/results/` or the relevant
workflow subfolder.
