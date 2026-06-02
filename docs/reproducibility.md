# Reproducibility

## Environment

Install the Python dependencies from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Minimal Reproducibility Path

These commands check that the Python translation is usable:

```bash
python MetRep_Python/scripts/01_validate_against_matlab.py
python MetRep_Python/scripts/02_run_baseline.py
```

`01_validate_against_matlab.py` verifies parameter count, core state order, and
initial conditions against the MATLAB reference metadata.

`02_run_baseline.py` runs a standard 60-day non-Dexa baseline simulation. This
is the first functional smoke test for the translated Python model: if it runs,
the ODE system, parameters, initial conditions, scenario forcing, result export,
and plotting path are all working.

## Analysis Reproduction

Sensitivity:

```bash
python MetRep_Python/scripts/04_run_sensitivity.py
```

SVD identifiability:

```bash
python MetRep_Python/scripts/05_run_identifiability.py \
  --outputs FSH PGF P4 E2 INH IGF1 Insulin Glucose Glucagon \
  --days 50 \
  --dt 2 \
  --prefix structid_50d_measurable
```

Profile likelihood confirmation:

```bash
python MetRep_Python/scripts/10_run_profile_likelihood.py \
  --from-identifiability MetRep_Python/results/tables/structid_50d_measurable_holistic_table.csv \
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

## Results Policy

The full working `MetRep_Python/results/` directory may contain exploratory,
smoke-test, and intermediate outputs. Public releases should use
`results_final/` for curated tables and figures only.
