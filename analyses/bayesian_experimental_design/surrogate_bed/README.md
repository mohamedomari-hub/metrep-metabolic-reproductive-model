# Surrogate BED Plan

This folder is a future-work scaffold for a reproducible Python surrogate BED
workflow. It is not part of the current curated result set.

## Recommended Strategy

Use the ODE model as the reference model and the surrogate only as an
accelerator. A surrogate result should be reported only after three checks:

1. The surrogate predicts held-out ODE simulations accurately for the measured
   species and sampling days used in BED.
2. The candidate mutual-information ranking is stable as the Monte Carlo sample
   size increases.
3. Simulated trajectories used by the surrogate workflow pass biological
   admissibility checks.

## Surrogate Choice

The recommended first implementation is a PCA-compressed multi-output emulator
with a tree ensemble regressor:

- PCA compresses correlated multi-species/time outputs.
- ExtraTrees or Random Forest provides a robust nonlinear baseline.
- Held-out ODE predictions remain the required validation target.

Gaussian processes can be useful for smaller training sets and uncertainty
quantification, but they scale poorly for large Monte Carlo designs. Neural
surrogates should only be used if enough ODE simulations are available for
training and validation.

## Required Inputs

The scaffold script expects:

- prior parameter samples, one row per parameter set;
- ODE outputs for the same samples, one row per sample and one column per
  output feature;
- optional biological admissibility flags for each sample.

The script does not rerun simulations. Generate these input tables separately
from the ODE model, then use the scaffold to train and validate the surrogate
and compute BED summaries.

## Reporting Rule

Do not include surrogate BED figures in `results_final/` until the validation
metrics, convergence diagnostics, and admissibility filtering are documented.
