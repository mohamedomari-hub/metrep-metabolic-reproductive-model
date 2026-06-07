# Bayesian Inference And BED Integration

This folder contains the final selected-parameter Bayesian inference layer for
the MetRep portfolio workflow. It uses the fixed representative 3x3
profile-likelihood parameter set:

- practically identifiable: `insulin_glucose_threshold`,
  `inhibin_clearance`, `hp_p4_follicle_scale`
- boundary-limited: `blood_to_liver_glucose_threshold`, `gnrh_clearance`,
  `hp_iof_threshold`
- weak/flat: `insulin_igf_threshold`, `feed_direct_blood_fraction`,
  `lh_basal_release`

The workflow order is:

```text
local sensitivity
-> biological admissibility filtering
-> global sensitivity on admissible ensemble
-> SVD identifiability
-> profile likelihood
-> uncertainty propagation
-> Bayesian inference and BED
```

## Scientific Inputs

The broad `+/-5%` ODE bank is the prior archive for posterior distributions.
The SMC+ML enriched bank is used only after ODE confirmation, and only for
stable BED ranking, smaller MI error bars, and broader admissible-ensemble
coverage. Surrogate-predicted candidates are never treated as scientific truth.

Observable biomarkers used downstream:

```text
FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon
```

Glucagon was excluded from the historical biological admissibility filter but
is retained as an observable biomarker for downstream uncertainty propagation,
global sensitivity, and Bayesian experimental design.

## Scripts

`select_bayesian_targets.py`

Builds the fixed 3x3 target table and links each parameter to global
sensitivity biomarkers and uncertainty-informed time windows.

`reduced_mcmc/run_reduced_mcmc.py`

Runs reduced selected-parameter posterior inference using likelihood weights
on the real broad-prior ODE archive. It does not train a surrogate and does
not run new ODE simulations; posterior support is limited to the saved broad
`+/-5%` archive rows.

`abc_smc/run_abc_smc.py`

Runs archive-based sequential ABC filtering over real ODE rows from the
broad-prior bank. It uses the same independent, global cumulative, high/low
day, and parameter-specific GSA+uncertainty+MI guided observation scenarios
where available.

`compare_bayesian_methods.py`

Combines posterior reweighting, reduced ODE-archive posterior, and
archive-based sequential ABC filtering summaries into a single
method-comparison table and figure.

## Commands

Run after global sensitivity, uncertainty propagation, and targeted BED outputs
exist:

```bash
python analyses/bayesian_inference/select_bayesian_targets.py \
  --profile-dir results_final/tables \
  --global-sensitivity-dir analyses/global_sensitivity/outputs_enriched \
  --uncertainty-dir analyses/uncertainty/outputs_smc_enriched_5pct_glucagon \
  --output-dir analyses/bayesian_inference/outputs
```

```bash
python analyses/bayesian_inference/reduced_mcmc/run_reduced_mcmc.py \
  --target-parameters analyses/bayesian_inference/outputs/bayesian_target_parameters.csv \
  --output-dir analyses/bayesian_inference/outputs/reduced_mcmc \
  --figure-dir analyses/bayesian_inference/figures/reduced_mcmc \
  --seed 42
```

```bash
python analyses/bayesian_inference/abc_smc/run_abc_smc.py \
  --target-parameters analyses/bayesian_inference/outputs/bayesian_target_parameters.csv \
  --output-dir analyses/bayesian_inference/outputs/abc_smc \
  --figure-dir analyses/bayesian_inference/figures/abc_smc \
  --seed 42
```

```bash
python analyses/bayesian_inference/compare_bayesian_methods.py
```

## Caveat

Different perturbation scales are used because each analysis answers a
different question: local sensitivity uses `+1%` one-at-a-time perturbations,
SVD identifiability uses small finite differences, profile likelihood explores
a wider parameter range, global sensitivity uses simulation-bank associations,
uncertainty propagation uses biologically admissible ensembles, and BED uses
prior-based information calculations.
