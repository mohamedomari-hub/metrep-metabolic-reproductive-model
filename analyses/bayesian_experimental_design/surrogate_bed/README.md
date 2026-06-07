# Bayesian Experimental Design

This folder contains the PhD-faithful Bayesian experimental design workflow and
the final targeted posterior-update layer used in the GitHub portfolio.

The scientific hierarchy is:

```text
broad +/-5% prior ODE bank
-> biological admissibility filter
-> SMC+ML enrichment with ODE confirmation
-> enriched admissible bank for stable MI ranking
-> broad-prior posterior reweighting for uncertainty reduction plots
```

Only ODE-confirmed rows are scientific truth. The enriched 12,721-sample bank is
used for MI ranking stability, smaller error bars, and broader admissible
coverage. Posterior density plots use the broad `+/-5%` prior bank so narrowing
is shown relative to the original thesis prior scale.

## Observable Biomarkers

Downstream BED uses all observable biomarkers stored in the bank:

```text
FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon
```

Glucagon was excluded from the biological admissibility filter but retained as
an observable biomarker for downstream uncertainty propagation, global
sensitivity, and Bayesian experimental design.

## Main Scripts

`prepare_phd_bed_bank.py`

Creates the broad `+/-5%` PhD prior simulation bank from ODE outputs. This is
the expensive ODE-generation step and should not be rerun unless a new bank is
needed.

`matlab_faithful_parameter_bed.py`

Runs the thesis-style density-ratio BED estimator using the ODE-confirmed
admissible bank. It estimates full-vector information `I(Theta; Y)` first, then
parameter-specific information for interpretation.

`run_targeted_bed_posterior_portfolio.py`

Builds final 3x3 target-parameter posterior updates from existing outputs:

- independent observation scenarios
- global cumulative biomarker acquisition, best 1 through best 9
- highest- versus lowest-information day comparison
- parameter-specific GSA + uncertainty + MI guided observation scenarios
- broad `+/-5%` prior to posterior reweighting
- MI stability comparison using the 12k enriched admissible bank

`debug_enriched_bed_posterior.py`

Audits posterior normalization, prior source, day consistency, cumulative
biomarker logic, and enriched-bank stability.

## Final Commands

Use the formatted 12,721-row ODE-confirmed bank for stable BED ranking:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/matlab_faithful_parameter_bed.py \
  --input-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon_smc_enriched \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_results_smc_enriched_5pct_glucagon \
  --seed 42
```

Build the targeted final posterior update portfolio:

```bash
python analyses/bayesian_experimental_design/surrogate_bed/run_targeted_bed_posterior_portfolio.py \
  --broad-prior-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --enriched-bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon_smc_enriched \
  --global-sensitivity-dir analyses/global_sensitivity/outputs_enriched \
  --uncertainty-dir analyses/uncertainty/outputs_smc_enriched_5pct_glucagon \
  --profile-dir results_final/tables \
  --output-dir analyses/bayesian_experimental_design/surrogate_bed/run_outputs_targeted \
  --figure-dir analyses/bayesian_experimental_design/surrogate_bed/run_figures_targeted \
  --seed 42
```

## Outputs

Targeted BED writes analysis-level outputs under:

```text
analyses/bayesian_experimental_design/surrogate_bed/run_outputs_targeted/
analyses/bayesian_experimental_design/surrogate_bed/run_figures_targeted/
```

and curated portfolio copies under:

```text
results_final/tables/
results_final/figures/
```

Key final tables include:

- `bed_guided_update_traceability.csv`
- `bed_guided_parameter_biomarker_mi.csv`
- `bed_guided_parameter_biomarker_day_mi.csv`
- `bed_guided_selected_observation_scenarios.csv`
- `bed_targeted_parameter_biomarker_time_links.csv`
- `bed_targeted_gsa_uncertainty_observation_scenarios.csv`
- `bed_targeted_posterior_narrowing_independent.csv`
- `bed_targeted_posterior_narrowing_cumulative.csv`
- `bed_targeted_high_vs_low_day_narrowing.csv`
- `bed_mi_stability_original_vs_12k.csv`

Key final figures include:

- `bed_mi_ranking_stability_12k.png`
- `bed_guided_parameter_biomarker_mi_bars.png`
- `bed_guided_day_mi_curves.png`
- `bed_targeted_independent_observation_posteriors_3x3.png`
- `bed_targeted_cumulative_biomarker_posteriors_3x3.png`
- `bed_targeted_high_vs_low_information_day_3x3.png`
- `bed_targeted_gsa_uncertainty_guided_posteriors_3x3.png`

The cumulative figure and the guided figure answer different questions. The
cumulative figure applies one global biomarker order to show best-1 through
best-9 acquisition. The guided figure is parameter-specific: each representative
parameter receives its own top 2-3 biomarkers and days selected from GSA links,
uncertainty windows, and MI/BED ranking.

## Caveat

The Python targeted posterior layer is a reproducible portfolio workflow over
saved ODE banks. It does not replace the original MATLAB thesis script; it
curates the same scientific logic into reusable tables and figures while
preserving the broad-prior posterior interpretation.
