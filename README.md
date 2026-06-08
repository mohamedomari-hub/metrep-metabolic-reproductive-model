# MetRep: Metabolic-Reproductive Mechanistic Model

MetRep is a mechanistic ordinary differential equation model of coupled
metabolic and reproductive endocrine regulation in cattle. The model connects
glucose-insulin-IGF metabolism with reproductive hormone regulation, ovarian
dynamics, and observable biomarkers including FSH, PGF, P4, E2, INH, IGF1,
insulin, glucose, and glucagon.

The repository presents a public scientific portfolio around the model:
sensitivity analysis, biological admissibility filtering, global sensitivity,
identifiability, uncertainty propagation, Bayesian inference, and Bayesian
experimental design.

## Why This Project Matters

Mechanistic endocrine-metabolic models encode biological feedback rather than
only statistical association. They can represent nonlinear regulation,
compensation between mechanisms, delayed responses, and time-dependent
biomarker dynamics.

The central modelling challenge is parameter uncertainty. Many parameters can
affect trajectories, but not all of them can be learned from the available
observable outputs. Compensation between mechanisms can make influential
parameters practically weak or boundary-limited. This makes identifiability and
experimental design essential parts of the scientific workflow.

The project links mechanistic ODE modelling with Bayesian experimental design
to ask which measurements are informative and whether new observations can
reduce uncertainty in biologically meaningful parameter directions.

## Scientific Questions

- Which mechanisms influence endocrine-metabolic dynamics?
- Which parameters are practically identifiable from observable biomarkers?
- Which parameters are boundary-limited, weakly identifiable, or flat?
- Which biomarkers and sampling windows are informative for parameter
  learning?
- Can Bayesian updating reduce uncertainty in representative parameter
  directions?

## Final Workflow

```text
local sensitivity
-> biological admissibility filtering
-> global sensitivity (admissible ensemble)
-> SVD identifiability
-> profile likelihood
-> uncertainty propagation
-> Bayesian inference & BED
```

Local sensitivity ranks one-at-a-time nominal parameter effects on biomarker
AUC endpoints. It identifies mechanisms that influence outputs near the
calibrated regime.

Biological admissibility filtering removes parameter sets that generate
non-physiological ODE trajectories. This keeps ensemble analyses focused on
plausible endocrine-metabolic behavior.

Global sensitivity uses the ODE-confirmed admissible ensemble to estimate
PRCC and Spearman parameter-biomarker AUC associations across observable
biomarkers.

SVD identifiability analyzes local output-informed parameter directions and
compensatory combinations.

Profile likelihood provides nonlinear practical-identifiability confirmation
for representative parameters.

Uncertainty propagation summarizes admissible trajectory variability and
identifies time windows where plausible trajectories diverge.

Bayesian inference and BED test whether selected biomarker-day observations
reduce parameter uncertainty, using broad prior ODE archives and
admissible-bank ranking for stability.

## Main Scientific Findings

Many parameters are locally influential, and several representative parameters
are practically identifiable. Compensation is also present: some influential
mechanisms remain difficult to estimate because other parameters can reproduce
similar output behavior.

Profile likelihood separates representative parameters into practically
identifiable, boundary-limited, and weak/flat classes. The admissible ensemble
improves biological realism by restricting global sensitivity, uncertainty,
and BED analyses to plausible ODE trajectories.

Global sensitivity links biomarkers to mechanisms. Uncertainty propagation
reveals informative time windows. BED proposes informative biomarker-day
observations and tests posterior narrowing under those measurements.

Bayesian updating confirms the identifiability diagnosis rather than
contradicting it. Practically identifiable parameters narrow strongly,
boundary-limited parameters show partial or one-sided learning, and weak/flat
parameters remain difficult even under guided observations.

## Repository Structure

```text
MetRep_Matlab/  Original MATLAB reference implementation
MetRep_Python/  Python model implementation and analysis scripts
docs/           Public methodology, model overview, and results summary
results_final/  Curated GitHub-facing figures and tables
```

Large generated ODE banks, particle pools, debug folders, and temporary outputs
are not part of the public repository. They should be regenerated or restored
locally when a full workflow rerun is needed.

### Code Organization

The repository keeps Python model code and Python analysis scripts under one
main folder:

- `MetRep_Python/model_definition/` contains reusable model components:
  parameters, ODE equations, simulation helpers, admissibility utilities,
  profile-likelihood utilities, and plotting helpers.
- `MetRep_Python/model_analysis/` contains runnable scripts for model
  validation, baseline simulation, local sensitivity, SVD identifiability,
  profile likelihood, uncertainty propagation, global sensitivity, BED,
  archive posterior analysis, ABC filtering, and method comparison.

This keeps all Python scripts that a reader may run in one clear place while
preserving the distinction between reusable model definitions and analysis
entry points.

## Reproducibility

The repository is organized so the public GitHub version stays lightweight
while preserving the scientific workflow. Large ODE banks, raw particle pools,
debug outputs, and temporary simulation products are intentionally ignored by
Git. Curated figures and tables needed to understand the final results are kept
under `results_final/`.

### Environment Setup

Use a clean Python environment and install the lightweight analysis
dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The MATLAB reference implementation is preserved in `MetRep_Matlab/`. The
Python workflow uses the model and parameter definitions under `MetRep_Python/`
and the final reproducibility scripts under `MetRep_Python/model_analysis/`.

### Repository Logic

The public repository separates three roles:

- `MetRep_Matlab/` and `MetRep_Python/` preserve the model implementation.
- `MetRep_Python/model_analysis/` contains the current reproducibility scripts
  used by the final portfolio.
- `docs/` and `results_final/` contain the public scientific narrative and
  curated outputs.

Local generated material should be written to ignored folders such as
`local_data/` and `local_outputs/`.

### Reproducibility Modes

Two reproducibility modes are supported:

- **Portfolio inspection:** read `docs/`, inspect curated tables and figures in
  `results_final/`, and run lightweight table/plot scripts.
- **Full regeneration:** rebuild ODE banks and rerun expensive ensemble,
  profile-likelihood, uncertainty, BED, and Bayesian workflows locally. This
  requires substantial runtime and is not needed for routine GitHub review.

### Rebuilding Large Simulation Banks

Large banks are not committed. The main broad-prior and enriched-bank paths used
by the reproducibility commands below are:

```text
local_data/phd_bed_bank_5pct_50k_glucagon/
local_data/phd_bed_bank_5pct_50k_glucagon_smc_enriched/
```

The broad `+/-5%` ODE bank is the prior archive for posterior distributions.
The enriched 12,721-row bank is used only after ODE confirmation and supports
global sensitivity, uncertainty propagation, MI/BED ranking stability, and
admissible ensemble coverage. Surrogate-predicted candidates are never treated
as scientific truth.

Example downstream commands, assuming the local ODE-confirmed banks exist:

```bash
python MetRep_Python/model_analysis/global_sensitivity_enriched_98x9.py \
  --enriched-bank-dir local_data/phd_bed_bank_5pct_50k_glucagon_smc_enriched \
  --output-dir local_outputs/global_sensitivity \
  --figure-dir results_final/figures \
  --table-dir results_final/tables
```

```bash
python MetRep_Python/model_analysis/uncertainty_propagation.py \
  --bank-dir local_data/phd_bed_bank_5pct_50k_glucagon_smc_enriched \
  --output-dir local_outputs/uncertainty
```

```bash
python MetRep_Python/model_analysis/targeted_bed_posterior_portfolio.py \
  --broad-prior-dir local_data/phd_bed_bank_5pct_50k_glucagon \
  --enriched-bank-dir local_data/phd_bed_bank_5pct_50k_glucagon_smc_enriched \
  --global-sensitivity-dir local_outputs/global_sensitivity \
  --uncertainty-dir local_outputs/uncertainty \
  --profile-dir results_final/tables \
  --output-dir local_outputs/bed_targeted \
  --figure-dir local_outputs/bed_figures \
  --seed 42
```

```bash
python MetRep_Python/model_analysis/select_bayesian_targets.py \
  --profile-dir results_final/tables \
  --global-sensitivity-dir local_outputs/global_sensitivity \
  --uncertainty-dir local_outputs/uncertainty \
  --output-dir local_outputs/bayesian_targets
```

```bash
python MetRep_Python/model_analysis/reduced_archive_posterior.py \
  --target-parameters local_outputs/bayesian_targets/bayesian_target_parameters.csv \
  --output-dir local_outputs/reduced_archive_posterior \
  --figure-dir local_outputs/reduced_archive_posterior_figures \
  --seed 42
```

```bash
python MetRep_Python/model_analysis/archive_abc_filtering.py \
  --target-parameters local_outputs/bayesian_targets/bayesian_target_parameters.csv \
  --output-dir local_outputs/archive_abc \
  --figure-dir local_outputs/archive_abc_figures \
  --seed 42
```

```bash
python MetRep_Python/model_analysis/compare_bayesian_methods.py
```

### Randomness And Determinism

Scripts that sample, resample, or split data expose a `--seed` argument. The
curated portfolio uses seed `42` unless noted otherwise. ODE-confirmed banks
are deterministic once the saved parameter rows exist, but Monte Carlo
summaries, posterior resampling, and plotting can vary slightly with seed or
with regenerated banks.

### Final Reproducibility Summary

Different perturbation scales are used because each analysis answers a
different question: local sensitivity uses `+1%` one-at-a-time perturbations,
SVD identifiability uses small finite differences, profile likelihood explores
a wider parameter range, global sensitivity uses simulation-bank associations,
uncertainty propagation uses biologically admissible ensembles, and BED uses
prior-based information calculations.

## Documentation

- [Methodology](docs/methodology.md)
- [Results summary](docs/results_summary.md)
- [Model overview](docs/model_overview.md)

## Limitations

- Archive-based sequential ABC filtering is an archive-based approximation,
  not full adaptive ABC-SMC.
- PRCC and Spearman summarize monotonic association and do not capture all
  nonlinear or non-monotonic relationships.
- The admissible ensemble is restricted to biologically plausible trajectories
  and should not be interpreted as full population variability.
- Local sensitivity is not identifiability; influential parameters can remain
  weakly identifiable when compensation exists.

## Future Work

- Full adaptive ABC-SMC with new ODE proposals.
- Richer observational panels and experimental schedules.
- Mechanistic-ML acceleration with strict ODE confirmation.
- Prospective experimental prioritization using BED-ranked biomarkers and
  sampling windows.
