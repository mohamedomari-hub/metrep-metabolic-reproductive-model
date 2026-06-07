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
analyses/       Analysis workflows for sensitivity, identifiability,
                uncertainty, Bayesian inference, and BED
docs/           Public methodology, results summary, and reproducibility notes
results_final/  Curated GitHub-facing figures and tables
MetRep_Matlab/  Original MATLAB reference implementation
MetRep_Python/  Python model implementation and core scripts
```

Legacy project notes and historical generated materials are preserved in
`archive_legacy_project_documentation/` for provenance.

## Reproducibility

Expensive ODE banks and profile-likelihood results were generated before final
portfolio curation. The public workflows operate on trusted ODE archives and
curated result tables, so routine inspection does not require recomputing large
ODE ensembles or profile likelihoods.

The final Bayesian/BED workflows use broad `+/-5%` prior ODE archives for
posterior distributions and the 12,721-row ODE-confirmed admissible bank for
MI/BED ranking stability and admissible ensemble coverage.

## Documentation

- [Methodology](docs/methodology.md)
- [Results summary](docs/results_summary.md)
- [Reproducibility notes](docs/reproducibility.md)
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
