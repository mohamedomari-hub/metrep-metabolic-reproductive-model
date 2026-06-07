# MetRep: Metabolic-Reproductive Mechanistic Model

MetRep is a mechanistic ODE model of coupled metabolic and reproductive
endocrine regulation in cattle. The model links glucose-insulin-IGF metabolism
with reproductive hormones, ovarian dynamics, and cycle-level biomarkers such
as FSH, PGF, P4, E2, INH, IGF1, insulin, glucose, and glucagon.

This repository presents a public, GitHub-facing analysis portfolio around the
model: sensitivity analysis, biological admissibility filtering, global
sensitivity, identifiability, uncertainty propagation, Bayesian inference, and
Bayesian experimental design (BED).

## Why This Project Matters

Mechanistic endocrine-metabolic models are useful because they encode biology,
not only statistical associations. They can represent feedback loops, delayed
responses, compensation between mechanisms, and time-dependent biomarker
dynamics.

The hard part is not only simulating the ODE model. The harder scientific
question is whether the model can be learned from available measurements. Some
parameters strongly affect outputs but can be compensated by others. Some are
practically identifiable, while others are boundary-limited or weakly informed.
That is why this repository connects mechanistic modelling with uncertainty
analysis, identifiability, and Bayesian experimental design.

The final workflow asks which parameters matter, which mechanisms are
learnable, which biomarkers and days are informative, and whether Bayesian
updating actually reduces uncertainty in the directions suggested by
identifiability analysis.

## Scientific Questions

- Which parameters influence metabolic and reproductive biomarker dynamics?
- Which mechanisms are practically identifiable from observable outputs?
- Which parameters are boundary-limited, weakly identifiable, or flat?
- Which biomarkers and sampling days are most informative for learning?
- Can Bayesian updating improve weak parameter directions?
- Does Bayesian design confirm or contradict profile-likelihood diagnosis?

## Final Analysis Workflow

```text
local sensitivity
-> biological admissibility filtering
-> global sensitivity (admissible ensemble)
-> SVD identifiability
-> profile likelihood
-> uncertainty propagation
-> Bayesian inference & BED
```

**Local sensitivity** uses one-at-a-time nominal perturbations to identify
parameters that strongly affect biomarker AUCs near the calibrated model.

**Biological admissibility filtering** removes simulated parameter sets that
produce non-physiological trajectories, so ensemble analyses focus on plausible
ODE behavior.

**Global sensitivity** uses the ODE-confirmed admissible ensemble to estimate
PRCC and Spearman parameter-biomarker AUC associations across all observable
biomarkers.

**SVD identifiability** evaluates local parameter directions, compensation, and
weakly informed combinations in the measurement-output space.

**Profile likelihood** confirms nonlinear practical identifiability for a fixed
3 x 3 representative parameter set.

**Uncertainty propagation** summarizes the median, 5th percentile, 95th
percentile, and nominal trajectory across admissible ODE-confirmed ensembles.

**Bayesian inference & BED** uses broad prior ODE archives and admissible-bank
ranking to test whether selected biomarker-day observations reduce parameter
uncertainty.

## Key Scientific Findings

Many mechanisms are locally influential and several representative parameters
are practically identifiable. Other parameters are boundary-limited or remain
weak/flat because their effects can be compensated or are poorly expressed in
the available observable panel.

Global sensitivity links specific biomarkers to mechanisms. Uncertainty
propagation highlights time windows where plausible trajectories diverge. BED
then proposes informative biomarker-day measurements, including
parameter-specific guided observation scenarios.

The central result is:

```text
Bayesian updating confirmed the identifiability diagnosis rather than
contradicting it.
```

Practically identifiable parameters narrow strongly after informative
observations. Boundary-limited parameters show partial or one-sided learning.
Weak/flat parameters often remain broad even after guided observations. BED
helps prioritize measurements where information exists, but it does not
magically rescue parameters that are practically non-identifiable under the
available output panel.

## Representative Results

### 1. Local Sensitivity And SVD Overview

![Sensitivity top parameters](results_final/figures/sensitivity_top_parameters.png)

Local sensitivity identifies nominally influential mechanisms, especially in
metabolic control and endocrine regulation. This is the first screen: it shows
which parameters can move outputs, but not whether those parameters can be
uniquely estimated.

![SVD singular values](results_final/figures/identifiability_singular_values.png)

The SVD screen shows how much information is available across parameter
directions. Rapidly decaying singular values indicate weak directions and
potential parameter compensation.

### 2. Compensation Network

![Compensation network](results_final/figures/identifiability_compensation_network_core.png)

The compensation network summarizes parameter pairs or groups that can trade
off while preserving similar output behavior. This explains why some sensitive
parameters are still difficult to estimate.

### 3. Profile Likelihood Representative 3 x 3

![Representative profile likelihood classes](results_final/figures/profile_likelihood_representative_3x3.png)

The representative profile likelihood figure separates parameters into
practically identifiable, boundary-limited, and weak/flat classes. This 3 x 3
parameter set is used consistently in the Bayesian/BED portfolio.

### 4. Global Sensitivity Heatmap

![Global sensitivity heatmap](results_final/figures/global_sensitivity_98x9_prcc_heatmap.png)

The PRCC heatmap links each model parameter to observable biomarker AUCs inside
the ODE-confirmed admissible ensemble. These associations are interpretable
screening signals, not formal Sobol decompositions.

### 5. Uncertainty Propagation

![Uncertainty propagation summary](results_final/figures/uncertainty_readme_summary.png)

Uncertainty propagation shows where admissible trajectories diverge over time.
Those windows define biologically meaningful candidate days for experimental
design.

### 6. BED Guided Posterior Updates

![BED guided posterior updates](results_final/figures/bed_targeted_gsa_uncertainty_guided_posteriors_3x3.png)

The guided BED posterior updates connect profile-likelihood class, global
sensitivity biomarkers, uncertainty windows, and MI-based day ranking. Strongly
identifiable parameters narrow most clearly.

### 7. Archive-Based ABC Parameter-Specific Updates

![Archive ABC parameter-specific updates](results_final/figures/abc_smc_parameter_specific_gsa_uncertainty_guided_posteriors_3x3.png)

Archive-based sequential ABC filtering provides a likelihood-free check using
precomputed ODE rows. It supports the same qualitative conclusion: posterior
learning is strongest for parameters already supported by identifiability
diagnostics.

### 8. Bayesian Method Comparison

![Bayesian method comparison](results_final/figures/bayesian_method_comparison_summary.png)

Posterior reweighting, reduced ODE-archive posterior weighting, and
archive-based sequential ABC filtering are compared on the same representative
parameters. Agreement across methods strengthens the final interpretation.

## Repository Structure

```text
analyses/                         Analysis scripts and workflow-specific outputs
docs/                             Methodology, results summary, reproducibility notes
results_final/                    Curated GitHub-facing figures and tables
MetRep_Matlab/                    Original MATLAB reference implementation
MetRep_Python/                    Python model implementation and core scripts
archive_legacy_project_documentation/
                                  Legacy documentation and provenance archive
```

## Reproducibility

The expensive ODE simulation banks and profile-likelihood calculations were
generated before final portfolio curation. Public-facing scripts operate on
trusted ODE archives and curated outputs; they should not require recomputing
large ODE banks or profile likelihoods for routine inspection.

Core scripts are organized by analysis type under `analyses/`. Final curated
figures and tables are copied to `results_final/` so the scientific story can
be inspected without traversing every intermediate output folder.

## Methods Documentation

- [Methodology](docs/methodology.md)
- [Results summary](docs/results_summary.md)
- [Reproducibility notes](docs/reproducibility.md)
- [Model overview](docs/model_overview.md)

## Limitations

- Archive-based sequential ABC filtering is not full adaptive ABC-SMC: it
  filters precomputed ODE rows and does not perturb particles or rerun ODEs.
- The admissible ensemble is restricted to biologically plausible regions and
  should not be interpreted as full population variability.
- Local sensitivity is not global identifiability; an influential parameter can
  still be weakly identifiable if compensation exists.
- PRCC and Spearman summarize monotonic association and do not capture all
  nonlinear or non-monotonic relationships.
- BED identifies informative measurements within the available model outputs;
  it cannot rescue parameters that are structurally weak under the observable
  panel.

## Future Work

- Full adaptive ABC-SMC with new ODE proposals.
- Richer observational panels and experimental schedules.
- Mechanistic-ML acceleration with strict ODE confirmation.
- Prospective experimental prioritization using BED-ranked biomarkers and
  sampling windows.
