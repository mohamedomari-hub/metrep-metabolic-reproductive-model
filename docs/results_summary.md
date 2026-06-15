# Results Summary

This document summarizes the curated results of the MetRep endocrine-metabolic
model portfolio. The analysis follows one diagnostic chain:

```text
local sensitivity
-> biological admissibility filtering
-> global sensitivity
-> SVD identifiability
-> profile likelihood
-> uncertainty propagation
-> Bayesian experimental design
-> Bayesian posterior updating
```

The central conclusion is that Bayesian updating confirms the identifiability
diagnosis: identifiable parameters narrow strongly, boundary-limited parameters
update partially, and weak/flat parameters remain weakly informed.

## 1. Local Sensitivity Analysis

![Top local sensitivity parameters](../results_final/figures/sensitivity_top_parameters.png)

AUC-based local sensitivity ranks parameters by cumulative trajectory response
rather than by a single time point. This is important for an oscillatory
endocrine-metabolic model because a single day can fall inside a follicular
phase, luteal phase, transition period, or pulse event.

The strongest local responses involve glucose-insulin regulation, IGF dynamics,
and reproductive endocrine mechanisms. High-ranking examples include
`insulin_glucose_threshold`, `blood_usage_glucose_threshold`,
`blood_to_liver_glucose_threshold`, `igf_clearance`, `cl_to_p4_scale`,
`cl_formation_scale`, `insulin_secretion_max`, and `insulin_clearance`.

This section identifies which mechanisms move the outputs, not which
parameters can be uniquely estimated. Blue bars indicate positive AUC
sensitivity and orange bars indicate negative AUC sensitivity.

## 2. Biological Admissibility Filtering

Biological admissibility filtering restricts ensemble analyses to plausible
ODE trajectories. This prevents global sensitivity, uncertainty propagation,
and BED from being dominated by non-physiological simulations.

Only ODE-confirmed admissible rows are treated as scientific truth. The
enriched admissible ensemble contains 12,721 simulations. The historical
admissibility filter used FSH, PGF, P4, E2, INH, IGF1, insulin, and glucose.
Glucagon was excluded from that filter but retained downstream as an
observable biomarker for uncertainty propagation, global sensitivity, and BED.

## 3. Global Sensitivity Analysis

PRCC and Spearman summarize monotonic parameter-biomarker AUC associations
across the 12,721-row ODE-confirmed admissible ensemble. They identify
parameter-biomarker links; they are not causal effects and not Sobol variance
decompositions. These links provide biomarker candidates for BED.

![Strongest PRCC links by biomarker](../results_final/figures/global_sensitivity_heatmap_top_parameters.png)

The compact heatmap highlights the strongest PRCC links by biomarker. Metabolic
threshold and transport parameters tend to associate with glucose, insulin,
IGF1, and glucagon AUCs. Reproductive endocrine parameters tend to associate
with FSH, PGF, P4, E2, and INH AUCs.

![Representative parameter-biomarker PRCC heatmap](../results_final/figures/global_sensitivity_representative_identifiability_parameters.png)

The representative-parameter heatmap connects the 3 x 3 profile-likelihood
parameter set to observable biomarkers. Full 98 x 9 PRCC and Spearman maps are
retained as supplementary diagnostic outputs.

## 4. SVD Identifiability

The SVD analysis separates informed from weak local parameter directions. The
singular-value spectrum shows how much information is available in different
parameter combinations, while the decision map separates estimate, fix-anchor,
and low-information candidates.

![Singular values](../results_final/figures/identifiability_singular_values.png)

![Identifiability decision map](../results_final/figures/identifiability_decision_map.png)

Influential parameters are not automatically identifiable. A parameter may
move model outputs locally but still be difficult to estimate if other
parameters can compensate for its effect.

The decision map combines output sensitivity with nullspace participation:
estimate candidates are sensitive and weakly compensated; fixed-anchor
candidates remain influential but strongly compensatory; fixed low-information
candidates have low output information in this setting. Representative
estimate candidates include `insulin_glucose_threshold`,
`blood_usage_glucose_threshold`, `inhibin_clearance`,
`blood_to_liver_glucose_threshold`, `hp_p4_follicle_scale`, and
`hm_ih_2_threshold`.

## 5. Profile Likelihood

![Representative profile likelihood classes](../results_final/figures/profile_likelihood_representative_3x3.png)

Profile likelihood is the nonlinear confirmation step: it tests whether the
SVD diagnosis persists when nuisance parameters can re-adjust.

The representative 3 x 3 panel contains:

- Practically identifiable: `insulin_glucose_threshold`,
  `inhibin_clearance`, `hp_p4_follicle_scale`.
- Boundary-limited: `blood_to_liver_glucose_threshold`, `gnrh_clearance`,
  `hp_iof_threshold`.
- Weak/flat: `insulin_igf_threshold`, `feed_direct_blood_fraction`,
  `lh_basal_release`.

The top row shows clear profile curvature and bounded minima. The middle row
shows partial or one-sided support. The bottom row remains broad, shallow, or
flat, indicating limited learnability under the current observable panel.

## 6. Uncertainty Propagation

![Uncertainty propagation summary](../results_final/figures/uncertainty_readme_summary.png)

Uncertainty propagation summarizes admissible trajectory variability using the
ensemble median, nominal trajectory, and 5th-95th percentile range. Informative
windows are times where plausible trajectories separate meaningfully.

BED uses these uncertainty-rich windows to restrict candidate sampling days,
linking plausible trajectory variability to experimental timing.

## 7. Bayesian Experimental Design

BED links the diagnostic chain:

```text
profile likelihood class
-> GSA-linked biomarkers
-> uncertainty-selected time windows
-> MI-ranked biomarker/day candidates
-> guided posterior update
```

BED is not a magic rescue for non-identifiable parameters. It improves
learning where the model-output system contains information. Weak/flat
parameters may remain broad even under guided observations.

![Guided biomarker MI bars](../results_final/figures/bed_guided_parameter_biomarker_mi_bars.png)

The MI biomarker ranking identifies which observable biomarkers are most
informative for each representative parameter. This connects global
sensitivity links to candidate measurements.

![Guided posterior updates](../results_final/figures/bed_targeted_gsa_uncertainty_guided_posteriors_3x3.png)

The guided posterior update applies the full chain parameter by parameter.
Practically identifiable parameters show stronger posterior contraction;
boundary-limited parameters show partial or one-sided updates; weak/flat
parameters remain broad when selected observations do not resolve compensating
directions.

![High versus low information day](../results_final/figures/bed_targeted_high_vs_low_information_day_3x3.png)

The high-versus-low information day comparison checks whether MI-ranked days
produce more interpretable posterior movement than less informative days. The
independent observation, cumulative biomarker acquisition, and MI day-curve
plots are retained as supplementary traceability outputs.

## 8. Bayesian Posterior Updating

![Reduced archive posterior](../results_final/figures/mcmc_prior_vs_posterior_3x3.png)

The reduced ODE-archive posterior uses likelihood weights over precomputed
broad-prior ODE rows. No surrogate-generated trajectories are used as posterior
truth, and the posterior support is limited to the saved archive.

![Archive ABC posteriors](../results_final/figures/abc_smc_parameter_specific_gsa_uncertainty_guided_posteriors_3x3.png)

Archive-based sequential ABC filtering reduces a distance tolerance over
precomputed ODE rows. It is not full adaptive ABC-SMC because particles are
not perturbed and ODEs are not rerun. It is a likelihood-free archive check
using only ODE-confirmed simulations.

| Parameter class | Expected posterior behavior | Interpretation |
|---|---|---|
| Practically identifiable | Strong narrowing | Existing outputs contain constraining information |
| Boundary-limited | Partial or one-sided update | Observations constrain only part of the explored range |
| Weak/flat | Broad posterior remains | Measurements do not resolve compensating directions |

The posterior figures compare contraction behavior across archive-based
approximations. They should be interpreted as consistency checks, not as a
competition for a single best method.

## 9. Integrated Scientific Story

Local sensitivity identifies mechanisms that affect outputs. Biological
admissibility filtering restricts ensemble analyses to plausible ODE behavior.
Global sensitivity maps parameters to biomarkers. SVD separates informed from
weak parameter directions. Profile likelihood confirms identifiable,
boundary-limited, and weak/flat cases. Uncertainty propagation identifies
informative time windows. BED selects biomarker-day observations. Bayesian
posterior updating tests whether the selected observations reduce uncertainty.

The main conclusion is constructive: BED improves learning when identifiable
information exists, but it does not artificially rescue weak or compensatory
parameters. This is exactly why identifiability analysis is needed before
Bayesian experimental design.
