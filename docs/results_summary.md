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

The central conclusion is that Bayesian updating is broadly consistent with
the identifiability diagnosis: identifiable parameters tend to narrow more
strongly, boundary-limited parameters update partially, and weak/flat
parameters remain comparatively weakly informed.

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

![Biologically admissible ensemble trajectories](../results_final/figures/uncertainty_readme_summary.png)

Biological admissibility filtering restricts ensemble analyses to plausible
ODE trajectories. The filter removes simulations with implausible
endocrine-metabolic behavior before ensemble-level analyses are summarized.

Only ODE-confirmed admissible rows are treated as scientific truth. The
enriched admissible ensemble contains 12,721 simulations. The historical
admissibility filter used FSH, PGF, P4, E2, INH, IGF1, insulin, and glucose.
Glucagon was excluded from that filter but retained downstream as an
observable biomarker for uncertainty propagation, global sensitivity, and BED.

## 3. Global Sensitivity Analysis

PRCC and Spearman summarize monotonic parameter-biomarker AUC associations
across the 12,721-row ODE-confirmed admissible ensemble. They identify
association structures within the admissible ensemble; they are not causal
effects and not Sobol variance decompositions. These links provide biomarker
candidates for BED.

![Strongest PRCC links by biomarker](../results_final/figures/global_sensitivity_heatmap_top_parameters.png)

The compact heatmap highlights the strongest PRCC links by biomarker. Metabolic
threshold and transport parameters preferentially associate with glucose,
insulin, IGF1, and glucagon AUCs within the admissible ensemble. Reproductive
mechanisms preferentially associate with FSH, PGF, P4, E2, and INH AUCs.
Coupled signals can also appear because the model links metabolic and
reproductive endocrine dynamics.

![Representative parameter-biomarker PRCC heatmap](../results_final/figures/global_sensitivity_representative_identifiability_parameters.png)

The representative-parameter heatmap connects the 3 x 3 profile-likelihood
parameter set to observable biomarkers. It supports the BED step by indicating
which biomarkers are plausible candidates for parameter-specific observation
scenarios. Full 98 x 9 PRCC and Spearman maps are retained as supplementary
diagnostic outputs in the Appendix.

## 4. SVD Identifiability

The SVD analysis separates informed from weak local parameter directions. The
singular-value spectrum shows how much information is available in different
parameter combinations, while the compensation and regional maps translate
that structure into practical calibration guidance.

![Singular values](../results_final/figures/identifiability_singular_values.png)

![Dominant compensatory parameter relationships](../results_final/figures/identifiability_compensation_network_core.png)

![Sensitivity-compensation regional map](../results_final/figures/identifiability_sensitivity_vs_nullspace_all_parameters.png)

Influential parameters are not automatically identifiable. A parameter may
move model outputs locally but still be difficult to estimate if other
parameters can compensate for its effect.

The regional map combines output sensitivity with nullspace participation.
Sensitive and weakly compensated parameters are estimation candidates.
Sensitive but highly compensated parameters are influential but risky to
estimate jointly. Weakly sensitive parameters contain limited information
under the current output panel, and strongly compensated parameters indicate
trade-offs among model mechanisms. Representative estimate candidates include
`insulin_glucose_threshold`, `inhibin_clearance`,
`blood_to_liver_glucose_threshold`, and `hp_p4_follicle_scale`.

## 5. Profile Likelihood

![Representative profile likelihood classes](../results_final/figures/profile_likelihood_representative_3x3.png)

Profile likelihood is the nonlinear confirmation step: it tests whether the
SVD diagnosis persists when nuisance parameters can re-adjust. This distinction
is important because local identifiability structure can change when other
parameters compensate over a wider parameter range.

The representative 3 x 3 panel contains:

- Practically identifiable: `insulin_glucose_threshold`,
  `inhibin_clearance`, `hp_p4_follicle_scale`.
- Boundary-limited: `blood_to_liver_glucose_threshold`, `gnrh_clearance`,
  `hp_iof_threshold`.
- Weak/flat: `insulin_igf_threshold`, `feed_direct_blood_fraction`,
  `lh_basal_release`.

The top row shows bounded minima on both sides of the optimum, consistent with
practical identifiability under the selected observables. The middle row shows
local information but incomplete constraint over the explored range, indicating
boundary-limited behavior. The bottom row remains broad, shallow, or flat,
indicating limited learnability under the current observable panel.

## 6. Uncertainty Propagation

![Reproductive uncertainty propagation](../results_final/figures/uncertainty_reproductive.png)

![Metabolic uncertainty propagation](../results_final/figures/uncertainty_metabolic.png)

Uncertainty propagation summarizes admissible trajectory variability using the
ensemble median, nominal trajectory, and 5th-95th percentile range across all
observable biomarkers. Some biomarkers exhibit narrow uncertainty regions,
whereas others show broader admissible variability under the retained
parameterizations.

Informative measurement windows correspond to periods where trajectories
separate sufficiently across plausible parameterizations. BED uses these
uncertainty-rich windows to restrict candidate sampling days, linking
trajectory variability to experimental timing without overstating precision.

## 7. Bayesian Experimental Design

BED links the diagnostic chain:

```text
profile likelihood class
-> GSA-linked biomarkers
-> uncertainty-selected time windows
-> MI-ranked biomarker/day candidates
-> guided posterior update
```

Bayesian experimental design improves parameter learning when informative
model-output relationships exist, but limited posterior contraction may
persist for weakly identifiable or highly compensatory parameters.

![Guided biomarker MI bars](../results_final/figures/bed_guided_parameter_biomarker_mi_bars.png)

The MI biomarker ranking identifies which observable biomarkers are most
informative for each representative parameter. Different parameters favor
different biomarker/day combinations because information content varies across
biomarkers and sampling times.

![Guided posterior updates](../results_final/figures/bed_targeted_gsa_uncertainty_guided_posteriors_3x3.png)

The guided posterior update applies the full chain parameter by parameter.
Practically identifiable parameters show stronger posterior contraction;
boundary-limited parameters show partial or one-sided updates; weak/flat
parameters remain broad when selected observations do not resolve compensating
directions.

![High versus low information day](../results_final/figures/bed_targeted_high_vs_low_information_day_3x3.png)

The high-versus-low information day comparison checks whether MI-ranked days
produce greater or more interpretable posterior movement than less informative
periods. Additional traceability analyses, including independent observation
scenarios, cumulative biomarker acquisition, and candidate-day information
curves, are provided in the Appendix.

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
approximations. They should be interpreted as consistency checks rather than
as a competition for a single best method.

## 9. Overall Scientific Interpretation

Local sensitivity identifies mechanisms that affect outputs. Biological
admissibility filtering restricts ensemble analyses to plausible ODE behavior.
Global sensitivity maps parameters to biomarkers. SVD separates informed from
weak parameter directions. Profile likelihood confirms identifiable,
boundary-limited, and weak/flat cases. Uncertainty propagation identifies
informative time windows. BED selects biomarker-day observations. Bayesian
posterior updating tests whether the selected observations reduce uncertainty.

Overall, Bayesian updating is broadly consistent with the identifiability
diagnosis. Practically identifiable parameters tend to exhibit stronger
posterior contraction, boundary-limited parameters show partial learning, and
weak or compensatory parameters may remain broadly distributed. BED therefore
improves learning when identifiable information exists, but it does not
artificially resolve parameter directions that remain weak under the available
observable panel.

# Appendix — Supporting Diagnostic Results

## Global Sensitivity Supplementary Figures

![Full 98 x 9 PRCC heatmap](../results_final/figures/global_sensitivity_98x9_prcc_heatmap.png)

The full PRCC map preserves the complete parameter-biomarker association
screen across all analyzed parameters and observable biomarkers. It is useful
for auditability and for identifying secondary associations that are not shown
in the compact main-text heatmap.

![Full 98 x 9 Spearman heatmap](../results_final/figures/global_sensitivity_98x9_spearman_heatmap.png)

The Spearman map shows direct rank associations without partialing out the
other parameters. Agreement with PRCC supports robust monotonic relationships;
differences indicate associations that may depend on covariance among sampled
parameters.

## BED Supporting Traceability Figures

![Independent observation posteriors](../results_final/figures/bed_targeted_independent_observation_posteriors_3x3.png)

Independent observation scenarios show how individual selected measurements
reshape representative parameter priors. These panels separate the effect of
single observation scenarios from cumulative or guided acquisition strategies.

![Cumulative biomarker posteriors](../results_final/figures/bed_targeted_cumulative_biomarker_posteriors_3x3.png)

The cumulative acquisition analysis adds biomarkers sequentially according to
their ranking. It provides traceability for how posterior distributions change
as additional measurements are incorporated.

![Guided day MI curves](../results_final/figures/bed_guided_day_mi_curves.png)

The MI day curves show candidate-day information scores within
uncertainty-selected windows. They document why different representative
parameters can receive different biomarker-day observation scenarios.
