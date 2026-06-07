# Results Summary

This document summarizes the curated scientific results for the MetRep
metabolic-reproductive model portfolio. The results are organized as a single
diagnostic chain:

```text
sensitivity
-> compensation
-> identifiability
-> uncertainty
-> BED
-> posterior learning
```

The central scientific conclusion is that Bayesian updating confirms the
identifiability diagnosis rather than contradicting it.

## 1. Local Sensitivity Analysis

![Top local sensitivity parameters](../results_final/figures/sensitivity_top_parameters.png)

Local sensitivity ranks parameters by one-at-a-time AUC response near the
calibrated parameter vector. The strongest effects involve metabolic control,
glucose-insulin regulation, IGF dynamics, and reproductive endocrine
mechanisms. This confirms that the selected observable panel responds to both
metabolic and reproductive subsystems.

Instead of evaluating sensitivity at a single arbitrary time point, the
AUC-based score summarizes the cumulative trajectory response over the
simulation window. This is important for an oscillatory endocrine-metabolic
model because a single time point can be misleading if the system is in a
follicular phase, luteal phase, transition period, or pulse event.

The strongest AUC-level responses are concentrated in glucose-insulin
regulation, IGF dynamics, and reproductive steroid-related parameters. The
highest-ranking parameters include `insulin_glucose_threshold`,
`blood_usage_glucose_threshold`, `blood_to_liver_glucose_threshold`,
`igf_clearance`, `cl_to_p4_scale`, `cl_formation_scale`,
`insulin_secretion_max`, and `insulin_clearance`.

The result is local by design. It identifies influential mechanisms near the
nominal model, but it does not establish whether those parameters can be
estimated uniquely. Compensation and weak output directions are addressed by
the identifiability analyses.

In the sensitivity plot, blue bars indicate positive sensitivity and orange
bars indicate negative sensitivity. The sign identifies the direction of the
AUC response; the magnitude identifies the strength of the local response.

## 2. Biological Admissibility Filtering

Biological admissibility filtering restricts ensemble analyses to plausible ODE
trajectories. This matters because nonlinear endocrine-metabolic systems can
produce unrealistic oscillations, unstable biomarker levels, or failed cycle
structure even when sampled parameters remain inside a nominal prior range.

The historical admissibility panel used FSH, PGF, P4, E2, INH, IGF1, insulin,
and glucose. Glucagon was excluded from the admissibility rule but retained as
an observable biomarker for downstream uncertainty propagation, global
sensitivity, and Bayesian experimental design.

The enriched ODE-confirmed admissible ensemble contains 12,721 biologically
admissible simulations. Only ODE-confirmed rows are treated as scientific
truth. The admissible ensemble provides the basis for global sensitivity,
uncertainty propagation, and MI/BED ranking stability.

<img width="3520" height="1650" alt="image" src="https://github.com/user-attachments/assets/86fd9eaf-dced-48f7-86e1-6bda72f8be3e" />
<img width="3520" height="1650" alt="image" src="https://github.com/user-attachments/assets/9b6f6ad7-e22f-4178-bea6-ec20f9b14c0c" />
<img width="3520" height="1650" alt="image" src="https://github.com/user-attachments/assets/526e6f1d-12aa-410a-857d-b9b4c30b934b" />


## 3. Global Sensitivity Analysis

The global sensitivity heatmaps should be read as parameter-biomarker maps.
Each cell summarizes how one model parameter is associated with one observable
biomarker AUC across the 12,721-row ODE-confirmed admissible ensemble. PRCC
and Spearman provide complementary monotonic association summaries: PRCC
adjusts for the other sampled parameters, while Spearman shows the direct rank
association. Strong links indicate biomarkers that carry interpretable
information about a parameter within the biologically admissible model region.

<img width="2640" height="3763" alt="image" src="https://github.com/user-attachments/assets/c3531d43-e1a5-454e-98ca-83e103dc983f" />
<img width="2640" height="3763" alt="image" src="https://github.com/user-attachments/assets/682e01e0-1ebd-4fee-921e-a5e976a2293a" />
<img width="2420" height="1408" alt="image" src="https://github.com/user-attachments/assets/dd6c8f1e-992d-4ec7-b8e5-3450a103627c" />
<img width="2420" height="1430" alt="image" src="https://github.com/user-attachments/assets/927e6a94-8bfa-4949-8fc8-e1c9440ad1f0" />
<img width="2640" height="1560" alt="image" src="https://github.com/user-attachments/assets/0211ca59-edab-4bc7-acd2-55091a8191ca" />


Global sensitivity evaluates parameter-biomarker AUC associations across the
ODE-confirmed admissible ensemble. PRCC and Spearman heatmaps show which
parameters are monotonically associated with each observable biomarker.

Metabolic thresholds and transport parameters tend to associate with glucose,
insulin, IGF1, and glucagon AUC features. Reproductive endocrine parameters
link to FSH, PGF, P4, E2, and INH features. These links define biologically
interpretable parameter-biomarker relationships for downstream BED targeting.

The associations are not causal effects and are not formal Sobol variance
decompositions. They are monotonic screening signals inside the admissible
simulation ensemble.

## 4. SVD Identifiability

The SVD figures summarize whether the measured outputs can distinguish
parameter directions. The singular-value spectrum separates well-informed
directions from weak directions. The compensation network shows parameters
that can trade off inside poorly informed directions. The decision map then
translates these diagnostics into practical estimate, fix, or anchor
candidates for downstream inference.

![Singular values](../results_final/figures/identifiability_singular_values.png)

![Compensation network](../results_final/figures/identifiability_compensation_network_core.png)


<img width="3545" height="2005" alt="image" src="https://github.com/user-attachments/assets/1de6b7f8-37ee-4602-a24b-39816415d127" />


The SVD screen evaluates which local parameter directions are expressed in the
measured output space. Rapidly decaying singular values indicate weakly
informed directions. Nullspace participation identifies parameters that
contribute strongly to poorly informed directions.

The singular-value spectrum summarizes the information content of the local
trajectory sensitivity matrix. Large singular values correspond to
well-informed parameter directions, while very small singular values correspond
to weakly informed or practically non-identifiable directions. The spectrum
spans many orders of magnitude, indicating that the selected measurement panel
strongly informs some parameter combinations while leaving other combinations
weakly determined.

The interpretation is not that the entire model is identifiable or
non-identifiable. Instead, the model contains a mixture of well-informed and
poorly informed parameter directions. This is the main reason the analysis
separates "influential" from "identifiable": a parameter can move an output
locally while still being difficult to estimate if another mechanism can
produce a similar trajectory change.

The compensation network shows how parameters can trade off while preserving
similar output behavior. Biologically, this reflects the fact that multiple
mechanisms can influence the same endocrine-metabolic trajectories. This is why
local influence does not automatically imply practical identifiability.

The decision map combines two pieces of information: whether a parameter
affects the selected outputs, and whether it participates in weak or
compensatory directions. A parameter is a strong estimation candidate only if
it affects the outputs and is not strongly confounded with other parameters.

Using measurable outputs FSH, PGF, P4, E2, INH, IGF1, insulin, glucose, and
glucagon, the SVD screen classified the analyzed parameters into 60 estimate
candidates, 12 fixed-anchor candidates, and 25 fixed low-information
candidates. The estimate group contains parameters that are sufficiently
sensitive and have low enough nullspace participation to be reasonable
calibration or profile-likelihood candidates. Examples include
`insulin_glucose_threshold`, `blood_usage_glucose_threshold`,
`inhibin_clearance`, `blood_to_liver_glucose_threshold`,
`hp_p4_follicle_scale`, and `hm_ih_2_threshold`.

The fixed-anchor group contains parameters that remain influential but also
participate strongly in weak or compensatory directions. Examples include
`insulin_clearance`, `insulin_secretion_max`, `fsh_syn_scale`,
`insulin_fsh_threshold`, `blood_usage_max`, `glucagon_clearance`, and
`glucagon_secretion_max`. These parameters are biologically important, but the
current output panel does not separate them cleanly enough to estimate all of
them freely together.

The fixed low-information group should be interpreted only within this
analysis setting. It does not mean biologically irrelevant. It means that, for
the selected outputs and simulation window, the parameter has low local output
influence and should not be estimated from this dataset.

The all-parameter decision view separates four practical cases: sensitive and
weakly compensated parameters are the best estimate candidates; sensitive but
highly compensated parameters are important but risky to estimate freely; low
sensitivity parameters are poor calibration targets in the current experiment;
and compensatory but weakly sensitive parameters are poor candidates unless a
new experiment is designed to excite them.

## 5. Profile Likelihood

![Representative profile likelihood classes](../results_final/figures/profile_likelihood_representative_3x3.png)

Profile likelihood provides nonlinear confirmation of practical
identifiability. The representative 3 x 3 panel contains:

- Practically identifiable parameters:
  `insulin_glucose_threshold`, `inhibin_clearance`,
  `hp_p4_follicle_scale`.
- Boundary-limited parameters:
  `blood_to_liver_glucose_threshold`, `gnrh_clearance`,
  `hp_iof_threshold`.
- Weak or flat parameters:
  `insulin_igf_threshold`, `feed_direct_blood_fraction`,
  `lh_basal_release`.

Practically identifiable profiles show clear curvature around the optimum.
Boundary-limited profiles show partial or one-sided support. Weak/flat profiles
remain broad or shallow, indicating limited learnability under the available
observable panel.

This nonlinear check complements the SVD screen: SVD identifies weak local
directions and compensation structure, while profile likelihood tests whether
those weaknesses remain after nuisance parameters are allowed to re-adjust over
a wider parameter range.

The top row shows representative practically identifiable parameters
(`insulin_glucose_threshold`, `inhibin_clearance`, and
`hp_p4_follicle_scale`) with clear profile minima and likelihood increases on
both sides of the optimum. These parameters appear sufficiently constrained by
the current measurement panel.

The middle row illustrates boundary-limited cases
(`blood_to_liver_glucose_threshold`, `gnrh_clearance`, and
`hp_iof_threshold`). These profiles show local curvature but remain
incompletely bounded within the explored range. This indicates that the
parameters are locally informative but not fully constrained under the current
experimental setup.

The bottom row demonstrates weak or flat failure modes.
`insulin_igf_threshold` has some curvature but remains broad, whereas
`feed_direct_blood_fraction` and `lh_basal_release` show flat or weakly
bounded profiles. These cases demonstrate that local sensitivity or SVD
estimate classification does not automatically guarantee nonlinear practical
identifiability.

This classification becomes the reference diagnosis for BED and Bayesian
posterior comparisons.

## 6. Uncertainty Propagation

<img width="3300" height="2772" alt="image" src="https://github.com/user-attachments/assets/1e730211-c338-4336-bfa0-31134330b2c3" />


Uncertainty propagation summarizes the 5th-95th percentile range, ensemble
median, and nominal trajectory across admissible ODE-confirmed simulations.
The informative windows are periods where biologically plausible trajectories
separate meaningfully.

These windows are used by BED to restrict candidate sampling days. The
uncertainty results therefore connect admissible model variability to
experimental timing.

## 7. Bayesian Experimental Design

Bayesian experimental design is the constructive follow-up to sensitivity and
identifiability analysis. Sensitivity, SVD, compensation, and profile
likelihood diagnose which parameter directions are influential, weakly
informed, or confounded. BED then asks which future measurements are most
likely to reduce those weaknesses.

The main BED interpretation is based on the parameter-specific guided
posterior update. The independent observation, cumulative acquisition,
high-versus-low day, and MI day-curve plots are supporting traceability and
sanity-check outputs.

### MI Biomarker Ranking

![Guided biomarker MI bars](../results_final/figures/bed_guided_parameter_biomarker_mi_bars.png)

The biomarker MI bars identify which observable biomarkers are most informative
for each representative parameter. This decomposes the BED result by measured
species and explains why different parameters receive different guided
observation scenarios.

These bars connect GSA to BED. GSA first identifies biomarkers whose AUCs are
associated with each parameter across the admissible ensemble. MI then asks
which of those biomarkers would be most informative as measurements for
posterior updating. Biomarkers that rank highly in both analyses are strong
candidates for targeted experimental observation.

The interpretation is mechanistic. Metabolic threshold parameters are expected
to link most strongly to insulin, glucose, IGF1, and glucagon observations.
Reproductive endocrine parameters are expected to link most strongly to PGF,
FSH, P4, E2, and INH observations. A mismatch between these expectations and
the MI ranking would indicate either compensation, weak excitation, or a
biomarker that carries indirect information through coupled model dynamics.

### Guided Parameter-Specific Updates

![Guided posterior updates](../results_final/figures/bed_targeted_gsa_uncertainty_guided_posteriors_3x3.png)

The guided update is parameter-specific. For each representative parameter,
the workflow selects biomarkers from global sensitivity links, restricts days
to uncertainty-rich windows, ranks biomarker-day candidates by MI, and updates
the posterior using the selected scenario.

This creates the intended thesis-style chain:

```text
profile likelihood class
-> GSA-linked biomarkers
-> uncertainty-selected time windows
-> MI-ranked biomarker/day observations
-> posterior update
```

The posterior plots should therefore be interpreted relative to the profile
likelihood class. Practically identifiable parameters show the strongest
posterior contraction because the existing observable system already contains
information that can constrain them. Boundary-limited parameters show partial
or one-sided updates because the selected observations constrain only part of
the explored range. Weak or flat parameters remain broad when the chosen
biomarkers and times do not resolve compensating model directions.

This result is not a claim that BED can rescue every non-identifiable
parameter. It shows where additional measurements improve learning and where
the model remains weakly informed even after targeted observation selection.

### Supporting Traceability: Independent Observation Scenarios

![Independent observation posteriors](../results_final/figures/bed_targeted_independent_observation_posteriors_3x3.png)

Independent observation scenarios show how individual selected measurements
reshape representative parameter priors. These plots separate direct
single-scenario learning from cumulative acquisition effects.

The scientific purpose is to test whether a single biologically plausible
measurement can move the broad +/-5% prior in the direction expected from the
profile-likelihood diagnosis. Practically identifiable parameters are expected
to respond more strongly when the observation targets a linked biomarker.
Boundary-limited parameters may move mainly toward one side of the prior.
Weak/flat parameters may remain broad because the selected observation does
not isolate that parameter from compensating mechanisms.

### Supporting Traceability: Cumulative Biomarker Acquisition

![Cumulative biomarker posteriors](../results_final/figures/bed_targeted_cumulative_biomarker_posteriors_3x3.png)

The cumulative acquisition analysis adds biomarkers from best 1 through best 9
using a global biomarker ranking. It shows how posterior narrowing changes as
additional measurements are included. This is intentionally different from the
parameter-specific guided update.

This panel asks whether information accumulates as additional biomarkers are
measured. The expected thesis-style behavior is progressive contraction for
parameters whose profile likelihood already shows practical learnability.
For boundary-limited parameters, added biomarkers may improve one side of the
posterior more than the other. For weak/flat parameters, adding biomarkers can
change the posterior shape without producing a sharply bounded distribution,
which is consistent with limited identifiability.

### Supporting Traceability: High Versus Low Information Day

![High versus low information day](../results_final/figures/bed_targeted_high_vs_low_information_day_3x3.png)

The high-versus-low day comparison tests whether BED-ranked sampling days
produce stronger posterior narrowing than less informative days. The result
connects the MI ranking to posterior behavior.

This comparison links the time component of BED to posterior learning. A
high-information day should generally produce stronger or more interpretable
posterior movement than a low-information day for the same biomarker set. When
the difference is small, it indicates that either the parameter is weakly
learnable under that observation set or the selected biomarker is informative
over a broader window rather than at one sharply optimal day.

### Supporting Traceability: MI Day Curves

![Guided day MI curves](../results_final/figures/bed_guided_day_mi_curves.png)

The day curves show candidate-day information scores within
uncertainty-selected windows. Because each parameter has its own linked
biomarkers and uncertainty windows, the selected days can differ across
parameters.

These curves complete the chain from parameter class to observation timing.
Profile likelihood identifies the parameter class, GSA selects biologically
linked biomarkers, uncertainty propagation restricts attention to windows
where admissible trajectories separate, and MI ranks the remaining
biomarker-day candidates. The resulting posterior update is therefore based on
parameter-specific measurement logic rather than random observation selection.

The day-curve figure is a traceability plot rather than the main posterior
result. It explains why the guided posterior update uses particular
biomarker-day combinations.

## 8. Bayesian Inference

### Reduced ODE-Archive Posterior

![Reduced archive posterior](../results_final/figures/mcmc_prior_vs_posterior_3x3.png)

The reduced ODE-archive posterior uses likelihood weights over broad-prior ODE
archive rows. No new ODE simulations are proposed or accepted, and no
surrogate-generated trajectories are used as posterior truth. The resulting
posterior curves provide a transparent check on whether the selected
observations concentrate probability in parameter space.

This is an archive-based posterior approximation, not live ODE MCMC. The
posterior support is limited to the precomputed broad-prior simulation archive,
so roughness or multimodality can reflect both biological information and the
finite density of available archive rows.

### Archive ABC Filtering

![Archive ABC posteriors](../results_final/figures/abc_smc_parameter_specific_gsa_uncertainty_guided_posteriors_3x3.png)

Archive-based sequential ABC filtering reduces a distance tolerance over
precomputed ODE rows. It is not full adaptive ABC-SMC because particles are not
perturbed and ODEs are not rerun. The method cannot discover posterior regions
absent from the archive, but it provides a likelihood-free check using trusted
ODE simulations.

### Bayesian Method Comparison

![Bayesian method comparison](../results_final/figures/bayesian_method_comparison_summary.png)

The method comparison summarizes posterior behavior across reweighting,
reduced archive posterior updates, and archive-based sequential ABC filtering.
The plot compares posterior contraction behavior across archive-based
approximations; it is not intended to select a single "best" method. Agreement
across methods strengthens the interpretation that posterior narrowing is
driven by the selected observations rather than by one approximation scheme.
Differences across methods highlight sensitivity to likelihood weighting,
ABC tolerance, distance metric, and archive coverage.

The comparison should therefore be interpreted as a consistency check across
complementary archive-based approximations: posterior reweighting, reduced
ODE-archive posterior updating, and archive-based sequential ABC filtering.

## 9. Integrated Scientific Story

The results form a coherent scientific chain:

1. Local sensitivity identifies mechanisms that affect outputs near the nominal
   model.
2. SVD identifies weak directions and compensatory parameter combinations.
3. Profile likelihood confirms which representative parameters are practically
   identifiable, boundary-limited, or weak/flat.
4. Biological admissibility filtering restricts ensemble analyses to plausible
   ODE trajectories.
5. Global sensitivity links parameters to observable biomarkers.
6. Uncertainty propagation identifies informative sampling windows.
7. BED selects biomarker-day observations for parameter learning.
8. Bayesian posterior updates test whether those observations reduce
   uncertainty.

Bayesian updating confirms the identifiability diagnosis rather than
contradicting it. Practically identifiable parameters narrow strongly because
the selected observable panel contains information that constrains them.
Boundary-limited parameters show partial or one-sided updates because the
available observations constrain only part of the explored range. Weak or flat
parameters remain broad even under guided observations because the selected
measurements do not fully resolve compensating model directions.

The final conclusion is therefore constructive rather than negative. BED
improves learning where information exists, but it does not magically rescue
parameters that are practically non-identifiable under the available output
panel. Identifiability analysis diagnoses weak directions; BED identifies
measurement choices that can reduce those weaknesses when the model and
observable biomarkers contain sufficient information.
