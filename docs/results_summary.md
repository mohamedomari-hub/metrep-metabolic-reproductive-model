# Results Summary

This repository curates the Python sensitivity, identifiability, profile-likelihood, and surrogate Bayesian experimental design outputs for the MetRep metabolic-reproductive model.

The full biological scenario simulations are reported in the MetRep model paper, the dexamethasone perturbation scenario is reported in the Dexa paper, and the full Bayesian experimental design analysis is reported in the PhD thesis. The purpose of this repository is to show the reproducible Python workflow: model validation, local sensitivity screening, SVD-based identifiability analysis, compensation diagnostics, nonlinear profile-likelihood confirmation, and a surrogate BED pilot.

For technical definitions and formulas, see:

- docs/sensitivity_identifiability.md
- docs/plot_interpretation_guide.md
- docs/bayesian_experimental_design.md

---

## Result Set

| Step | Main question | Main output | Appendix / diagnostic output |
|---|---|---|---|
| Python validation | Does the Python translation preserve model metadata and simulate correctly? | results_final/tables/validation_report.csv, results_final/figures/baseline_selected_states.png | results_final/figures/baseline_all_states.png |
| Sensitivity | Which parameters most strongly affect selected outputs? | results_final/figures/sensitivity_top_parameters.png, results_final/tables/sensitivity_top_parameters.csv | — |
| SVD identifiability | Which parameter directions are informed or weakly informed? | results_final/figures/identifiability_singular_values.png | results_final/figures/identifiability_svd_ranking.png |
| Estimate/fix decision | Which parameters should be estimated, anchored, or fixed? | results_final/figures/identifiability_decision_map.png, results_final/tables/structid_50d_measurable_holistic_table.csv | results_final/figures/identifiability_class_counts.png |
| Compensation | Which parameters compensate each other? | results_final/figures/identifiability_compensation_network_core.png, results_final/figures/identifiability_parameter_scenario_map.png | results_final/figures/identifiability_compensation_edges.png, results_final/figures/identifiability_nullspace_participation_all_parameters.png, results_final/figures/identifiability_sensitivity_vs_nullspace_all_parameters.png, results_final/figures/identifiability_sensitivity_ranked_by_class.png |
| Profile likelihood | Do selected parameters remain identifiable under nonlinear profiling? | results_final/figures/profile_likelihood_representative_3x3.png, results_final/tables/profile_50d_balanced_relaxed_summary.csv | results_final/figures/profile_50d_balanced_synthetic_outputs.png |
| Surrogate BED pilot | Which sampling day/species panel is informative, and how much faster is the surrogate? | analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_thesis_style_summary.png | analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png |

---

# Main Results

## 1. Baseline Python Model Validation

![Baseline selected states](../results_final/figures/baseline_selected_states.png)

The baseline simulation is used as a functional check of the translated Python model before sensitivity, identifiability, profile likelihood, or BED analyses are interpreted. It is not itself an identifiability test. Its purpose is to confirm that the translated model produces biologically plausible trajectories for the selected metabolic and reproductive states.

Baseline selected states

The selected-state baseline panel should be used in the main report because it gives a readable overview of the key model outputs. The full 22-state panel is useful as a technical completeness check, but it is too dense for the main narrative and is therefore better placed in the appendix.

---

## 2. Local Sensitivity Screening

Top local sensitivity parameters

![Top local sensitivity parameters](../results_final/figures/sensitivity_top_parameters.png)

The local sensitivity screen ranks parameters by their effect on output AUCs. Instead of evaluating sensitivity at a single arbitrary time point, the AUC-based score summarizes the cumulative trajectory response over the simulation window. This is useful for an oscillatory endocrine-metabolic model because a single time point can be misleading if the system is in a follicular phase, luteal phase, transition period, or pulse event.

The strongest AUC-level responses are concentrated in glucose-insulin regulation, IGF dynamics, and reproductive steroid-related parameters. The highest-ranking parameters include insulin_glucose_threshold, blood_usage_glucose_threshold, blood_to_liver_glucose_threshold, igf_clearance, cl_to_p4_scale, cl_formation_scale, insulin_secretion_max, and insulin_clearance.

This result indicates that the selected outputs are responsive to both metabolic control parameters and reproductive endocrine parameters. However, sensitivity alone is not sufficient for parameter estimation. A parameter can strongly affect the model outputs and still be difficult to estimate if another parameter can compensate for it.

---

## 3. SVD Identifiability Spectrum

Singular values

![Singular values](../results_final/figures/identifiability_singular_values.png)

The singular-value spectrum summarizes the information content of the local trajectory sensitivity matrix. Large singular values correspond to well-informed parameter directions, while very small singular values correspond to weakly informed or practically non-identifiable directions.

The spectrum spans many orders of magnitude, with the largest singular value around 1.27e3 and the smallest values around 9.04e-14. This large numerical spread indicates that the selected measurement panel informs some combinations of parameters strongly, while leaving other combinations weakly determined.

The interpretation is not that the entire model is identifiable or non-identifiable. Instead, the model contains a mixture of well-informed and poorly informed parameter directions. This motivates the next step: identifying which individual parameters are sensitive, which participate in weak directions, and which should be estimated, anchored, or fixed.

---

## 4. Estimate, Anchor, and Irrelevant Parameter Classes

Decision map

![Decision map](../results_final/figures/identifiability_decision_map.png)

The decision map combines two pieces of information:

- Sensitivity: does the parameter affect the selected outputs?
- Nullspace involvement: does the parameter participate in weak or compensatory directions?

This distinction is essential because high sensitivity does not automatically mean good estimability. A parameter is a good estimation candidate only if it affects the outputs and is not strongly confounded with other parameters.

Using measurable outputs FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon, the SVD screen classified the analyzed parameters as:

- 60 Estimate
- 12 Fix (anchor)
- 25 Fix (irrelevant)

The Estimate group contains parameters that are sufficiently sensitive and have low enough nullspace participation to be reasonable calibration or profile-likelihood candidates. Examples include insulin_glucose_threshold, blood_usage_glucose_threshold, inhibin_clearance, blood_to_liver_glucose_threshold, hp_p4_follicle_scale, and hm_ih_2_threshold.

The Fix (anchor) group contains parameters that remain influential but also participate strongly in weak or compensatory directions. Examples include insulin_clearance, insulin_secretion_max, fsh_syn_scale, insulin_fsh_threshold, blood_usage_max, glucagon_clearance, and glucagon_secretion_max. These parameters are biologically important, but the current output panel does not separate them cleanly enough to estimate all of them freely together.

The Fix (irrelevant) group should be interpreted only within this analysis setting. It does not mean biologically irrelevant. It means that, for the selected outputs and simulation window, the parameter has low local output influence and should not be estimated from this dataset. Several lactation and storage threshold parameters fall into this group because they are weakly excited in the non-lactating baseline setting.

---

## 5. Core Compensation Network

Core compensation network

![Core compensation network](../results_final/figures/identifiability_compensation_network_core.png)

The compensation network visualizes the strongest parameter trade-offs extracted from weak SVD directions. Each node is a parameter, and each edge indicates that two parameters can compensate for each other in the local sensitivity structure.

Node color shows the estimate/fix recommendation:

- blue = estimate candidate
- orange = fixed anchor
- gray = fixed irrelevant

Node size shows sensitivity, and edge width shows compensation strength. Large, highly connected nodes are especially important because they are influential but also involved in parameter trade-offs.

Strong compensation relationships include pairs such as:

- fsh_syn_scale with insulin_fsh_scale
- blood_usage_max with insulin_clearance
- insulin_clearance with insulin_fsh_threshold
- hp_pg_iof_scale with hp_enz_pg_scale
- milk_fat_threshold with milk_storage_threshold
- lactation_oxt_decay with fat_mobilization_storage_threshold

These pairings indicate directions where changes in one parameter can be partly offset by changes in another. Such parameters should not be estimated together without additional information, stronger priors, anchoring, or a more informative experimental design.

---

## 6. Parameter Scenario Map

Parameter scenario map

![Parameter scenario map](../results_final/figures/identifiability_parameter_scenario_map.png)

The parameter scenario map is the main all-parameter diagnostic. It replaces the unreadable all-parameter network because it shows all parameters without forcing every parameter into a crowded graph.

Each point is one analyzed parameter. The x-axis shows nullspace or compensation involvement, and the y-axis shows local sensitivity. This separates four practical cases:

1. Sensitive and weakly compensated  
   These are the best estimate candidates.

2. Sensitive but highly compensated  
   These parameters are important but risky to estimate freely because other parameters can mimic their effects.

3. Low sensitivity  
   These parameters are poor calibration targets in the current experiment and should usually remain fixed.

4. Compensatory but weakly sensitive  
   These parameters participate in weak directions but do not strongly affect the selected outputs. They are poor candidates for estimation unless a new experiment is designed to excite them.

This figure is especially useful because it summarizes the practical modelling decision: estimate, fix as anchor, or fix as low-information.

---

## 7. Profile Likelihood Confirmation

Representative profile likelihood classes

![Representative profile likelihood classes](../results_final/figures/profile_likelihood_representative_3x3.png)

Profile likelihood was used as a nonlinear confirmation step after the local SVD screen. The SVD analysis is local and linear, while profile likelihood checks whether selected parameters remain practically identifiable when the model is refitted over a parameter range.

The 60 selected Estimate parameters were profiled against synthetic measurable outputs while nuisance parameters were allowed to compensate. The profile-likelihood classes were:

- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable
- 2 flat/non-identifiable

The 3x3 panel shows representative examples from each class instead of plotting all 60 profiles in one crowded figure.

The top row shows practically identifiable examples with clear profile minima. These support the SVD decision that a substantial subset of selected parameters is estimable from the current output panel.

The middle row shows boundary-limited examples. Across the full table, the boundary-limited parameters were blood_to_liver_glucose_threshold, igf_lh_sensitivity_threshold, hp_enz_pg_threshold, gnrh_clearance, hp_p4_enz_scale, and hp_iof_threshold. Their best profile point occurred at the edge of the tested range, so they should not be claimed as fully bounded without wider profile ranges or additional information.

The bottom row shows weakly identifiable and flat/non-identifiable cases. The weakly identifiable parameter was insulin_igf_threshold. Its profile had some curvature, but a broad part of the tested range remained acceptable. The flat/non-identifiable parameters were feed_direct_blood_fraction and lh_basal_release; their tested changes did not increase the loss enough to support reliable estimation under the current output panel.

---

## 8. Bayesian Experimental Design Link

Sensitivity, SVD, compensation analysis, and profile likelihood identify where the current measurement panel is informative and where it leaves uncertainty. Bayesian experimental design is the constructive follow-up: it asks which future sampling days and measured species are expected to reduce uncertainty most efficiently.

Surrogate BED thesis-style summary

![Surrogate BED thesis-style summary](../analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_thesis_style_summary.png)

![ODE versus surrogate speed](../analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png)

The surrogate BED pilot uses the translated Python ODE model to generate a local Monte Carlo prior around the nominal parameter set, retains biologically admissible trajectories, trains an ExtraTrees surrogate on ODE output features, and evaluates candidate measurements using mutual information and posterior reweighting.

In this pilot run, 389 biologically admissible samples remained from 500 ODE simulations. For the target parameter insulin_glucose_threshold, the highest all-species candidate was day 88, while the lowest all-species candidate was day 56. The top trajectory panel places these candidate days in the nominal follicle/P4 cycle. The middle panel ranks candidate sampling days by estimated mutual information. The lower-left panel decomposes the day-88 information by species. The lower-right panel shows how different measurement choices reshape the prior into posterior densities by likelihood-weighted reweighting.

This BED result should be labelled as a surrogate pilot, not as the final thesis-level BED result. It is useful for reproducing the structure of the PhD analysis in Python, but a final public BED claim should use more ODE training samples and demonstrate stable surrogate validation, mutual-information convergence, and candidate-ranking robustness.

The scientific workflow is therefore sequential:

1. Sensitivity identifies parameters that affect outputs.
2. SVD identifies informed and weak directions.
3. Compensation analysis explains which parameters are confounded.
4. Profile likelihood confirms nonlinear identifiability behaviour.
5. BED proposes future measurements to reduce uncertainty.

---

# Appendix / Diagnostic Figures

The following figures are useful for technical checking and transparency, but they are not needed in the main narrative.

## Appendix A. Full Baseline State Panel

Baseline all model states

![Baseline all model states](../results_final/figures/baseline_all_states.png)

The all-state baseline panel shows every model state in the translated 22-state Python implementation. This figure is useful for confirming that no state is silently missing, unstable, or incorrectly mapped. It is placed in the appendix because it is too dense for the main results.

---

## Appendix B. SVD Sensitivity Ranking

SVD ranking

![SVD ranking](../results_final/figures/identifiability_svd_ranking.png)

The SVD ranking summarizes which parameters contribute most strongly to the local trajectory sensitivity matrix after scaling by parameter and output magnitude. This ranking is related to, but not identical to, the AUC-based local sensitivity ranking.

The AUC sensitivity ranking asks: which parameters change cumulative output exposure?

The SVD ranking asks: which parameters change the full time-resolved trajectory sensitivity matrix?

Therefore, parameters that change pulse timing, oscillation shape, or trajectory patterns can rank highly in the SVD score even if their AUC effect is less dominant. The highest-ranked parameters should not automatically be selected for estimation; they must also be checked for nullspace participation and compensation.

---

## Appendix C. Identifiability Class Counts

Class counts

![Class counts](../results_final/figures/identifiability_class_counts.png)

This figure summarizes the number of parameters in each recommendation class. It is useful as a quick overview, but the same information is already stated in the main text. Therefore, it is better treated as an appendix figure rather than a main result.

---

## Appendix D. Full Nullspace Participation Ranking

All-parameter nullspace participation

![All-parameter nullspace participation](../results_final/figures/identifiability_nullspace_participation_all_parameters.png)

This plot ranks all parameters by participation in weak or non-identifiable SVD directions. High values indicate parameters that appear strongly in poorly informed directions.

The figure is useful for model debugging and detailed parameter review, but it contains many parameter names and is visually dense. The parameter scenario map in the main results summarizes the same concept more clearly.

---

## Appendix E. Full Sensitivity versus Nullspace Scatter

Sensitivity versus nullspace

![Sensitivity versus nullspace](../results_final/figures/identifiability_sensitivity_vs_nullspace_all_parameters.png)

This plot shows all parameters in sensitivity-nullspace space. It is useful as a raw diagnostic, but the improved parameter scenario map is clearer and easier to interpret. This version is therefore kept as a supplementary diagnostic.

---

## Appendix F. Strongest Compensation Pairs

Compensation pairs

![Compensation pairs](../results_final/figures/identifiability_compensation_edges.png)

This barplot ranks the strongest compensation pairs extracted from weak SVD directions. It provides a quantitative counterpart to the compensation network. Since the network is more intuitive for the main narrative, this pair-ranking plot is kept in the appendix.

---

## Appendix G. Sensitivity Ranked by Class

Sensitivity ranked by class

![Sensitivity ranked by class](../results_final/figures/identifiability_sensitivity_ranked_by_class.png)

This plot shows all analyzed parameters ranked by sensitivity and colored by recommendation class. It is useful for checking whether a fixed parameter was fixed because it is low-sensitivity or because it is sensitive but confounded. It is placed in the appendix because the full list of parameters is too dense for the main results.

---

## Appendix H. Synthetic Profile-Likelihood Outputs

Synthetic measurable outputs

![Synthetic measurable outputs](../results_final/figures/profile_50d_balanced_synthetic_outputs.png)

This figure shows the synthetic measurable outputs used for profile likelihood. It supports transparency of the profiling setup but is methodological rather than interpretive, so it is better placed in the appendix.

---

## Appendix I. ODE versus Surrogate Speed

ODE versus surrogate speed

![ODE versus surrogate speed](../analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png)

The speed comparison was measured locally using 12 ODE simulations and 200 repeated surrogate predictions on the same number of parameter samples. The ODE solver required 0.327 seconds per sample, while surrogate prediction required 0.00152 seconds per sample, giving an observed prediction speedup of about 215x.

This figure is useful for the computational engineering story because it explains why a surrogate is useful for BED candidate exploration. For a biological modelling paper, it can remain in the appendix. For a GitHub or AI-engineering portfolio, it can also be shown in the main text.
