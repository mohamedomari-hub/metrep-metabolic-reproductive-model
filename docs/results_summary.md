# Results Summary

This repository curates the Python sensitivity and identifiability outputs for
the MetRep metabolic-reproductive model. The full scenario simulations are
reported in the MetRep model paper, the Dexa perturbation scenario is reported
in the Dexa paper, and the full Bayesian experimental design analysis is
reported in the PhD thesis.

For technical definitions and formulas, see:

- `docs/sensitivity_identifiability.md`
- `docs/plot_interpretation_guide.md`
- `docs/bayesian_experimental_design.md`

## Result Set

| Step | Main question | Curated output |
|---|---|---|
| Python validation | Does the Python translation preserve MATLAB metadata? | `results_final/tables/validation_report.csv` |
| Baseline run | Does the translated model produce baseline trajectories? | `results_final/figures/baseline_selected_states.png`, `results_final/figures/baseline_all_states.png` |
| Sensitivity | Which parameters most strongly affect selected outputs? | `results_final/tables/sensitivity_top_parameters.csv`, `results_final/figures/sensitivity_top_parameters.png` |
| SVD identifiability | Which parameter directions are informed or weak? | `results_final/figures/identifiability_singular_values.png`, `results_final/figures/identifiability_svd_ranking.png` |
| SVD classes | Which parameters should be estimated, anchored, or fixed? | `results_final/tables/structid_50d_measurable_holistic_table.csv`, `results_final/figures/identifiability_decision_map.png`, `results_final/figures/identifiability_class_counts.png` |
| Compensation | Which parameters can compensate each other? | `results_final/tables/structid_50d_measurable_compensation_edges.csv`, `results_final/figures/identifiability_compensation_edges.png`, `results_final/figures/identifiability_compensation_network_core.png`, `results_final/figures/identifiability_parameter_scenario_map.png` |
| Profile likelihood | Which selected parameters remain practically identifiable after nonlinear profiling? | `results_final/tables/profile_50d_balanced_relaxed_summary.csv`, `results_final/figures/profile_likelihood_representative_3x3.png` |
| Surrogate BED pilot | Which sampling day/species panel is expected to be informative, and how much faster is the surrogate than ODE simulation? | `analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_thesis_style_summary.png`, `analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png` |

## Baseline And Sensitivity

The baseline simulation is used as a functional check of the translated Python
model before analysis. It is not itself an identifiability test; it establishes
that the model produces trajectories for the selected metabolic and
reproductive states.

![Baseline all model states](../results_final/figures/baseline_all_states.png)

The all-state baseline panel shows every model state in the translated
22-state Python implementation. This view is useful for checking that no state
is silently missing or numerically unstable before sensitivity,
identifiability, profile likelihood, or BED analyses are interpreted.

![Top local sensitivity parameters](../results_final/figures/sensitivity_top_parameters.png)

The local sensitivity screen shows that the strongest AUC-level responses are
concentrated in insulin/glucose regulation, IGF dynamics, and reproductive
steroid-related parameters. The largest absolute relative sensitivity in the
curated table is `insulin_glucose_threshold` acting on insulin AUC
(`4.43`), followed by `blood_usage_glucose_threshold` acting on insulin AUC
(`3.08`). Other high-ranking parameters include
`blood_to_liver_glucose_threshold`, `igf_clearance`, `cl_to_p4_scale`,
`cl_formation_scale`, `insulin_secretion_max`, and `insulin_clearance`.

This result indicates that the selected outputs are responsive to both
metabolic control parameters and reproductive endocrine parameters. However,
sensitivity alone is not sufficient for estimation. A parameter can produce a
large output response and still be difficult to estimate if another parameter
can compensate for it.

## SVD Identifiability

![SVD ranking](../results_final/figures/identifiability_svd_ranking.png)

The SVD ranking summarizes which parameters contribute most strongly to the
local trajectory sensitivity matrix after scaling by parameter and output
magnitude. The highest-ranked parameters are not automatically selected for
estimation; they must also be checked for nullspace participation.

![Singular values](../results_final/figures/identifiability_singular_values.png)

The singular-value spectrum shows a large separation between the strongest
informed directions and the weakest directions. In the curated singular-value
table, the largest singular value is about `1.27e3`, while the smallest values
fall to about `9.04e-14`. This numerical spread supports the conclusion that
the current measurement set informs some parameter combinations strongly but
leaves other combinations weakly determined.

The biological interpretation is not that the whole model is identifiable or
non-identifiable. Instead, the selected outputs provide information about some
directions in parameter space, while other directions remain ambiguous under
the local linear approximation.

## Estimate, Anchor, And Irrelevant Classes

![Decision map](../results_final/figures/identifiability_decision_map.png)

![Class counts](../results_final/figures/identifiability_class_counts.png)

Using measurable outputs `FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose,
Glucagon`, the SVD screen classified the analyzed parameters as:

- 60 `Estimate`
- 12 `Fix (anchor)`
- 25 `Fix (irrelevant)`

The `Estimate` group contains parameters that are sensitive enough and have low
enough nullspace participation to be reasonable calibration/profile-likelihood
candidates. Examples near the top of the holistic table include
`insulin_glucose_threshold`, `blood_usage_glucose_threshold`,
`inhibin_clearance`, `blood_to_liver_glucose_threshold`,
`hp_p4_follicle_scale`, and `hm_ih_2_threshold`.

The `Fix (anchor)` group contains parameters that retain output influence but
also participate strongly in weak or compensatory directions. Examples include
`insulin_clearance`, `insulin_secretion_max`, `fsh_syn_scale`,
`insulin_fsh_threshold`, `blood_usage_max`, `glucagon_clearance`, and
`glucagon_secretion_max`. These parameters are scientifically important, but
the current output panel does not separate them cleanly enough to estimate all
of them freely together.

The `Fix (irrelevant)` label should be interpreted only within this analysis
setting. It does not mean biologically irrelevant. It means that, for the
selected outputs and simulation window, the parameter has low local output
influence and therefore should not be estimated from this dataset. Several
lactation/storage threshold parameters fall into this group because their
sensitivity is near zero in this non-lactating baseline setting while their
nullspace participation is high.

## Nullspace And Compensation

![All-parameter nullspace participation](../results_final/figures/identifiability_nullspace_participation_all_parameters.png)

![Sensitivity versus nullspace](../results_final/figures/identifiability_sensitivity_vs_nullspace_all_parameters.png)

The all-parameter nullspace plot shows which parameters participate most in
weak directions. The sensitivity/nullspace map then separates two cases that
need different scientific treatment:

- sensitive parameters with low nullspace participation are supported by the
  current outputs and are better candidates for estimation;
- sensitive parameters with high nullspace participation are influential but
  confounded, so they are better treated as anchors or constrained by prior
  knowledge;
- low-sensitivity parameters with high nullspace participation are weakly
  informed under the current output panel and should remain fixed unless a new
  experiment is designed to target them.

![Compensation pairs](../results_final/figures/identifiability_compensation_edges.png)

![Core compensation network](../results_final/figures/identifiability_compensation_network_core.png)

The compensation analysis identifies pairs of parameters that appear together
in weak SVD directions. Strong pairings in the curated table include:

- `lactation_oxt_decay` with `fat_mobilization_storage_threshold`
- `milk_fat_threshold` with `milk_storage_threshold`
- `fsh_syn_scale` with `insulin_fsh_scale`
- `hp_pg_iof_scale` with `hp_enz_pg_scale`
- `blood_usage_max` with `insulin_clearance`
- `insulin_clearance` with `insulin_fsh_threshold`

These pairings indicate directions where changes in one parameter can be
partly offset by changes in another. The core network keeps only the strongest
compensation edges so the parameter trade-offs remain readable. Node color
shows the estimate/fix recommendation, node size shows sensitivity, and edge
width shows compensation strength. Parameters that are both large and strongly
connected are the most important targets for anchoring, prior constraints, or
Bayesian experimental design.

![Parameter scenario map](../results_final/figures/identifiability_parameter_scenario_map.png)

The scenario map replaces the unreadable all-parameter network as the main
all-parameter diagnostic. Each point is one analyzed parameter. The x-axis
shows nullspace or compensation involvement, and the y-axis shows local
sensitivity. This makes the four practical cases visible at once:

- sensitive and weakly compensated parameters are the best estimate candidates;
- sensitive but highly compensated parameters are risky to estimate freely;
- low-sensitivity parameters are better fixed in this experiment;
- compensatory but weakly sensitive parameters are poor calibration targets
  unless a new experiment is designed to excite them.

![Sensitivity ranked by class](../results_final/figures/identifiability_sensitivity_ranked_by_class.png)

The ranked sensitivity-by-class plot shows all analyzed parameters in one
ordered list. It is useful for checking whether a fixed parameter is fixed
because it is genuinely low-sensitivity in this analysis setting or because it
is sensitive but confounded and therefore better treated as an anchor.

## Profile Likelihood Confirmation

![Synthetic measurable outputs](../results_final/figures/profile_50d_balanced_synthetic_outputs.png)

![Representative profile likelihood classes](../results_final/figures/profile_likelihood_representative_3x3.png)

Profile likelihood was used as a nonlinear confirmation step after the local
SVD screen. The 60 selected `Estimate` parameters were profiled against
synthetic measurable outputs while nuisance parameters were allowed to
compensate.

The profile-likelihood classes were:

- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable
- 2 flat/non-identifiable

The 3x3 panel shows representative examples from each class rather than all 60
profiles in one crowded figure. The top row shows practically identifiable
examples with clear profile minima. These examples support the SVD decision
that a substantial subset of the selected parameters is estimable from the
current output panel.

The middle row shows boundary-limited examples. Across the full table, the 6
boundary-limited parameters were `blood_to_liver_glucose_threshold`,
`igf_lh_sensitivity_threshold`, `hp_enz_pg_threshold`, `gnrh_clearance`,
`hp_p4_enz_scale`, and `hp_iof_threshold`. Their best profile point occurred at
the edge of the tested range, so they should not be claimed as fully bounded
without either a wider profile range or additional information.

The bottom row shows the weakly identifiable and flat/non-identifiable cases.
The weakly identifiable parameter was `insulin_igf_threshold`. Its profile had
some curvature, but a broad part of the tested range remained acceptable. The
flat/non-identifiable parameters were `feed_direct_blood_fraction` and
`lh_basal_release`; their tested changes did not increase the loss enough to
support reliable estimation under the current output panel.

## Bayesian Experimental Design Link

BED should be read as the planned constructive follow-up to the
identifiability analysis. Sensitivity, SVD, compensation, and profile
likelihood identify which parameters or directions are insufficiently
informed; BED evaluates which candidate sampling days and measured species are
expected to reduce that uncertainty.

![Surrogate BED thesis-style summary](../analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_thesis_style_summary.png)

The surrogate BED pilot uses the translated Python ODE model to generate a
local Monte Carlo prior around the nominal parameter set, retains biologically
admissible trajectories, trains an ExtraTrees surrogate on the ODE output
features, and then evaluates candidate measurements using mutual information
and posterior reweighting.

In this pilot run, 389 biologically admissible samples remained from 500 ODE
simulations. For the target parameter `insulin_glucose_threshold`, the highest
all-species candidate was day 88, while the lowest all-species candidate was
day 56. The top trajectory panel places these candidate days in the nominal
follicle/P4 cycle. The middle panel ranks candidate sampling days by estimated
mutual information. The lower-left panel decomposes the day-88 information by
species; insulin contributed the largest single-species MI in this pilot. The
lower-right panel shows how different measurement choices reshape the prior
into posterior densities by likelihood-weighted reweighting of the prior
samples.

![ODE versus surrogate speed](../analyses/bayesian_experimental_design/surrogate_bed/run_outputs/figures/surrogate_bed_ode_vs_surrogate_speed.png)

The speed comparison was measured locally using 12 ODE simulations and 200
repeated surrogate predictions on the same number of parameter samples. The ODE
solver required `0.327` seconds per sample, while surrogate prediction required
`0.00152` seconds per sample, giving an observed prediction speedup of about
`215x`. This is the practical reason for using the surrogate: candidate
ranking, posterior curves, and convergence checks can be explored much more
quickly after the ODE training library has been generated.

This BED result should still be labeled as a surrogate pilot, not as the final
thesis result. The current pilot is useful for reproducing the structure of the
PhD analysis in Python, but a final public BED claim should use more ODE
training samples and demonstrate stable surrogate validation, mutual
information convergence, and candidate ranking.

Scientifically, this makes the workflow sequential: the classical analysis
defines the information gap, and BED proposes how future experiments could
reduce that gap.

## Publication-Reported Results

- MetRep scenario simulations: reported in the MetRep model paper.
- Dexa perturbation simulation/validation: reported in the Dexa paper.
- Bayesian experimental design: reported in the PhD thesis; the Python
  surrogate BED figures here are pilot reproductions of the thesis-style
  workflow and should be interpreted with the validation caveat above.
