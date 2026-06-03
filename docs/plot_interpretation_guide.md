# Plot Interpretation Guide

This guide explains the curated sensitivity and identifiability plots in
`results_final/figures/`. The goal is to make the figures understandable without
requiring the reader to inspect the full analysis code first.

## Baseline Simulation

Figure:

- `baseline_selected_states.png`

What it shows:

- A standard Python simulation of selected MetRep state variables.
- This is a basic functional check that the translated Python ODE model runs
  and produces trajectories for important metabolic and reproductive states.

Technical terms:

- **ODE model**: ordinary differential equation model; a system of equations
  describing how state variables change over time.
- **State variable**: a dynamic quantity in the model, such as glucose,
  insulin, IGF-1, P4, E2, follicle, or corpus luteum.
- **Baseline simulation**: a model run under the standard/reference scenario,
  before perturbing parameters or changing experimental design.

Interpretation:

- This figure is not an identifiability result. It shows that the Python model
  can produce the expected type of trajectories before analysis begins.

## Sensitivity: Top Parameters

Figure:

- `sensitivity_top_parameters.png`

Table:

- `sensitivity_top_parameters.csv`

What it shows:

- Parameters ranked by their largest absolute local sensitivity across selected
  model outputs.
- Higher bars indicate that a small parameter perturbation causes a larger
  relative change in at least one output metric.
- The x-axis is the maximum absolute relative sensitivity across outputs. In
  simple terms, it asks: "when this parameter is slightly changed, what is the
  largest relative response seen in the measured model outputs?"

Technical terms:

- **Sensitivity analysis**: tests how much model outputs change when parameters
  are perturbed.
- **Local sensitivity**: sensitivity near the nominal parameter value, not over
  the entire possible parameter range.
- **Relative sensitivity**: output change scaled relative to the baseline
  output and parameter perturbation.
- **AUC**: area under the curve; here it summarizes the total trajectory
  magnitude over time for a model output.
- **Absolute sensitivity**: ignores the direction of change and focuses on
  magnitude.

Interpretation:

- Sensitive parameters are influential, but sensitivity alone does not mean a
  parameter is identifiable.
- A parameter can strongly affect outputs and still be hard to estimate if
  another parameter can compensate for it.
- For example, a parameter may strongly change glucose or insulin trajectories,
  but if another parameter can be adjusted to undo that change, the data may not
  be able to separate the two parameters.
- This plot motivates the next step: identifiability analysis.

## SVD: Singular Values

Figure:

- `identifiability_singular_values.png`

Table:

- `structid_50d_measurable_singular_values.csv`

What it shows:

- The singular values of the local sensitivity matrix, plotted on a log scale.
- Large singular values correspond to strongly informed parameter directions.
- Very small singular values correspond to weakly informed directions.

Technical terms:

- **Sensitivity matrix**: a matrix whose entries describe how outputs change
  when parameters change.
- **SVD**: singular value decomposition; a linear algebra method that separates
  the sensitivity matrix into independent directions of information.
- **Singular value**: a number measuring how strongly the outputs respond along
  one independent parameter direction.
- **Log scale**: a plot scale where equal vertical distances represent
  multiplicative changes; useful when values span many orders of magnitude.
- **Weak direction**: a parameter combination that changes the outputs very
  little and is therefore hard to estimate from the chosen measurements.

Interpretation:

- A steep drop from large to tiny singular values indicates that the current
  measurement set informs some parameter combinations much better than others.
- This supports the idea that the model has identifiable and weakly
  identifiable directions, rather than all parameters being equally estimable.
- This plot is not a parameter ranking. It does not say which parameter is weak
  by itself. Instead, it says how many independent parameter combinations are
  well supported by the selected outputs.
- The large singular values represent combinations that the data can "see".
  The small tail represents combinations that the data can barely distinguish.
- Those weak combinations are then inspected in the nullspace and compensation
  plots.

## SVD: Ranking

Figure:

- `identifiability_svd_ranking.png`

Table:

- `structid_50d_measurable_ranking_and_participation.csv`

What it shows:

- Parameters ranked by their local SVD sensitivity-ranking score.
- These are the parameters with the largest contribution to output variation in
  the SVD screen.
- The score is normalized, so the highest ranked parameter is shown as 1 and
  the others are shown relative to it.

Technical terms:

- **Ranking score**: a scalar score summarizing how strongly each parameter
  contributes to the sensitivity matrix.
- **Output variation**: change in measurable model outputs caused by parameter
  perturbation.
- **Measurable outputs**: the selected outputs assumed to be experimentally
  observable: FSH, PGF, P4, E2, INH, IGF1, insulin, glucose, and glucagon.

Interpretation:

- Parameters high in this plot are important for the selected outputs.
- This plot should be read together with nullspace participation. A parameter
  can be high-impact but still compensatory.
- If a parameter is high in this ranking and low in nullspace participation, it
  is a good estimation candidate.
- If a parameter is high in this ranking and also high in nullspace
  participation, it is influential but difficult to estimate alone.

## SVD: Nullspace Participation

Figure:

- `identifiability_nullspace_participation.png`
- `identifiability_nullspace_participation_all_parameters.png`
- `identifiability_sensitivity_vs_nullspace_all_parameters.png`

Table:

- `structid_50d_measurable_ranking_and_participation.csv`
- `structid_50d_measurable_holistic_table.csv`

What it shows:

- Parameters ranked by how strongly they participate in weak or nullspace
  directions.
- High values indicate that the parameter is involved in combinations that are
  difficult to distinguish from the selected outputs.
- `identifiability_nullspace_participation.png` focuses on the strongest
  weak-direction participants.
- `identifiability_nullspace_participation_all_parameters.png` shows the full
  parameter list, so the reader can see both problematic and well-supported
  parameters.
- `identifiability_sensitivity_vs_nullspace_all_parameters.png` puts every
  parameter on one map: sensitivity on the x-axis and nullspace participation
  on the y-axis.

Technical terms:

- **Nullspace**: parameter directions that produce little or no change in the
  outputs under the local linear approximation.
- **Nullspace participation**: how much a parameter contributes to those weak
  directions.
- **Compensatory parameter**: a parameter whose effect can be offset by changes
  in another parameter, making separate estimation difficult.
- **Local linear approximation**: the approximation that small parameter changes
  affect outputs approximately linearly near the nominal parameter set.
- **Sensitivity axis**: how strongly a parameter affects outputs.
- **Nullspace axis**: how strongly a parameter belongs to weak or compensatory
  directions.

Interpretation:

- High nullspace participation means the parameter may be hard to estimate
  independently.
- These parameters are good candidates for fixing, anchoring, or targeting with
  improved experimental design.
- A parameter can be sensitive and still appear in the nullspace. This means it
  affects the model, but its effect is not unique enough under the current
  measurements.
- In the all-parameter sensitivity/nullspace map, high sensitivity and low
  nullspace means a favorable estimation candidate.
- High sensitivity and high nullspace means an important but compensatory
  parameter, usually better treated as an anchor or constrained by prior
  knowledge.
- Low sensitivity and high nullspace means the current outputs do not provide
  useful information for estimating that parameter.
- Low sensitivity and low nullspace means the parameter is not driving this
  analysis strongly.
- This is the clearest answer to the question "are they sensitive?" Some of the
  weak/nullspace parameters are sensitive, but others are weak because the
  selected outputs barely respond to them.

## SVD: Compensation Edges

Figure:

- `identifiability_compensation_edges.png`
- `identifiability_compensation_network.png`
- `identifiability_compensation_network_sensitivity.png`

Table:

- `structid_50d_measurable_compensation_edges.csv`
- `structid_50d_measurable_holistic_table.csv`

What it shows:

- Parameter pairs that strongly compensate each other in weak/nullspace
  directions.
- A strong edge means the two parameters appear together in a compensatory
  direction.
- The network plot shows the same idea as a graph: each node is a parameter,
  and each edge is a compensation relationship.
- `identifiability_compensation_network.png` emphasizes the compensation
  structure and SVD class.
- `identifiability_compensation_network_sensitivity.png` adds the sensitivity
  information: larger nodes are more sensitive, and red-outlined nodes are
  locally sensitive parameters.

Technical terms:

- **Compensation edge**: a pairwise relationship indicating that increasing one
  parameter while decreasing another can preserve similar model outputs.
- **Edge strength**: a score representing how strongly that pair appears in the
  nullspace structure.
- **Parameter pair**: two model parameters considered together.
- **Node**: one model parameter in the compensation graph.
- **Node color**: the SVD class of the parameter: `Estimate`,
  `Fix (anchor)`, or `Fix (irrelevant)`.
- **Node size**: in the sensitivity-aware network, scaled by local sensitivity,
  so larger nodes affect the selected outputs more strongly.
- **Red node outline**: a sensitive parameter in the compensation network
  using a normalized local sensitivity threshold.
- **Edge width**: scaled by compensation strength; thicker edges indicate
  stronger compensatory relationships.

Interpretation:

- This plot helps explain why some parameters should not all be freely
  estimated at once.
- Strong compensation pairs suggest where additional measurements or prior
  constraints may be needed.
- Clusters of connected nodes indicate groups of parameters that may move
  together in weakly informed directions.
- A node connected to many strong edges is a warning sign: that parameter may
  be influential, but its separate value is difficult to recover unless one of
  its partners is fixed, constrained, or targeted by a better experiment.
- In the sensitivity-aware network, the most important cases are large nodes
  with red outlines that are also connected by thick edges. These parameters
  matter to the outputs, but their effects can be masked by compensation.
- Small nodes without a red outline are less urgent for calibration in the
  current measurement setting, even if they appear in the compensation network.
- The graph is also a BED guide. Candidate sampling times/species should be
  judged by whether they reduce ambiguity within these connected compensation
  groups.

## SVD: Decision Map

Figure:

- `identifiability_decision_map.png`

Table:

- `structid_50d_measurable_holistic_table.csv`

What it shows:

- Each parameter positioned by normalized sensitivity and normalized nullspace
  involvement.
- Points are colored by the three-class decision: `Estimate`, `Fix (anchor)`,
  or `Fix (irrelevant)`.

Technical terms:

- **Normalized sensitivity**: sensitivity rescaled to a 0..1 range for
  comparison across parameters.
- **Normalized nullspace involvement**: nullspace participation rescaled to a
  0..1 range.
- **Estimate**: a parameter selected as suitable for estimation/profile
  likelihood.
- **Fix (anchor)**: a high-impact but compensatory parameter that should be
  fixed or constrained to stabilize estimation.
- **Fix (irrelevant)**: a parameter with low impact under the selected outputs,
  so estimating it is not useful in this analysis setting.

Interpretation:

- Parameters with high sensitivity and low nullspace involvement are the most
  favorable estimation candidates.
- Parameters with high sensitivity and high nullspace involvement are important
  but compensatory, so they are better treated as anchors.
- Parameters with low sensitivity are fixed as irrelevant for this selected
  measurement setting.
- `Fix (irrelevant)` does not mean biologically irrelevant. It means the
  current measurable outputs do not support estimating that parameter in this
  analysis.
- The decision map gives the operational rule: estimate the separable
  parameters, anchor the high-impact compensatory parameters, and leave
  low-information parameters fixed unless BED or a new dataset targets them.

## SVD: Class Counts

Figure:

- `identifiability_class_counts.png`

Table:

- `identifiability_class_counts.csv`

What it shows:

- The number of parameters assigned to each SVD class.

Technical terms:

- **Class count**: the number of parameters in each recommendation category.
- **Recommendation category**: the analysis decision assigned to a parameter:
  estimate, fix as anchor, or fix as irrelevant.

Interpretation:

- This is the headline summary of the SVD screen:
  - 60 `Estimate`
  - 12 `Fix (anchor)`
  - 25 `Fix (irrelevant)`
- The model is not simply identifiable or non-identifiable as a whole. Instead,
  different parameter groups have different levels of information support.

## Profile Likelihood: Representative Classes

Figure:

- `profile_likelihood_representative_3x3.png`

Tables:

- `profile_50d_balanced_relaxed_summary.csv`
- `profile_likelihood_class_counts.csv`

What it shows:

- Representative profile likelihood curves for the four profile classes after
  the SVD screen.
- Each curve shows how model fit changes when one parameter is fixed at
  different values while nuisance parameters are allowed to compensate.
- The panel does not show all 60 profiles. It shows selected examples for:
  practically identifiable, boundary-limited, weakly identifiable, and
  flat/non-identifiable classes.

Technical terms:

- **Profile likelihood**: a practical identifiability method where one
  parameter is fixed across a grid and the remaining nuisance parameters are
  re-optimized.
- **Nuisance parameter**: a parameter allowed to vary during profiling to test
  whether it can compensate for the fixed parameter.
- **Delta loss**: increase in loss relative to the best fit on the profile.
- **Practical identifiability**: whether a parameter can be estimated within a
  finite range given the data, noise assumptions, and model structure.
- **Boundary-limited**: the profile suggests an interval may extend beyond the
  tested parameter range.
- **Weakly identifiable**: the profile contains some information but not enough
  for a strong finite interval.
- **Flat/non-identifiable**: changing the parameter does not meaningfully worsen
  the fit.

Interpretation:

- SVD is used as a fast local screen.
- Profile likelihood is used as a nonlinear confirmation step.
- The current profile result summarizes the selected parameters as:
  - 51 practically identifiable
  - 6 boundary-limited
  - 1 weakly identifiable
  - 2 flat/non-identifiable
- The 3x3 panel should be read together with the profile summary table. The
  table gives the complete count; the figure provides interpretable examples
  of each class.

## How The Plots Fit Together

The figures should be read in this order:

```text
sensitivity_top_parameters.png
-> identifiability_singular_values.png
-> identifiability_nullspace_participation.png
-> identifiability_nullspace_participation_all_parameters.png
-> identifiability_sensitivity_vs_nullspace_all_parameters.png
-> identifiability_compensation_edges.png
-> identifiability_compensation_network.png
-> identifiability_compensation_network_sensitivity.png
-> identifiability_decision_map.png
-> identifiability_class_counts.png
-> profile_likelihood_representative_3x3.png
```

Together they show:

1. Which parameters affect outputs.
2. Which parameter directions are weakly informed.
3. Which parameters compensate for each other.
4. Which parameters should be estimated or fixed.
5. Which selected parameters are practically identifiable after nonlinear
   profile-likelihood confirmation.
6. Which weak or compensatory directions should motivate Bayesian experimental
   design.
