# Sensitivity And Identifiability

The analysis workflow uses three complementary stages.

These stages are the bridge between the model and Bayesian experimental design:
they reveal which outputs and parameters are well informed, and which parameter
directions need better experimental measurements.

## 1. Sensitivity Analysis

Sensitivity analysis asks which parameters have the largest local effect on
selected model outputs. It is useful as a first screen, but sensitivity alone is
not identifiability: a parameter can strongly affect outputs and still be hard
to estimate if another parameter can compensate for it.

Main script:

```bash
python MetRep_Python/scripts/04_run_sensitivity.py
```

## 2. SVD Identifiability Screen

The SVD workflow builds a local sensitivity matrix using selected measurable
outputs and analyzes its singular values and nullspace directions.

The key SVD outputs are:

- `identifiability_singular_values.png`: shows the singular-value spectrum. A
  wide drop toward small singular values indicates directions in parameter space
  that are weakly informed by the selected outputs.
- `identifiability_nullspace_participation.png`: ranks parameters by how much
  they participate in weak/nullspace directions. High participation means a
  parameter is involved in compensatory combinations.
- `identifiability_compensation_edges.png`: shows strong parameter-pair
  compensation relationships inferred from the nullspace.
- `identifiability_compensation_network.png`: shows the same compensation
  structure as a node-link graph, where nodes are parameters and edges indicate
  compensatory relationships.
- `identifiability_decision_map.png`: combines normalized sensitivity and
  nullspace involvement to support the `Estimate`, `Fix (anchor)`, and
  `Fix (irrelevant)` classification.

For a plot-by-plot explanation of the technical terms and conclusions, see
`docs/plot_interpretation_guide.md`.

It classifies parameters into:

- `Estimate`: sensitive and sufficiently separable
- `Fix (anchor)`: high-impact but compensatory; useful as an anchor
- `Fix (irrelevant)`: low-impact in the selected output setting

The class should guide the next modeling step:

| Class | Practical decision |
|---|---|
| `Estimate` | Include in calibration/profile likelihood, because the selected outputs carry enough local information. |
| `Fix (anchor)` | Do not freely estimate together with its compensation partners; fix it, constrain it with prior knowledge, or use it as an anchor. |
| `Fix (irrelevant)` | Keep fixed for this dataset and scenario; it may need a different experiment or output panel to become informative. |

The compensation network is especially useful for the `Fix (anchor)` group. A
node is a parameter, and an edge means two parameters appear together in a weak
direction. Strong edges indicate that the model can preserve similar outputs by
moving those parameters together. This is why compensation is not a failure of
the model; it is a design signal that the current measurements cannot separate
those mechanisms cleanly.

This is a fast local-linear diagnostic and is used to select parameters for
profile likelihood.

## 3. Profile Likelihood Confirmation

Profile likelihood fixes one selected parameter over a grid of values and
allows nuisance parameters to compensate. The resulting loss curve provides a
nonlinear practical-identifiability confirmation.

Current profile-likelihood summary:

- 60 profiled parameters
- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable: `insulin_igf_threshold`
- 2 flat/non-identifiable: `feed_direct_blood_fraction`, `lh_basal_release`

The recommended interpretation is that the selected measurable outputs identify
a substantial subset of parameters, while some metabolic/endocrine feedback
parameters require fixing, anchoring, or more informative experimental design.

## Connection To BED

The identifiability result motivates Bayesian experimental design. Parameters
that are weak, flat, boundary-limited, or compensatory are not simply failures;
they indicate where the current measurement set is not informative enough.

BED addresses the next question:

```text
Given the model and its weakly informed directions,
which sampling times and measured species would add the most information?
```

This makes the workflow constructive: identifiability diagnoses the problem,
and BED proposes how future experiments could improve it.

In practical terms, BED should be aimed at:

- parameters with strong nullspace participation;
- compensation pairs or clusters in the node graph;
- profile-likelihood cases that are weak, flat, or boundary-limited;
- outputs and sampling times likely to distinguish between compensatory
  mechanisms.
