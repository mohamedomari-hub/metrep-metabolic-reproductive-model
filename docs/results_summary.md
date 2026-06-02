# Results Summary

This repository curates the Python sensitivity and identifiability outputs. The
full MetRep scenario simulations are reported in the MetRep model paper, the
Dexa perturbation scenario is reported in the Dexa paper, and the full Bayesian
experimental design analysis is reported in the PhD thesis.

## Curated Repository Results

For detailed explanations of the plots and technical terms, see
`docs/plot_interpretation_guide.md`.

| Step | Question | Main curated output |
|---|---|---|
| Python validation | Does Python metadata match MATLAB reference? | `results_final/tables/validation_report.csv` |
| Baseline run | Does the translated Python model run? | `results_final/figures/baseline_selected_states.png` |
| Sensitivity | Which parameters have the largest local output effect? | `results_final/tables/sensitivity_top_parameters.csv`, `results_final/figures/sensitivity_top_parameters.png` |
| SVD identifiability | Which parameters are estimable, anchors, or irrelevant? | `results_final/tables/identifiability_class_counts.csv`, `results_final/figures/identifiability_decision_map.png` |
| SVD nullspace | Which parameters compensate or sit in weak directions? | `results_final/figures/identifiability_singular_values.png`, `results_final/figures/identifiability_nullspace_participation.png`, `results_final/figures/identifiability_compensation_edges.png`, `results_final/figures/identifiability_compensation_network.png` |
| Profile likelihood | Which selected parameters are practically identifiable? | `results_final/tables/profile_50d_balanced_relaxed_summary.csv`, `results_final/figures/profile_50d_balanced_relaxed_combined_profiles.png` |

## Identifiability Headline

Using measurable outputs
`FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon`, the SVD screen
classified the analyzed parameters as:

- 60 `Estimate`
- 12 `Fix (anchor)`
- 25 `Fix (irrelevant)`

Profile likelihood was then used as a nonlinear confirmation step for the
selected `Estimate` parameters:

- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable
- 2 flat/non-identifiable

## How To Use The Identifiability Classes

The SVD result should be read as an estimation strategy, not only as a label
for each parameter.

| Class | Meaning | Recommended use |
|---|---|---|
| `Estimate` | The parameter is sensitive enough and sufficiently separated from weak/nullspace directions. | Candidate for calibration and profile likelihood. |
| `Fix (anchor)` | The parameter affects outputs, but it also participates strongly in compensation directions. | Keep fixed, constrain with prior knowledge, or use as an anchor while estimating other parameters. |
| `Fix (irrelevant)` | The parameter has low influence for the selected measurable outputs and scenario. | Do not estimate from this dataset; keep at reference value unless a different experiment targets it. |

This does not mean that `Fix (irrelevant)` parameters are biologically
unimportant. It means the selected outputs and simulation window do not provide
enough useful information about them in this analysis.

## Compensation And Node Graph

The compensation tables and network plot explain why some parameters cannot be
estimated together without additional information.

- `identifiability_compensation_edges.png` lists the strongest parameter pairs
  that can compensate for each other.
- `identifiability_compensation_network.png` shows those relationships as a
  graph: nodes are parameters, edges are compensation relationships, node color
  is the SVD class, and thicker edges indicate stronger compensation.

Examples of strong compensation relationships include
`fsh_syn_scale` with `insulin_fsh_scale`, `blood_usage_max` with
`insulin_clearance`, and lactation/storage terms such as
`lactation_oxt_decay` with `fat_mobilization_storage_threshold`. These pairs
are important because changing one parameter can be partly offset by changing
the other, producing similar measured outputs.

Operationally, compensation pairs should be handled in one of three ways:

1. Fix or tightly constrain one parameter in a compensatory pair.
2. Estimate only the better-supported member of the pair.
3. Use Bayesian experimental design to choose new sampling times or species
   that break the compensation.

## Link To Bayesian Experimental Design

The sensitivity and identifiability analyses define the problem: some
parameters are estimable, while others are weak, compensatory, or irrelevant
under the current measurement set. BED is the constructive next step. It asks
which candidate measurements would add the most information about the weakly
informed directions.

In this project story:

```text
Sensitivity finds influential parameters.
SVD identifies weak/nullspace and compensation structure.
Profile likelihood confirms practical identifiability for selected estimates.
BED proposes future measurements to improve the weak or compensatory directions.
```

## Publication-Reported Results

- MetRep scenario simulations: reported in the MetRep model paper.
- Dexa perturbation simulation/validation: reported in the Dexa paper.
- Bayesian experimental design: reported in the PhD thesis; selected figures
  are included in `results_final/figures/` for context.
