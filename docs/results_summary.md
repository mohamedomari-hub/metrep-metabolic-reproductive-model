# Results Summary

## Main Modelling Diagnostics

- Local sensitivity identifies nominal mechanistic influence.
- Profile likelihood evaluates practical estimability.
- PRCC and Spearman summarize global parameter-biomarker associations across
  biologically admissible simulations.
- Uncertainty propagation quantifies local robustness around the calibrated
  model.
- Bayesian experimental design identifies informative measurement choices.

## Key Interpretation

Parameters such as `insulin_glucose_threshold` and
`blood_to_liver_glucose_threshold` show strong sensitivity and global
association patterns, while weak or flat profile-likelihood parameters
generally show limited diagnostic signal. This supports a coherent distinction
between influential, boundary-limited, and weakly identifiable mechanisms.

The uncertainty bands are based on a narrow, biologically filtered `+/-0.5%`
ensemble and therefore describe local robustness rather than population-level
variability. BED results are maintained separately under
`analyses/bayesian_experimental_design/` while the expanded prior-based run is
being regenerated.
