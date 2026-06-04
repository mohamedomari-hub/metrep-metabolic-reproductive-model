# Sensitivity And Identifiability Analysis

This combined folder is retained as a legacy navigation entry. The curated
analysis structure is now split into:

- `analyses/local_sensitivity/`
- `analyses/identifiability/`

The executable scripts remain in `MetRep_Python/scripts/` to preserve imports
and reproducibility.

This module summarizes the staged Python analysis workflow:

```text
sensitivity analysis
-> SVD identifiability screen
-> profile likelihood confirmation
```

Curated public outputs are stored in `results_final/`.

Method summary:

- Sensitivity perturbs one parameter at a time and measures the output change.
- SVD identifiability decomposes the sensitivity matrix into informed and weak
  parameter directions.
- Compensation is inferred from parameters that appear together in weak
  directions.
- Profile likelihood fixes one parameter at different values and checks whether
  model fit gets worse after nuisance parameters are allowed to adjust.

The output of this module motivates BED: weak, flat, boundary-limited, or
compensatory directions identify where future experiments should be more
informative.
