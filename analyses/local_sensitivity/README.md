# Local Sensitivity

Local sensitivity measures the nominal one-at-a-time AUC response of observable
model outputs to a `+1%` parameter perturbation.

Executable workflow:

```bash
python MetRep_Python/scripts/04_run_sensitivity.py
```

Curated outputs:

- `results_final/tables/metrep_non_lactating_standard_local_sensitivity_allparams.csv`
- `results_final/figures/sensitivity_top_parameters.png`

See `docs/methodology.md` for mathematical definitions and interpretation.
