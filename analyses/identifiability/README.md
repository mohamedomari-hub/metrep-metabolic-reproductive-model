# Identifiability

Identifiability analysis combines local sensitivity derivatives, SVD/nullspace
screening, compensation diagnostics, and nonlinear profile likelihood.

Executable workflows:

```bash
python MetRep_Python/scripts/05_run_identifiability.py
python MetRep_Python/scripts/10_run_profile_likelihood.py
```

Curated outputs are stored in `results_final/figures/` and
`results_final/tables/`. See `docs/methodology.md` for mathematical definitions
and interpretation.
