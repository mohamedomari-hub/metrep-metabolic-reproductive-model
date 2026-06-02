# MATLAB BED Reference

This folder contains one clean MATLAB Bayesian experimental design (BED)
workflow for the GitHub repository.

Run:

```matlab
BED_1M_ALL
```

This script uses the published `v3` model equations, with Dexa PK/PD switched
off through `BovSys_run_v3_baseline.m`. In other words, BED is presented as a
baseline MetRep experimental-design analysis, not as a Dexa perturbation
simulation.

## Files

- `BED_1M_ALL.m`: clean GitHub-facing BED script using the published v3
  baseline model interface
- `BovSys_run_v3_baseline.m`: non-interactive v3 baseline runner used by the
  BED script; returns `[T,Y,OvT]` and keeps Dexa disabled
- `BovSys_para_dexa_v3.m`: v3 parameter vector used by the BED script
- `BovSys_Equa_dexa_v3.m`: v3 ODE right-hand side used by the baseline runner

The repository story uses v3 as the canonical MATLAB model line and exposes the
non-interactive baseline interface required by BED through
`BovSys_run_v3_baseline.m`.

The public project story uses BED as the experimental-design response to the
identifiability analysis:

```text
Sensitivity / identifiability analysis
-> identifies weak or non-identifiable parameter directions
-> Bayesian experimental design asks which sampling days and measured species
   would provide more information about those directions
```

A compact Python BED reproduction can be added later using the translated
Python MetRep model.
