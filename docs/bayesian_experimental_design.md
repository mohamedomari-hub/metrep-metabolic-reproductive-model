# Bayesian Experimental Design

The Bayesian experimental design analysis evaluates which sampling times and
measured species provide the most information about model quantities or
parameters.

In this project, BED is used as the experimental-design answer to the
identifiability analysis. Sensitivity and identifiability diagnose which model
directions are weakly informed by existing outputs; BED asks how future
measurements should be chosen to improve those directions.

## Current Repository Status

The original PhD BED implementation is MATLAB code. It is retained for
provenance and methodological traceability.

For a clearer GitHub presentation, the repository also includes a v3 baseline
port:

- `analyses/bayesian_experimental_design/matlab_original/BED_1M_ALL.m`
- `analyses/bayesian_experimental_design/matlab_original/BovSys_run_v3_baseline.m`

This port uses the published v3 model equations with Dexa PK/PD switched off.
It should be interpreted as BED on the baseline MetRep model, not as a Dexa
simulation.

The BED workflow is still not the recommended lightweight reproducibility path,
because:

- the original script is large and combines several analysis variants
- it uses MATLAB parallel loops
- it depends on Statistics and Machine Learning Toolbox functions
- it contains hardcoded historical cluster paths
- it saves intermediate `.mat` files in the working directory
- exact repeatability is limited because the original script does not define a
  clean public random-seed/output convention

## Recommended Public Presentation

Present v3 as the canonical model line. Use the single clean MATLAB BED script
`analyses/bayesian_experimental_design/matlab_original/BED_1M_ALL.m`, which is
the v3 baseline port.

The full BED result is reported in the PhD thesis. This repository includes
selected BED figures for context:

- mutual information by candidate sampling day
- per-species information ranking
- posterior narrowing for informative designs
- timing/surrogate comparison, if included

A compact Python reproduction can be added later. It should reproduce the main
BED story rather than port the original MATLAB file line by line.

## Interpretation In The Project

The intended message is:

```text
Classical identifiability analysis shows the current measurement set does not
equally inform all model parameters.

Bayesian experimental design ranks candidate sampling times and measured
species by expected information gain, suggesting how future experiments can be
made more informative.
```

This connects the identifiability and BED parts into one workflow rather than
two unrelated analyses.

## Link Back To Identifiability

The BED section should be presented as the answer to the specific weaknesses
found by sensitivity and identifiability:

| Identifiability finding | BED interpretation |
|---|---|
| Sensitive and separable parameters | Current output panel is informative enough; these can be estimated and checked by profile likelihood. |
| High-impact compensation pairs | Future designs should target sampling times/species that separate the paired mechanisms. |
| Strong nullspace participation | Add measurements expected to reduce uncertainty in those weak directions. |
| `Fix (irrelevant)` parameters | Do not spend estimation effort on them unless BED suggests a different output/time window can make them informative. |
| Weak, flat, or boundary-limited profiles | Use BED to propose more informative observations before claiming precise estimates. |

This framing keeps BED from looking like an unrelated extra analysis. The
classical analysis says where the model is under-informed; BED says how a
future experiment could improve that information.
