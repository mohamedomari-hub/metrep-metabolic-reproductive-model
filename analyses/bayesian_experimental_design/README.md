# Bayesian Experimental Design

The MATLAB BED material is retained outside the Python model workflow as
methodological provenance.

For public GitHub presentation, use the single clean BED script in
`matlab_original/BED_1M_ALL.m`. It calls
`matlab_original/BovSys_run_v3_baseline.m`, which uses the published v3 model
equations with Dexa PK/PD switched off. This keeps v3 as the canonical model
line while presenting BED as a baseline MetRep design analysis.

The current curated story should focus on:

- BED as the experimental-design response to weak or compensatory
  identifiability findings
- the original MATLAB implementation as PhD provenance
- the future Python surrogate workflow as a reproducibility extension

Method summary:

- Sample possible parameter sets.
- Simulate the model for each parameter set.
- Compare candidate sampling days and measured species.
- Use mutual information to rank which measurements are expected to be most
  informative.
- Use posterior narrowing to show how informative measurements reduce
  uncertainty.

The current RF/surrogate comparison figures are not included in
`results_final/`. They should be reported only after the surrogate validation,
mutual-information convergence, and biological admissibility checks are
documented.

Posterior summary:

- Draw prior samples for the selected uncertain parameters.
- Simulate the model for each prior sample at the candidate sampling day.
- Generate one fixed synthetic observation from the nominal/reference
  simulation plus Gaussian measurement noise.
- Compute a Gaussian likelihood for each prior sample.
- Normalize the likelihoods into weights.
- Estimate the posterior for the target parameter as a weighted KDE of the
  prior samples.

Mutual-information summary:

- Build Monte Carlo pairs of target quantity and candidate measurement.
- Estimate marginal and joint densities for the target, measurement, and their
  joint distribution.
- Rank candidate designs by the average log-density ratio
  `log(p(target, measurement) / (p(target) p(measurement)))`.

Within the full project logic, BED is the proposed way to improve information
about parameters or model directions that classical identifiability analysis
shows are weakly informed.

A compact Python BED reproduction can still be treated as future work. The
recommended scaffold is in `surrogate_bed/`.
