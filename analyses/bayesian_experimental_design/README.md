# Bayesian Experimental Design

The MATLAB BED material is retained outside the Python model workflow as
methodological provenance.

For public GitHub presentation, use the single clean BED script in
`matlab_original/BED_1M_ALL.m`. It calls
`matlab_original/BovSys_run_v3_baseline.m`, which uses the published v3 model
equations with Dexa PK/PD switched off. This keeps v3 as the canonical model
line while presenting BED as a baseline MetRep design analysis.

The current curated story should focus on:

- mutual information across candidate sampling days
- informative species ranking
- posterior narrowing under informative measurements
- timing or surrogate comparison, where relevant

Method summary:

- Sample possible parameter sets.
- Simulate the model for each parameter set.
- Compare candidate sampling days and measured species.
- Use mutual information to rank which measurements are expected to be most
  informative.
- Use posterior narrowing to show how informative measurements reduce
  uncertainty.

Within the full project logic, BED is the proposed way to improve information
about parameters or model directions that classical identifiability analysis
shows are weakly informed.

A compact Python BED reproduction can still be treated as future work.
