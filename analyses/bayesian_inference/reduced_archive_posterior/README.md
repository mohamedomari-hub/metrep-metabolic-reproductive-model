# Reduced Archive Posterior

This workflow performs selected-parameter posterior inference by weighting and
resampling rows from the saved broad `+/-5%` ODE archive.

It does not train a surrogate and does not run new ODE simulations. Posterior
support is therefore limited to precomputed ODE-confirmed archive rows.

