# Model Overview

The BovSys/MetRep model is a mechanistic ODE model for dairy cow
metabolic-reproductive dynamics. It couples estrous-cycle endocrine regulation
with glucose-insulin-glucagon metabolism and feeding or lactation inputs.

The Python implementation in this repository currently translates the 22-state
non-Dexa core model. The original MATLAB implementation also includes a
dexamethasone PK/PD extension with three additional states:

- `A_dep`: intramuscular depot amount
- `A_cent`: central/systemic dexamethasone amount
- `C_e`: effect-site concentration

## State Groups

- Reproductive/endocrine states: GnRH, FSH, LH, follicle, PGF, CL, P4, E2,
  inhibin, enzyme proxy, OXT, IOF, IGF-1
- Metabolic states: insulin, blood glucose, fat, liver glucose, glucose
  storage, glucagon
- Dexa extension states in MATLAB reference: depot, central, and effect-site
  compartments

## Project Logic

The model is used as the center of a complete systems-modeling workflow:

```text
mechanistic model
-> Python translation
-> sensitivity and identifiability analysis
-> Bayesian experimental design to improve weakly informed parameters
-> Dexa perturbation as validation/extension
```

The classical analyses identify which parts of the parameter space are well
informed by the current outputs. Bayesian experimental design then asks what
additional measurements would make the model more informative. The Dexa
extension provides an external pharmacological perturbation scenario for model
validation.

## Implementation Roles

- MATLAB is the original reference implementation.
- Python is the open translated implementation for users without MATLAB.
- Sensitivity, identifiability, and profile likelihood are implemented on the
  Python core model.
- MetRep scenario simulations are reported in the MetRep model paper.
- BED is reported in the PhD thesis. The MATLAB BED folder contains one clean
  v3 baseline BED port, while surrogate-assisted BED figures are excluded from
  `results_final` until the validation workflow is documented.
- Dexa perturbation behavior is reported in the Dexa paper. The MATLAB
  reference code is retained here; Python Dexa implementation is planned.
