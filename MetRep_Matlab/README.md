# MATLAB Reference Implementation

This folder contains the original MATLAB reference implementation of the
BovSys/MetRep model, including the dexamethasone PK/PD extension.

The Python implementation in `../MetRep_Python/` is the open translated core
model for users without a MATLAB license. MATLAB remains the reference for the
Dexa extension until the Python Dexa module is implemented and validated.

## Files

- `BovSys_run_dexa_v3.m`: interactive MATLAB runner and plotting workflow
- `BovSys_Equa_dexa_v3.m`: 25-state ODE right-hand side
- `BovSys_para_dexa_v3.m`: 98-parameter vector used by the model
- `DM.mat`, `ML.mat`: DMI and milk forcing data
- `Data_P4_Hol.mat`, `Data_IGF_Hol.mat`: hormone comparison data

## Notes

The MATLAB runner uses interactive prompts and MATLAB plotting. It is included
for provenance and reference behavior, not as the primary reproducible Python
workflow.
