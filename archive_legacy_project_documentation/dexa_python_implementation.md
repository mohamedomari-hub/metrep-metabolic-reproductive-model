# Dexa Python Implementation

The Python model now includes an optional dexamethasone perturbation workflow
translated from the MATLAB v3 files:

- `MetRep_Matlab/BovSys_Equa_dexa_v3.m`
- `MetRep_Matlab/BovSys_para_dexa_v3.m`
- `MetRep_Matlab/BovSys_run_dexa_v3.m`

The standard Python simulations still run the 22-state non-Dexa core model.
Dexa is activated only when `MetRep_Python/model_running/07_run_dexa_scenarios.py`
is used.

## State Extension

The 22 MetRep states are preserved. Dexa adds three states:

- `A_dep`: intramuscular depot amount of Dexa, ng
- `A_cent`: central/systemic Dexa amount, ng
- `C_e`: effect-site Dexa concentration, ng/mL

The extended model therefore has 25 states. The Dexa states are appended after
the core state vector, so the original state order remains unchanged for
baseline simulations, sensitivity analysis, SVD identifiability, and profile
likelihood.

## PK/PD Equations

For a dose `D` placed into the intramuscular depot, the MATLAB reference uses:

```text
dA_dep/dt = -ka A_dep
dA_cent/dt = F ka A_dep - ke A_cent
C = A_cent / Vd
dC_e/dt = keo (C - C_e)
```

with:

```text
ka = 13.4352
ke = 2.7086
F = 0.72
keo = 0.7
Vd = 1.105 L/kg * 600 kg
D = 0.02 mg/kg * 600 kg
```

The effect-site concentration modifies the metabolic model using the same
pharmacodynamic multipliers as the MATLAB code:

```text
Effect_gluca = 1 + 3 * C_e^10 / (C_e^10 + 1.8^10)
Effect_bt    = 1 -     C_e^7  / (C_e^7  + 1.8^7)
```

`Effect_gluca` stimulates glucagon secretion. `Effect_bt` inhibits selected
glucose storage, liver-to-fat, and milk-related blood glucose use terms.

## Running Dexa In Python

Run the short non-lactating Dexa perturbation matching the MATLAB day-0 dose
logic:

```bash
python MetRep_Python/model_running/07_run_dexa_scenarios.py \
  --scenario baseline_non_lactating \
  --days 3 \
  --dose-day 0 \
  --figure-dir results_final/figures \
  --table-dir results_final/tables \
  --prefix dexa_non_lactating_standard_3d
```

Run a lactating perturbation with dose at day 50:

```bash
python MetRep_Python/model_running/07_run_dexa_scenarios.py \
  --scenario lactating_c0_20 \
  --dose-day 50 \
  --figure-dir results_final/figures \
  --table-dir results_final/tables \
  --prefix dexa_lactating_c0_20_day50
```

The script saves:

- a Dexa versus no-Dexa response figure
- a response summary table with min/max values and maximum Dexa-baseline
  differences for plotted states

Add `--save-trajectories` if full simulation CSV/NPZ outputs are needed.

## Separation From Identifiability Analysis

The sensitivity, SVD identifiability, and profile-likelihood scripts analyze
only the core parameter list in `model_definition.parameters.PARAMETERS`.
That list contains the 98 metabolic-reproductive parameters from the MATLAB
model. Dexa PK/PD constants live in `model_definition.dexa.DexaConfig` and are
not included in the baseline identifiability workflow.

This separation is deliberate:

```text
core MetRep analysis = 22 states, 98 core parameters, no Dexa PK/PD constants
Dexa perturbation    = optional 25-state simulation using fixed PK/PD constants
```

This allows users to reproduce the sensitivity/SVD/profile-likelihood results
without introducing Dexa-specific parameters, while still allowing the model to
be challenged with Dexa as an external validation scenario.
