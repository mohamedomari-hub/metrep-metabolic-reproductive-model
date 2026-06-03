# Dexa Extension Plan

The MATLAB reference implementation includes a dexamethasone PK/PD module. The
Python implementation currently covers the 22-state non-Dexa core model.

In the project logic, Dexa is the perturbation-validation endpoint. The Dexa
simulation result is reported in the Dexa paper. This repository retains the
MATLAB reference implementation and documents the planned Python extension.

## MATLAB Reference Dexa Module

The Dexa module adds three states:

- `A_dep`: intramuscular depot amount
- `A_cent`: central/systemic amount
- `C_e`: effect-site concentration

The module uses first-order absorption, first-order elimination, and an
effect-site compartment. Pharmacodynamic effects alter glucagon secretion and
selected glucose-utilization terms.

## Python Implementation Plan

1. Extend the Python state vector from 22 to 25 states.
2. Add Dexa PK parameters and dose configuration.
3. Implement dose-event handling for day 0 and delayed dosing.
4. Add PD multipliers for glucagon secretion and glucose utilization.
5. Add non-lactating and lactating Dexa scenario runners.
6. Validate Python Dexa trajectories against MATLAB reference outputs.
7. Add Dexa-specific figures and curated result tables.

Until these steps are complete, the Python Dexa script should be described as
planned/in progress rather than completed.

## Validation Goal

The intended Dexa validation question is:

```text
When the model is challenged with dexamethasone administration,
does it reproduce the expected direction and timing of metabolic responses?
```

This makes Dexa different from an ordinary scenario simulation. It is a
stress-test of the model structure after the baseline model has been analyzed
and the information limitations have been characterized.
