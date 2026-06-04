# Model Overview

## Introduction

The MetRep model is a mechanistic systems-biology model describing the coupled metabolic and reproductive physiology of dairy cattle using a system of nonlinear ordinary differential equations (ODEs).

The model was originally developed during a PhD project in applied mathematics and systems biology to investigate the interaction between:

- energy metabolism
- glucose–insulin regulation
- reproductive endocrine dynamics
- nutritional status
- pharmacological perturbation (dexamethasone)

The framework integrates physiological processes across multiple biological scales and allows simulation of both baseline endocrine–metabolic regulation and stress/perturbation scenarios such as negative energy balance or glucocorticoid administration.

The repository contains:

1. the original MATLAB implementation (scientific reference),
2. a translated Python implementation for open and reproducible modeling,
3. a complete model-analysis workflow, including:
   - local sensitivity analysis,
   - structural/practical identifiability diagnostics,
   - profile likelihood analysis,
   - uncertainty propagation,
   - global parameter screening,
   - Bayesian experimental design (BED).

---

# Biological Scope

The model describes how metabolic state influences reproductive dynamics in dairy cattle.

Biologically, the system couples two major subsystems:

## 1. Reproductive Endocrine Axis

The reproductive component describes hormonal interactions controlling follicular development and luteal dynamics during the estrous cycle.

It includes mechanistic feedback loops involving:

- hypothalamic signaling,
- pituitary hormone secretion,
- ovarian follicle growth,
- corpus luteum dynamics,
- steroid hormone regulation.

Key modeled hormones include:

- GnRH (gonadotropin releasing hormone)
- FSH (follicle stimulating hormone)
- LH (luteinizing hormone)
- P4 (progesterone)
- E2 (estradiol)
- Inhibin
- PGF (prostaglandin F2α)
- IGF-1 (insulin-like growth factor 1)

These components regulate:

- follicle maturation,
- ovulation timing,
- luteolysis,
- hormonal feedback control.

The reproductive subsystem contains nonlinear threshold and feedback mechanisms to reproduce oscillatory estrous-cycle behavior.

---

## 2. Metabolic Regulation

The metabolic component represents glucose homeostasis and energy allocation.

It includes interactions among:

- blood glucose
- insulin
- glucagon
- hepatic glucose production
- glucose storage
- body fat reserves

The metabolic subsystem mechanistically represents:

- glucose uptake,
- glycogen dynamics,
- endocrine control of metabolism,
- energetic stress,
- nutritional inputs.

A key modeling motivation is understanding how altered metabolic state influences reproductive performance.

---

# State Variables

The baseline MetRep model contains 22 dynamic states.

## Reproductive / Endocrine States

| State | Description |
|--------|-------------|
| GnRH | Gonadotropin releasing hormone |
| FSH | Follicle stimulating hormone |
| LH | Luteinizing hormone |
| Follicle | Dominant follicle development |
| PGF | Prostaglandin F2α |
| CL | Corpus luteum |
| P4 | Progesterone |
| E2 | Estradiol |
| INH | Inhibin |
| Enzyme proxy | Enzymatic regulation term |
| OXT | Oxytocin |
| IOF | Intra-ovarian factor |
| IGF1 | Insulin-like growth factor 1 |

## Metabolic States

| State | Description |
|--------|-------------|
| Insulin | Plasma insulin |
| Glucose | Blood glucose |
| Fat | Fat reserve compartment |
| Liver glucose | Hepatic glucose pool |
| Glucose storage | Glycogen/glucose storage |
| Glucagon | Plasma glucagon |

Additional states represent intermediary physiological dynamics and coupling processes between endocrine and metabolic regulation.

---

# Dexamethasone Extension

The repository also contains an optional dexamethasone (Dexa) PK/PD extension, increasing the model from 22 to 25 states.

This extension was developed to investigate how glucocorticoid administration perturbs metabolic and reproductive regulation.

The Dexa component introduces a mechanistic pharmacokinetic–pharmacodynamic (PK/PD) model consisting of three additional states:

| State | Description |
|--------|-------------|
| A_dep | Intramuscular depot amount |
| A_cent | Central/systemic dexamethasone amount |
| C_e | Effect-site concentration |

The PK model includes:

- absorption from intramuscular depot,
- systemic circulation,
- delayed effect compartment dynamics.

The PD component mechanistically perturbs:

- glucose metabolism,
- glucagon signaling,
- insulin regulation,
- endocrine responses.

This extension allows simulation of:

- acute glucocorticoid administration,
- delayed metabolic responses,
- endocrine disruption,
- recovery trajectories.

The Dexa extension is optional and separated from the baseline physiological model.

---

# Parameters

The baseline model contains 98 mechanistic parameters describing:

- secretion rates,
- degradation rates,
- Hill-type regulatory effects,
- feedback strengths,
- metabolic conversion rates,
- threshold functions,
- hormonal coupling dynamics.

These parameters were originally calibrated using experimental observations and physiological assumptions.

Because mechanistic models of this size can contain partially informed parameter combinations, the repository includes a full parameter-diagnostics workflow.

---

# Simulation Scenarios

The model can simulate multiple physiological scenarios.

## Baseline Non-Lactating Scenario

This is the main reference configuration used for:

- local sensitivity analysis,
- identifiability diagnostics,
- profile likelihood analysis,
- global parameter screening,
- uncertainty propagation,
- Bayesian experimental design.

It represents a physiologically stable non-lactating animal under standard nutritional conditions.

---

## Lactating Scenario

Lactation introduces additional energetic demands and altered metabolic regulation.

These scenarios are used to investigate how energetic burden affects reproductive dynamics.

---

## Acute Negative Energy Balance

Short-term energetic stress is simulated to study transient disruption of metabolism and endocrine regulation.

---

## Chronic Negative Energy Balance

Long-term energetic deficits are simulated to assess persistent physiological consequences.

---

## Dexamethasone Perturbation

Pharmacological perturbation scenarios simulate the effect of exogenous glucocorticoid administration on system dynamics.

These simulations are used as an external validation and perturbation framework.

---

# Model Analysis Workflow

The repository implements a complete systems-modeling workflow built around the mechanistic model.

text Mechanistic ODE Model         ↓ Python Translation         ↓ Local Sensitivity Analysis         ↓ Identifiability Diagnostics         ↓ Global Sensitivity Screening         ↓ Uncertainty Propagation         ↓ Bayesian Experimental Design         ↓ Dexa Perturbation / Validation 

---

## Local Sensitivity Analysis

Local sensitivity quantifies how small perturbations in parameters affect model outputs around a nominal physiological state.

Characteristics:

- 98 parameters
- one-at-a-time perturbation
- +1% parameter change
- AUC-based trajectory response
- baseline non-lactating scenario

This analysis identifies parameters with strong local influence on measurable outputs.

---

## Identifiability Analysis

Identifiability diagnostics determine whether model parameters can theoretically or practically be estimated from available measurements.

Implemented approaches include:

### SVD-Based Structural Screening

Uses trajectory sensitivities to evaluate:

- rank deficiency,
- parameter compensation,
- singular-value structure,
- poorly informed parameter directions.

### Profile Likelihood

Evaluates practical identifiability by:

- fixing one parameter,
- re-optimizing nuisance parameters,
- evaluating objective deterioration.

This identifies:

- identifiable parameters,
- weakly informed parameters,
- parameter compensation effects.

---

## Global Sensitivity Screening

Global parameter screening investigates parameter–output relationships across an ensemble of biologically plausible model realizations.

Unlike local sensitivity, this evaluates variability across the parameter space.

Current implementation uses:

- stored Monte Carlo simulation bank,
- biological admissibility filtering,
- PRCC/Spearman screening.

The current workflow is variance-based screening, not strict Sobol decomposition, because the simulation bank uses ordinary Monte Carlo sampling rather than Saltelli/Sobol sampling.

---

## Uncertainty Propagation

Uncertainty propagation evaluates how admissible parameter uncertainty translates into prediction uncertainty.

Outputs include:

- median trajectories,
- prediction bands,
- uncertainty envelopes,
- nominal-reference trajectories.

This quantifies model robustness under biologically realistic variability.

---

## Bayesian Experimental Design (BED)

BED evaluates:

> Which future measurements would maximally reduce parameter uncertainty?

The PhD workflow computes expected information gain from candidate measurements.

It identifies:

- informative biomarkers,
- optimal measurement timing,
- weakly informed model regions.

BED therefore complements identifiability analysis by proposing experiments that improve parameter estimation.

---

# MATLAB and Python Roles

## MATLAB

MATLAB contains the original scientific implementation used during model development and thesis work.

It serves as the reference implementation for:

- published simulations,
- Dexa studies,
- original BED workflow.

---

## Python

Python provides an open and reproducible implementation designed for:

- transparency,
- reproducibility,
- accessibility without MATLAB,
- modern scientific workflows.

It additionally supports:

- structured analyses,
- reproducible outputs,
- GitHub documentation,
- extensible scientific workflows.

---

# Model Flowchart

The figure below summarizes the physiological architecture of the MetRep model.

<img width="3872" height="1990" alt="MetRep flowchart" src="https://github.com/user-attachments/assets/9e2b9585-b882-446a-8113-839e2704f695" />

---

# Scientific Context

This repository aims to demonstrate not only a mechanistic physiological model, but also a complete systems-modeling workflow, combining:

- mechanistic biology,
- nonlinear dynamical systems,
- parameter diagnostics,
- uncertainty quantification,
- experimental design,
- pharmacological perturbation modeling.

The goal is to support interpretable and physiologically grounded modeling for dairy-cow endocrine–metabolic systems.
