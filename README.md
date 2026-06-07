# MetRep: Metabolic–Reproductive Mechanistic Model  
### Mechanistic Systems Biology, Bayesian Inference, and Experimental Design for Complex Endocrine–Metabolic ODE Models

## Overview

MetRep is a mechanistic systems-biology framework for studying endocrine–metabolic regulation using a nonlinear ordinary differential equation (ODE) model of bovine physiology.

The model integrates metabolic and reproductive pathways, including glucose–insulin regulation, energy balance, ovarian dynamics, and hormonal feedback mechanisms. It was originally developed during my PhD work and has since been extended into a reproducible computational workflow for sensitivity analysis, identifiability diagnostics, uncertainty propagation, Bayesian inference, and Bayesian experimental design (BED).

This repository focuses on a key scientific challenge in mechanistic modelling:

> Which parameters are learnable, which biomarkers are most informative, and when should measurements be collected to reduce uncertainty in a complex nonlinear biological system?

The workflow combines mechanistic modelling, biological admissibility filtering, global sensitivity analysis, Bayesian posterior updating, and experimental design to study parameter learning in a large nonlinear endocrine–metabolic model.

---

## Scientific Motivation

Large mechanistic ODE models are difficult to calibrate because they often contain:

- many parameters,
- nonlinear feedback loops,
- threshold and Hill-type mechanisms,
- correlated effects,
- practical non-identifiability,
- expensive experimental measurements.

In such systems, it is not enough to ask:

> Which parameters influence the model?

We also need to ask:

> Which biomarkers are informative for parameter estimation?  
> At what time points should measurements be taken?  
> Can experimental design improve parameter learning?

MetRep addresses these questions through a biologically constrained Bayesian workflow.

---

## Model Scope

The model describes coupled endocrine–metabolic regulation and includes interactions among:

### Reproductive system
- Follicular growth
- Corpus luteum dynamics
- Ovarian regulation
- Estrous-cycle hormonal feedback

### Metabolic system
- Glucose regulation
- Insulin signaling
- IGF-1 interactions
- Energy partitioning
- Liver–blood nutrient exchange

### Hormonal pathways
- FSH
- LH
- Progesterone (P4)
- Estradiol (E2)
- Prostaglandin F2α (PGF)
- Inhibin (INH)
- IGF1
- Insulin
- Glucose
- Glucagon

---

## Final Workflow

The final workflow implemented in this repository is:

text Mechanistic ODE model     ↓ Local sensitivity analysis     ↓ Biological admissibility filtering     ↓ Global sensitivity analysis (admissible ensemble)     ↓ SVD identifiability     ↓ Profile likelihood     ↓ Uncertainty propagation     ↓ Bayesian inference & experimental design (BED)         ├── PhD-style posterior reweighting         ├── Reduced archive-based Bayesian posterior update         ├── Archive-based sequential ABC filtering         └── SMC+ML admissible-bank enrichment                 for:                 - stable MI/BED ranking                 - smaller MI error bars                 - improved admissible ensemble coverage 

---

## Biological Admissibility Filtering

A central feature of this work is the use of biological admissibility filtering.

Instead of treating all parameter samples as equally plausible, model simulations are filtered according to biological constraints and expected physiological behavior.

Only ODE-confirmed biologically plausible simulations are retained.

### Result

- Initial admissible models: 2,957
- Final enriched admissible ensemble: 12,721 ODE-confirmed simulations

This admissible ensemble is used to improve:

- global sensitivity robustness,
- uncertainty propagation,
- mutual-information (MI) stability,
- Bayesian experimental design ranking.

Importantly:

> The enriched 12,721-model ensemble is not used as the plotted Bayesian prior.

Instead:

- Broad ±5% prior → used for prior-to-posterior Bayesian updates
- 12,721 admissible ensemble → used for robust ranking, MI stability, and ensemble coverage

---

## Global Sensitivity on the Admissible Ensemble

Global sensitivity analysis was performed using the biologically admissible ensemble.

Methods:

- Partial Rank Correlation Coefficient (PRCC)
- Spearman correlation

across:

- 98 parameters
- 9 observable biomarkers

### Purpose

This step identifies:

> Which biomarkers are most associated with each parameter within biologically plausible model behavior

rather than over the unconstrained parameter space.

Example use:

text Parameter     ↓ Top associated biomarkers     ↓ Candidate measurements for Bayesian design 

Representative heatmaps are included in:

text results_final/figures/ 

---

## Identifiability Diagnostics

Parameter learnability was assessed using:

### SVD identifiability

to assess local structural identifiability trends.

### Profile likelihood

to evaluate practical identifiability and uncertainty structure.

The analysis revealed three representative parameter classes:

### Practically identifiable
Strong posterior learning expected.

Examples:
- insulin_glucose_threshold
- inhibin_clearance
- hp_p4_follicle_scale

### Boundary-limited
Partial learning or one-sided narrowing.

Examples:
- blood_to_liver_glucose_threshold
- gnrh_clearance
- hp_iof_threshold

### Weak / flat
Limited learning despite measurements.

Examples:
- insulin_igf_threshold
- feed_direct_blood_fraction
- lh_basal_release

These representative parameter classes are later used in Bayesian inference and BED.

---

## Uncertainty Propagation

Uncertainty propagation was used to identify:

> When the model is most informative

by analyzing ensemble variability over time.

Outputs include:

- predictive uncertainty bands,
- time-varying variability,
- candidate informative windows for measurement.

This step narrows the search space for Bayesian experimental design.

---

## Bayesian Inference & Experimental Design (BED)

The Bayesian workflow combines:

### 1. Posterior reweighting (PhD-style)

Broad ±5% prior simulations are reweighted according to observational likelihood.

Used to study:

text prior → posterior update 

under different measurement scenarios.

### 2. Guided Bayesian updates

Measurement scenarios are selected through:

text Profile likelihood     ↓ Global sensitivity     ↓ Uncertainty propagation     ↓ Mutual information (MI)     ↓ Bayesian posterior update 

For each representative parameter:

1. top biomarkers are selected from global sensitivity,
2. informative time windows are selected from uncertainty analysis,
3. MI/BED ranks candidate biomarker–day observations,
4. the broad prior is updated accordingly.

This produces:

### Parameter-specific guided posterior updates

rather than one generic observation strategy.

### 3. Archive-based sequential ABC filtering

Approximate Bayesian posterior evidence is generated using:

> archive-based sequential ABC filtering

Important caveat:

This is not full live ABC-SMC.

Instead:

- precomputed ODE simulations are reused,
- no ODE reruns are performed,
- no live particle perturbation is used.

The method provides an efficient likelihood-free approximation for expensive ODE systems.

---

## Main Scientific Finding

The key result of this repository is:

> Bayesian posterior learning confirms the identifiability diagnosis rather than rescuing weak parameters.

Observed behavior:

### Practically identifiable parameters
Strong posterior contraction.

### Boundary-limited parameters
Partial or one-sided learning.

### Weak/flat parameters
Remain broad despite guided observations.

This suggests:

> Bayesian experimental design improves learning where information exists, but does not magically solve practical non-identifiability.

---

## Repository Structure

text analyses/ │── local_sensitivity/ │── global_sensitivity/ │── uncertainty/ │── bayesian_inference/ │── bayesian_experimental_design/ │── surrogate_admissible_prior_enrichment/  docs/ │── methodology.md │── results_summary.md  results_final/ │── figures/ │── tables/ 

---

## Reproducibility

Example commands:

### Global sensitivity (98 × 9)

bash python analyses/global_sensitivity/run_global_sensitivity_enriched_98x9.py 

### Bayesian target selection

bash python analyses/bayesian_inference/select_bayesian_targets.py 

### Bayesian/BED workflow

bash python analyses/bayesian_experimental_design/surrogate_bed/run_targeted_bed_posterior_portfolio.py 

### Archive-based ABC filtering

bash python analyses/bayesian_inference/abc_smc/run_abc_smc.py 

---

## Methodological Caveats

- The 12,721-model admissible ensemble is not the Bayesian prior.
- Broad ±5% priors are used for prior-to-posterior visualizations.
- Archive-based ABC filtering is not full live ABC-SMC.
- Posterior learning is constrained by model identifiability.
- Profile likelihood results are reused from prior work and are not recalculated.

---

## Documentation

Additional details are available in:

- docs/methodology.md
- docs/results_summary.md
- results_final/README.md

---

## Citation / Context

This repository extends a mechanistic endocrine–metabolic systems model originally developed during my PhD work and expands it with a reproducible workflow for sensitivity analysis, identifiability, uncertainty propagation, Bayesian inference, and Bayesian experimental design.
