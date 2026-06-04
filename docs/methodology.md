# Model Diagnostics, Identifiability, Uncertainty, and Bayesian Experimental Design

This document describes the analysis workflow used to evaluate the BovSys/MetRep metabolic-reproductive model.

The workflow follows a modelling narrative:

text Mechanistic ODE model → local sensitivity → SVD identifiability screen → profile likelihood confirmation → global sensitivity / admissible-bank association → uncertainty propagation → Bayesian experimental design 

The goal is not to use one single diagnostic, but to combine complementary methods:

- Sensitivity analysis asks which parameters influence outputs.
- Identifiability analysis asks whether influential parameters can be estimated uniquely.
- Global sensitivity screening asks which parameters drive variability across simulation ensembles.
- Uncertainty propagation asks how parameter uncertainty affects model predictions.
- Bayesian experimental design asks which measurements would reduce uncertainty most.

Different perturbation scales are used because each method answers a different question.

---

# 1. Local Sensitivity Analysis

Local sensitivity analysis asks:

> Which parameters have the largest local effect on selected model outputs near the nominal parameter set?

It is a useful first screen, but sensitivity alone is not identifiability. A parameter can strongly affect outputs and still be difficult to estimate if another parameter can compensate for it.

## Calculation

The workflow is:

1. Run the model with the reference parameter vector.
2. Perturb one parameter while keeping all others fixed.
3. Re-run the model.
4. Compare the output change.
5. Repeat for each parameter.

The public sensitivity table uses one-at-a-time relative sensitivity of output AUC values.

For output (y_k(t;\theta)), the area under the curve is

[
\mathrm{AUC}k(\theta)
=
\int{t_0}^{t_f} y_k(t;\theta),dt .
]

For parameter (\theta_j), the forward relative sensitivity is

[
S_{kj}
=
\frac{
\left[
\mathrm{AUC}k(\theta^{(j,+)})
-
\mathrm{AUC}k(\theta)
\right]
/
\mathrm{AUC}k(\theta)
}{h},
]

where

[
\theta^{(j,+)}j = (1+h)\theta_j,
\qquad
\theta^{(j,+)}\ell = \theta\ell
\quad \text{for } \ell \neq j.
]

In the current script,

[
h = 0.01,
]

corresponding to a (+1%) one-at-a-time perturbation.

This is a local diagnostic because it evaluates model response around the reference parameter set, not across the full biologically plausible parameter space.

Main script:

bash python MetRep_Python/scripts/04_run_sensitivity.py 

---

# 2. SVD Identifiability Screen

The SVD identifiability screen asks:

> Which parameter directions are locally informed by the selected measurable outputs, and which directions are weak or compensatory?

This is a fast local-linear diagnostic used before profile likelihood analysis.

## Sensitivity Matrix

The workflow builds a stacked trajectory sensitivity matrix. Rows correspond to output features across time, and columns correspond to parameters.

For parameter (\theta_j), selected output vector (Y(\theta)), and step size (h_j), the central-difference column is

[
S{:,j}
=
\frac{
Y(\theta^{(j,+)})
-
Y(\theta^{(j,-)})
}{2h_j},
]

where

[
\theta^{(j,+)}j = \theta_j + h_j,
\qquad
\theta^{(j,-)}j = \theta_j - h_j.
]

The finite-difference step is

[
h_j =
\max
\left(
h{\min},
h{\mathrm{rel}}|\theta_j|
\right).
]

In the default identifiability script,

[
h{\mathrm{rel}} = 10^{-3}.
]

If a nominal parameter value is zero, the implementation uses a nonzero fallback before applying the absolute minimum step.

The vector (Y(\theta)) is formed by stacking all selected output trajectories over sampled time points.

## Singular Value Decomposition

The sensitivity matrix is decomposed as

[
S = U\Sigma V^\top,
]

where:

- (S) is the stacked sensitivity matrix,
- (U) contains output-space directions,
- (\Sigma) contains singular values,
- (V^\top) contains parameter-space directions.

The numerical rank is defined using

[
\tau = c\sigma_1,
]

[
r =
\left|
\left{
\sigma_i : \sigma_i > \tau
\right}
\right|,
]

where (\sigma_1) is the largest singular value and (c) is the relative tolerance. The default tolerance is

[
c = 10^{-8}.
]

Large singular values indicate parameter combinations that strongly affect the selected outputs. Small singular values indicate weakly informed or nearly compensatory parameter directions.

## Nullspace Participation

The approximate nullspace is taken from the right singular vectors after the numerical rank:

[
\mathcal{N}
=
\left{
v_i^\top : i > r
\right}.
]

The nullspace participation score for parameter (j) is

[
n_j
=
\left(
\sum_q \mathcal{N}{qj}^2
\right)^{1/2}.
]

A high (n_j) means that parameter (j) participates strongly in weak or compensatory directions.

## Local SVD Ranking Score

The curated result uses a relative column-norm score:

[
s_j
=
\left|
\frac{\theta_j}{\max(|Y|,\epsilon)}
S{:,j}
\right|2 .
]

This scales the sensitivity column by parameter magnitude and output scale, making parameters more comparable.

## Parameter Classification

The SVD screen combines normalized sensitivity and normalized nullspace participation.

[
\bar{s}j
=
\frac{s_j-\min(s)}{\max(s)-\min(s)},
]

[
\bar{n}j
=
\frac{n_j-\min(n)}{\max(n)-\min(n)}.
]

Quantile thresholds are then defined as

[
s{\mathrm{hi}} = Q{0.75}(\bar{s}),
\qquad
s{\mathrm{lo}} = Q_{0.25}(\bar{s}),
]

[
n_{\mathrm{hi}} = Q_{0.75}(\bar{n}),
\qquad
n_{\mathrm{lo}} = Q_{0.25}(\bar{n}).
]

The practical rule is:

text if sensitivity <= sens_lo:     Fix (irrelevant) elif nullspace >= null_hi:     Fix (anchor) else:     Estimate 

The three classes are:

| Class | Interpretation | Practical decision |
|---|---|---|
| Estimate | Sensitive and sufficiently separable | Include in calibration/profile likelihood |
| Fix (anchor) | Influential but compensatory | Fix, constrain, or anchor using prior knowledge |
| Fix (irrelevant) | Weakly influential in this output setting | Keep fixed unless a different experiment makes it informative |

A parameter can therefore be sensitive but still classified as an anchor if it lies strongly in compensatory directions.

## Key Outputs

The SVD workflow produces:

- identifiability_singular_values.png
- identifiability_nullspace_participation.png
- identifiability_compensation_edges.png
- identifiability_compensation_network.png
- identifiability_decision_map.png

These plots diagnose rank, weak parameter directions, and compensation structure.

For plot-level interpretation, see:

text docs/plot_interpretation_guide.md 

---

# 3. Profile Likelihood Confirmation

Profile likelihood asks:

> Does the model fit become worse when a parameter is moved away from its best value, after allowing other parameters to compensate?

This is a nonlinear practical-identifiability diagnostic.

## Calculation

For observation (z_i), model prediction (m_i(\theta)), and standard deviation (\sigma_i), the weighted least-squares loss is

[
L(\theta)
=
\frac{1}{2}
\sum_i
\left(
\frac{m_i(\theta)-z_i}{\sigma_i}
\right)^2 .
]

For profiled parameter (\theta_j), the parameter is fixed at a grid value

[
\theta_j = a\theta_{j,\mathrm{ref}},
]

and selected nuisance parameters (\eta) are re-optimized:

[
P_j(a)
=
\min_{\eta}
\left[
L(\theta_j=a\theta_{j,\mathrm{ref}},\eta)
+
A(\theta)
\right],
]

where (A(\theta)) is the biological admissibility penalty.

The plotted profile is the increase from the best profile value:

[
\Delta L_j(a)
=
P_j(a)
-
\min_a P_j(a).
]

The profile cutoff used in the plots is

[
c_{\mathrm{profile}} = 1.92,
]

an approximate 95% cutoff for one profiled parameter. Some plots use (\log(1+\Delta L)) for readability, but classification is based on the untransformed loss.

## Profile Classes

| Profile class | Meaning |
|---|---|
| practically identifiable | The profile has a finite optimum and sufficient curvature |
| boundary-limited | The best point lies at the tested range boundary |
| weakly identifiable | The profile has curvature, but a wide acceptable range |
| flat/non-identifiable | Moving the parameter does not worsen the fit enough |

Current curated profile-likelihood summary:

- 60 profiled parameters
- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable: insulin_igf_threshold
- 2 flat/non-identifiable: feed_direct_blood_fraction, lh_basal_release

The interpretation is that many selected parameters are practically estimable under the synthetic-output setup, while some metabolic/endocrine feedback parameters require caution, anchoring, or better experimental design.

---

# 4. Global Sensitivity and Admissible-Bank Association

Global sensitivity screening asks:

> Across a simulation ensemble, which parameters are associated with biomarker variability?

This complements local sensitivity. Local sensitivity measures nominal one-at-a-time response, while global screening summarizes variation across many parameter combinations.

The observable biomarker panel is:

text FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose 

For biomarker (b), each simulation is summarized by an AUC endpoint:

[
\mathrm{AUC}b(\theta)
=
\int{t_0}^{t_f}
y_b(t;\theta),dt.
]

In practice, AUC is evaluated numerically over stored simulation time points using the trapezoidal rule.

## Full-Prior Variance-Based Screening

The unfiltered parameter bank contains ordinary independent Monte Carlo samples. Since the bank was not generated using a Saltelli/Sobol design, this analysis is not reported as strict Sobol sensitivity.

A screening statistic is computed as

[
S^{\mathrm{screen}}{i,b}
=
\frac{
\operatorname{Var}
\left(
\operatorname{E}
[
\mathrm{AUC}b
\mid
\theta_i
]
\right)
}{
\operatorname{Var}
(
\mathrm{AUC}b
)
}.
]

The conditional expectation is approximated using parameter bins.

This statistic gives a variance-based importance screen for the full prior ensemble, but it should not be interpreted as a formal Sobol index.

## Spearman Association

For the biologically admissible bank, Spearman rank correlation is used to measure monotonic association between a parameter and biomarker AUC:

[
\rho^{S}{i,b}
=
\operatorname{Corr}
\left(
\operatorname{rank}(\theta_i),
\operatorname{rank}(\mathrm{AUC}b)
\right).
]

Positive values mean that larger parameter values tend to be associated with larger biomarker AUC. Negative values mean the opposite.

## PRCC Association

Partial rank correlation coefficient (PRCC) measures the association between parameter (\theta_i) and biomarker AUC after accounting for the ranked effects of other parameters.

Let (r{\theta_i}) be the residual from regressing (\operatorname{rank}(\theta_i)) on the ranks of the other parameters. Let (r_b) be the residual from regressing (\operatorname{rank}(\mathrm{AUC}b)) on the ranks of the other parameters. Then

[
\mathrm{PRCC}{i,b}
=
\operatorname{Corr}
\left(
r{\theta_i},
r_b
\right).
]

PRCC and Spearman values indicate direction and strength of association, not formal variance decomposition.

---

# 5. Uncertainty Propagation

Uncertainty propagation asks:

> How much do predicted biomarker trajectories vary across admissible parameter sets?

The current uncertainty analysis uses the biologically admissible ( \pm 0.5% ) Monte Carlo bank. Therefore, it should be interpreted as local robustness around the calibrated model, not full population-level uncertainty.

For biomarker (b), time (t), and admissible parameter samples (\theta_1,\ldots,\theta_N), the ensemble median is

[
\tilde{y}b(t)
=
Q{0.50}
\left(
\left{
y_b(t;\theta_i)
\right}{i=1}^{N}
\right).
]

The lower and upper uncertainty bands are

[
y^{\mathrm{low}}b(t)
=
Q{0.05}
\left(
\left{
y_b(t;\theta_i)
\right}{i=1}^{N}
\right),
]

[
y^{\mathrm{high}}b(t)
=
Q{0.95}
\left(
\left{
y_b(t;\theta_i)
\right}{i=1}^{N}
\right).
]

Figures show:

- shaded region: 5th–95th percentile interval,
- solid line: ensemble median,
- dashed line: nominal trajectory.

Because the ensemble is narrow and biologically filtered, narrow uncertainty bands are expected and should be interpreted as evidence of local robustness near the calibrated parameter regime.

---

# 6. Bayesian Experimental Design

Bayesian experimental design asks:

> Which candidate measurements are expected to reduce uncertainty the most?

It naturally follows sensitivity and identifiability analysis. Classical diagnostics identify weakly informed or compensatory model directions; BED proposes future measurements that could improve them.

## Mutual Information

For target quantity (W) and candidate future measurement (Z), mutual information is

[
I(W;Z)
=
\iint
p(w,z)
\log
\left[
\frac{p(w,z)}{p(w)p(z)}
\right]
,dw,dz.
]

Equivalently,

[
I(W;Z)
=
H(W)
-
H(W\mid Z).
]

A high value means that observing (Z) is expected to strongly reduce uncertainty about (W).

The Monte Carlo estimate is

[
\widehat{I}(W;Z)
=
\frac{1}{N}
\sum_{i=1}^{N}
\log
\left[
\frac{
p(W_i,Z_i)
}{
p(W_i)p(Z_i)
}
\right].
]

The MATLAB BED workflow estimates the required densities using KDE/Gaussian density tools such as ksdensity, mvksdensity, normpdf, and mvnpdf.

## Posterior Updating

For a hypothetical observation (z^\ast), Bayes' rule gives

[
p(w\mid z^\ast)
=
\frac{
p(z^\ast\mid w)p(w)
}{
p(z^\ast)
}.
]

In the Monte Carlo implementation, samples drawn from the prior are reweighted by their likelihood under the synthetic observation.

For sample (\theta_i), model prediction (y_i), observation (z_{\mathrm{obs}}), and observation standard deviation (\sigma),

[
\ell_i
=
-\frac{1}{2}
\sum_j
\left(
\frac{z_{\mathrm{obs},j}-y_{i,j}}{\sigma_j}
\right)^2.
]

Stabilized weights are

[
\tilde{w}i
=
\exp
\left(
\ell_i - \max_k \ell_k
\right),
]

[
w_i
=
\frac{
\tilde{w}i
}{
\sum_k \tilde{w}k
}.
]

The posterior for parameter (\theta_p) is approximated by a weighted kernel density estimate:

[
p(\theta_p\mid z{\mathrm{obs}})
\approx
\sum{i=1}^{N}
w_i
K_h(\theta_p-\theta{i,p}),
]

where (K_h) is a kernel with bandwidth (h).

The posterior narrows when the candidate measurement is informative for the target parameter.

## Current Repository Status

The original PhD BED implementation is retained in MATLAB for methodological provenance. The public repository contains a cleaner v3 baseline BED port using the published MetRep equations with Dexa PK/PD switched off.

Relevant files:

text analyses/bayesian_experimental_design/matlab_original/BED_1M_ALL.m analyses/bayesian_experimental_design/matlab_original/BovSys_run_v3_baseline.m 

The BED workflow should be interpreted as baseline metabolic-reproductive BED, not a Dexa perturbation simulation.

---

# 7. Cross-Method Perturbation Scales

Different perturbation scales are intentional.

| Analysis | Perturbation / ensemble | Purpose |
|---|---|---|
| Local sensitivity | (+1%) one-at-a-time | Local nominal influence |
| SVD identifiability | small central finite difference | Local separability / weak directions |
| Profile likelihood | wider profiling range | Practical nonlinear estimability |
| Global association | Monte Carlo banks | Ensemble-level parameter-output associations |
| Uncertainty propagation | admissible ( \pm 0.5% ) bank | Local robustness around calibrated model |
| BED | prior-based Monte Carlo ensembles | Expected information gain |

The methods are complementary and should not be interpreted as identical sensitivity measures.
