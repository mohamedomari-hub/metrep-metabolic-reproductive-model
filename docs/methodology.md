# Sensitivity And Identifiability

The analysis workflow uses three complementary stages.

These stages are the bridge between the model and Bayesian experimental design:
they reveal which outputs and parameters are well informed, and which parameter
directions need better experimental measurements.

## 1. Sensitivity Analysis

Sensitivity analysis asks which parameters have the largest local effect on
selected model outputs. It is useful as a first screen, but sensitivity alone is
not identifiability: a parameter can strongly affect outputs and still be hard
to estimate if another parameter can compensate for it.

How it is calculated:

1. Run the model once with the reference parameter values.
2. Change one parameter by a small amount while keeping the others fixed.
3. Run the model again.
4. Compare how much selected outputs change.
5. Repeat this for each parameter.

The result is a local sensitivity score. In simple terms:

```math
S =
\frac{\Delta O}{\Delta \theta}
```

where $O$ is a model output summary and $\theta$ is a parameter.

More specifically, the public sensitivity table uses one-at-a-time relative
sensitivity of output AUC values. For parameter $\theta_j$ and output $y_k(t)$:

```math
\mathrm{AUC}_k(\theta) =
\int y_k(t; \theta)\,dt
```

```math
S_{kj} =
\frac{
  \left(\mathrm{AUC}_k(\theta_j(1+h)) - \mathrm{AUC}_k(\theta_j)\right)
  / \mathrm{AUC}_k(\theta_j)
}{h}
```

where $h = 0.01$ in the sensitivity script by default. This is a forward local
perturbation. Each parameter is changed separately while all other parameters
are kept at their reference values.

The analysis is called local because it tests changes around the reference
parameter set, not across all possible biological values. Outputs can be
summarized by trajectory metrics such as AUC, peak value, mean value, or final
value. Parameters with high scores are influential for the selected outputs.

Main script:

```bash
python MetRep_Python/scripts/04_run_sensitivity.py
```

## 2. SVD Identifiability Screen

The SVD workflow builds a local sensitivity matrix using selected measurable
outputs and analyzes its singular values and nullspace directions.

How it is calculated:

1. Build a sensitivity matrix.
2. In that matrix, rows represent output features and columns represent
   parameters.
3. Each entry says how much one output feature changes when one parameter is
   perturbed.
4. Apply singular value decomposition, or SVD, to the matrix.

SVD separates the sensitivity matrix into informed directions and weak
directions.

The implemented SVD screen uses a stacked trajectory sensitivity matrix. For a
parameter $\theta_j$, selected output vector $Y(\theta)$, and step size $h_j$, the
central-difference column is:

```math
S_{:,j} =
\frac{
  Y(\theta_j + h_j) - Y(\theta_j - h_j)
}{2h_j}
```

```math
h_j =
\max\left(
  \mathrm{absoluteStepMin},
  \mathrm{relativeStep}\cdot |\theta_j|
\right)
```

$Y(\theta)$ is made by stacking all selected output trajectories over the
simulation time points. In the default identifiability script,
`relative_step = 1e-3`. If a nominal parameter value is zero, the code uses
`relative_step * 1.0` before applying the absolute minimum step.

The matrix is then decomposed as:

```math
S = U \Sigma V^T
```

where:

- `S` is the stacked sensitivity matrix.
- `U` contains output-space directions.
- $\Sigma$ contains the singular values.
- $V^T$ contains parameter-space directions.

The numerical rank is calculated using a relative threshold:

```math
\tau = c\,\sigma_1
```

```math
r =
\left|
  \{\sigma_i : \sigma_i > \tau\}
\right|
```

where $\tau$ is the rank threshold, $c$ is the tolerance value, and
$\sigma_1$ is the largest singular value. The default tolerance value is
`1e-8`.

Large singular values indicate parameter combinations that strongly affect the
selected outputs. Very small singular values indicate weak or near-null
directions: parameter combinations that can change while producing little
observable change.

The nullspace part is used to detect compensation. If two parameters appear
together in weak directions, they may compensate for each other. That means the
model can produce similar outputs by changing both parameters together, making
their separate values difficult to estimate from the current measurements.

The nullspace basis is taken from the rows of $V^T$ after the numerical rank:

```math
\mathcal{N} =
\{v_i^T : i > r\}
```

where $r$ is the numerical rank and $\mathcal{N}$ is the nullspace basis.

Nullspace participation for each parameter is summarized as the Euclidean norm
of that parameter's coefficients across all nullspace directions:

```math
n_j =
\sqrt{
  \sum_q \mathcal{N}_{qj}^2
}
```

where $n_j$ is the nullspace participation score for parameter $j$.

The local SVD ranking score used in the curated result is `rel2_colnorm`:

```math
s_j =
\left\|
  \frac{\theta_j}{\max(|Y|,\epsilon)}
  S_{:,j}
\right\|_2
```

where $s_j$ is the local SVD ranking score.

This makes the ranking relative to both parameter size and output scale.

The key SVD outputs are:

- `identifiability_singular_values.png`: shows the singular-value spectrum. A
  wide drop toward small singular values indicates directions in parameter space
  that are weakly informed by the selected outputs.
- `identifiability_nullspace_participation.png`: ranks parameters by how much
  they participate in weak/nullspace directions. High participation means a
  parameter is involved in compensatory combinations.
- `identifiability_compensation_edges.png`: shows strong parameter-pair
  compensation relationships inferred from the nullspace.
- `identifiability_compensation_network.png`: shows the same compensation
  structure as a node-link graph, where nodes are parameters and edges indicate
  compensatory relationships.
- `identifiability_compensation_network_all_parameters.png`: shows all
  analyzed parameters in one zoned network. Strongly compensating parameters
  are central, while parameters with weak or no compensation are still visible
  in peripheral estimate/fix zones.
- `identifiability_decision_map.png`: combines normalized sensitivity and
  nullspace involvement to support the `Estimate`, `Fix (anchor)`, and
  `Fix (irrelevant)` classification.

For a plot-by-plot explanation of the technical terms and conclusions, see
`docs/plot_interpretation_guide.md`.

It classifies parameters into:

- `Estimate`: sensitive and sufficiently separable
- `Fix (anchor)`: high-impact but compensatory; useful as an anchor
- `Fix (irrelevant)`: low-impact in the selected output setting

The three-class decision uses normalized sensitivity and normalized nullspace
scores. The high and low thresholds are quantiles of the analyzed parameter
set:

```math
\bar{s}_j =
\frac{s_j-\min(s)}{\max(s)-\min(s)}
```

```math
\bar{n}_j =
\frac{n_j-\min(n)}{\max(n)-\min(n)}
```

```math
\begin{aligned}
s_{hi} &= Q_{0.75}(\bar{s}) \\
s_{lo} &= Q_{0.25}(\bar{s}) \\
n_{hi} &= Q_{0.75}(\bar{n}) \\
n_{lo} &= Q_{0.25}(\bar{n})
\end{aligned}
```

The practical rule is:

```text
if sensitivity <= sens_lo:
    Fix (irrelevant)
elif nullspace >= null_hi:
    Fix (anchor)
else:
    Estimate
```

This is why a parameter can be sensitive but still fixed as an anchor: it can
affect outputs, but it also sits strongly in compensatory directions.

The class should guide the next modeling step:

| Class | Practical decision |
|---|---|
| `Estimate` | Include in calibration/profile likelihood, because the selected outputs carry enough local information. |
| `Fix (anchor)` | Do not freely estimate together with its compensation partners; fix it, constrain it with prior knowledge, or use it as an anchor. |
| `Fix (irrelevant)` | Keep fixed for this dataset and scenario; it may need a different experiment or output panel to become informative. |

The compensation network is especially useful for the `Fix (anchor)` group. A
node is a parameter, and an edge means two parameters appear together in a weak
direction. Strong edges indicate that the model can preserve similar outputs by
moving those parameters together. This is why compensation is not a failure of
the model; it is a design signal that the current measurements cannot separate
those mechanisms cleanly.

This is a fast local-linear diagnostic and is used to select parameters for
profile likelihood.

## 3. Profile Likelihood Confirmation

Profile likelihood fixes one selected parameter over a grid of values and
allows nuisance parameters to compensate. The resulting loss curve provides a
nonlinear practical-identifiability confirmation.

How it is calculated:

1. Choose one parameter to test.
2. Fix that parameter at a sequence of values, for example 80%, 90%, 100%,
   110%, and 120% of its reference value.
3. At each fixed value, allow selected nuisance parameters to adjust.
4. Re-run the model and calculate the loss, meaning the mismatch between the
   model outputs and the synthetic/reference outputs.
5. Plot the increase in loss relative to the best fit.

The key question is:

```text
Does the fit get clearly worse when this parameter moves away from its best value?
```

If yes, the parameter is practically identifiable. If the curve stays flat, the
parameter is not well identified because other parameters can compensate. If the
curve rises only at one end, the result is boundary-limited and the tested range
may not be wide enough.

The calculation uses synthetic observations from the reference simulation. For
observation $z_i$, model prediction $m_i(\theta)$, and assumed standard
deviation $\sigma_i$, the fit loss is:

```math
L(\theta) =
\frac{1}{2}
\sum_i
\left(
  \frac{m_i(\theta)-z_i}{\sigma_i}
\right)^2
```

The profile for parameter $\theta_j$ fixes $\theta_j$ on a grid of multipliers
and re-optimizes selected nuisance parameters $\eta$:

```math
P_j(a) =
\min_{\eta}
\left[
  L(\theta_j = a\theta_{j,\mathrm{ref}}, \eta)
  + A(\theta)
\right]
```

where $A(\theta)$ is the biological admissibility penalty.

The plotted profile is the increase from the best value:

```math
\Delta L_j(a)
= P_j(a) - \min_a P_j(a)
```

The horizontal cutoff used in the plots is:

```math
c_{\mathrm{profile}} = 1.92
```

This is an approximate 95% cutoff for one profiled parameter. The script can
plot `log1p(delta_loss)` to keep very large curves readable, but the
classification is based on the untransformed loss values.

The class definitions used by the code are:

| Profile class | Code rule | Meaning |
|---|---|---|
| `flat/non-identifiable` | The total delta-loss span is below the cutoff. This rule is checked first. | Moving the parameter across the tested grid does not worsen the fit enough; nuisance parameters and model structure can compensate. |
| `boundary-limited` | The best profile point is at the lowest or highest tested multiplier, after excluding flat profiles. | The optimum may lie outside the tested grid, so the parameter is not safely bounded by this profile range. |
| `weakly identifiable` | At least 80% of the profile grid remains below the cutoff, after excluding flat and boundary-limited profiles. | The profile has some curvature, but the acceptable range is broad. |
| `practically identifiable` | None of the above warning conditions apply. | The tested data/output setup gives a finite optimum for this parameter. |

Current profile-likelihood summary:

- 60 profiled parameters
- 51 practically identifiable
- 6 boundary-limited
- 1 weakly identifiable: `insulin_igf_threshold`
- 2 flat/non-identifiable: `feed_direct_blood_fraction`, `lh_basal_release`

This means that, among the 60 parameters selected from the SVD `Estimate`
group, 51 had a clear enough profile minimum under the synthetic-data setup.
The other 9 were not failures of the model, but they need caution:
boundary-limited parameters need a wider or better-supported profile range,
the weakly identifiable parameter needs more information, and the flat
parameters should not be interpreted as precisely estimable from this output
panel.

The recommended interpretation is that the selected measurable outputs identify
a substantial subset of parameters, while some metabolic/endocrine feedback
parameters require fixing, anchoring, or more informative experimental design.

# Bayesian Experimental Design

The Bayesian experimental design analysis evaluates which sampling times and
measured species provide the most information about model quantities or
parameters.

In this project, BED is used as the experimental-design answer to the
identifiability analysis. Sensitivity and identifiability diagnose which model
directions are weakly informed by existing outputs; BED asks how future
measurements should be chosen to improve those directions.

| Identifiability finding | BED interpretation |
|---|---|
| Sensitive and separable parameters | Current output panel is informative enough; these can be estimated and checked by profile likelihood. |
| High-impact compensation pairs | Future designs should target sampling times/species that separate the paired mechanisms. |
| Strong nullspace participation | Add measurements expected to reduce uncertainty in those weak directions. |
| `Fix (irrelevant)` parameters | Do not spend estimation effort on them unless BED suggests a different output/time window can make them informative. |
| Weak, flat, or boundary-limited profiles | Use BED to propose more informative observations before claiming precise estimates. |

The classical analysis says where the model is under-informed; BED says how a
future experiment could improve that information.

## Methodology

BED asks which future measurements would be most useful before collecting the
data. In this project, the useful measurement is the one expected to give the
most information about a target model quantity or parameter.

The calculation follows this logic:

1. Sample many possible parameter sets around the reference model.
2. Run the model for each sampled parameter set.
3. Store simulated outputs for candidate sampling days and measured species.
4. Treat those simulated outputs as possible future observations.
5. Estimate how much each candidate observation reduces uncertainty about the
   target.
6. Rank sampling days and species by expected information gain.

The main information measure is mutual information.

In simple terms, mutual information measures how much knowing a candidate
measurement tells us about the target.

Mathematically, for a target quantity $W$ and a candidate future measurement
$Z$, mutual information is:

```math
I(W; Z) =
\iint
p(w,z)
\log
\left(
  \frac{p(w,z)}{p(w)p(z)}
\right)
\,dw\,dz
```

Equivalently:

```math
I(W; Z) = H(W) - H(W \mid Z)
```

where $H(W)$ is the uncertainty before observing $Z$, and $H(W \mid Z)$ is the
remaining uncertainty after observing $Z$. Therefore, a high mutual information
value means the candidate measurement is expected to reduce uncertainty about
the target.

In the MATLAB BED script, this is estimated by Monte Carlo simulation and
density estimation:

```text
1. sample parameter sets
2. simulate model outputs for each parameter set
3. form simulated pairs (target W, candidate measurement Z)
4. estimate p(w), p(z), and p(w, z)
5. evaluate log( p(w, z) / (p(w) p(z)) )
6. average/rank this information over candidate designs
```

The script uses MATLAB density functions such as `ksdensity`, `mvksdensity`,
`normpdf`, and `mvnpdf`. Because this is the original PhD MATLAB workflow, the
repository presents it as methodological provenance and selected results, not
as a fully lightweight Python reproduction.

Here:

- the target can be ovulation time or a parameter;
- the measurement can be one species, several species, or a species pair at a
  candidate sampling day;
- higher mutual information means the design is more informative.

## Posterior Calculation By Importance Reweighting

The posterior plots were calculated by importance reweighting of prior samples.
In other words, the workflow first generated parameter samples from the prior,
then used a synthetic observation and a Gaussian observation model to assign a
likelihood weight to each sample.

The steps are:

1. Draw Monte Carlo parameter samples from the prior:

```math
\theta_i \sim p(\theta),
\qquad i=1,\ldots,N
```

In this workflow, the prior was uniform over the selected uncertain
parameters.

2. For each parameter sample, run the model and store the predicted measured
species at the candidate sampling day:

```math
y_i = y(\theta_i)
```

Here, $y_i$ can contain measured species such as FSH, PGF, P4, E2, INH, IGF1,
insulin, and glucose.

3. Generate one fixed synthetic observation vector from the nominal/reference
simulation. If $\mu$ is the nominal model prediction at that sampling day, then:

```math
z_{\mathrm{obs}}
= \mu + \sigma \odot \varepsilon,
\qquad
\varepsilon \sim \mathcal{N}(0,I)
```

The observation standard deviation is defined from a relative noise level:

```math
\sigma_j =
\mathrm{relSigma}\,|\mu_j|
```

with clipping away from zero so that nearly zero outputs do not give a zero
measurement error.

The same fixed observation vector can be used when comparing a full ODE model
and a surrogate calculation, so both workflows are evaluated against the same
synthetic data.

4. Compute the Gaussian likelihood for each prior sample:

```math
p(z_{\mathrm{obs}} \mid \theta_i)
\propto
\exp
\left[
  -\frac{1}{2}
  \sum_j
  \left(
    \frac{z_{\mathrm{obs},j}-y_{i,j}}{\sigma_j}
  \right)^2
\right]
```

Equivalently, the log-likelihood is:

```math
\ell_i =
-\frac{1}{2}
\sum_j
\left(
  \frac{z_{\mathrm{obs},j}-y_{i,j}}{\sigma_j}
\right)^2
```

5. Stabilize and normalize the likelihood weights:

```math
\tilde{w}_i =
\exp(\ell_i - \max_k \ell_k)
```

```math
w_i =
\frac{\tilde{w}_i}{\sum_k \tilde{w}_k}
```

The normalized weights satisfy:

```math
\sum_i w_i = 1
```

These weights measure how compatible each parameter sample is with the
synthetic observation under the assumed Gaussian measurement noise.

6. Estimate the posterior for the parameter of interest, for example
$\theta_p$, using a weighted distribution of the prior samples:

```math
p(\theta_p \mid z_{\mathrm{obs}})
\approx
\sum_i
w_i
K_h(\theta_p - \theta_{i,p})
```

where $K_h$ is a kernel density estimate with bandwidth $h$. Practically, this
means that samples with higher likelihood contribute more strongly to the
posterior density.

This is Bayes' rule written in an importance-sampling form:

```math
p(\theta \mid z_{\mathrm{obs}})
\propto
p(\theta)\,p(z_{\mathrm{obs}}\mid\theta)
```

Because the samples were already drawn from the prior, the likelihood becomes
the weight that reshapes the prior sample cloud into the posterior.

The posterior is narrower than the prior when the synthetic measurement is
informative for the target parameter. If the posterior looks similar to the
prior, that measurement does not strongly reduce uncertainty.

## Mutual Information Calculation

The mutual-information calculation uses the same Monte Carlo idea, but instead
of conditioning on one fixed observed vector, it evaluates how informative a
candidate measurement is on average.

For a candidate sampling day and species set:

1. Use the prior parameter samples to generate paired samples:

```math
(W_i, Z_i)
```

where $W_i$ is the target quantity for sample $i$ and $Z_i$ is the simulated
candidate measurement for the same sample.

2. Estimate the marginal and joint densities from the Monte Carlo cloud:

```math
p(w), \qquad p(z), \qquad p(w,z)
```

In the MATLAB workflow these densities are estimated with KDE/Gaussian density
tools such as `ksdensity`, `mvksdensity`, `normpdf`, and `mvnpdf`.

3. Compute the information contribution:

```math
\log
\left(
  \frac{p(W_i,Z_i)}
       {p(W_i)p(Z_i)}
\right)
```

4. Average this quantity across the Monte Carlo samples:

```math
\widehat{I}(W;Z)
=
\frac{1}{N}
\sum_{i=1}^{N}
\log
\left(
  \frac{p(W_i,Z_i)}
       {p(W_i)p(Z_i)}
\right)
```

A candidate sampling day/species combination receives a high mutual information
score when the simulated measurement $Z$ is strongly informative about the
target $W$. This is why mutual information is used to rank candidate designs
before collecting new data.

## Posterior Interpretation

The posterior is based on Bayes' rule:

```math
p(w \mid z^*) =
\frac{
  p(z^* \mid w)p(w)
}{p(z^*)}
```

where $z^*$ is a hypothetical or selected observation. In practical terms, the
BED result asks whether observing $z^*$ would make the distribution of $W$
narrower or more concentrated than the prior distribution.

This is why BED naturally follows identifiability analysis:

```text
identifiability finds weak or compensatory directions
BED asks which new measurements would reduce those weaknesses
```

## Current Repository Status

The original PhD BED implementation is retained in MATLAB for provenance and methodological traceability. For public presentation, the repository uses the cleaner v3 baseline BED port:

- `analyses/bayesian_experimental_design/matlab_original/BED_1M_ALL.m`
- `analyses/bayesian_experimental_design/matlab_original/BovSys_run_v3_baseline.m`

This version uses the published v3 MetRep equations with Dexa PK/PD switched off, so it should be interpreted as BED for the baseline metabolic-reproductive model, not as a Dexa simulation.

The MATLAB BED workflow is not presented as the lightweight reproducibility path because it combines historical analysis variants, uses MATLAB parallel/toolbox functions, contains hardcoded legacy paths, and does not define a clean public random-seed/output convention.

For reproducibility, the Python surrogate BED workflow is now used as the cleaner public pipeline. It is intended to reproduce the main BED logic from precomputed ODE simulation tables rather than port the original MATLAB file line by line. The full BED result remains reported in the PhD thesis; public BED figures should be interpreted as pilot reproductions unless surrogate validation, mutual-information convergence, and biological admissibility checks are documented.

# Surrogate BED Workflow

The preferred approach is to keep the ODE model as the reference model and use
a surrogate only as an accelerator. The implemented script and detailed
mathematical documentation are in:

`analyses/bayesian_experimental_design/surrogate_bed/`

The surrogate workflow should include:

- held-out ODE validation for the surrogate predictions;
- mutual-information convergence checks across increasing Monte Carlo sample
  sizes;
- biological admissibility filtering before posterior or mutual-information
  summaries are reported;
- repeated-seed or bootstrap uncertainty for candidate ranking stability.

The recommended first surrogate is a PCA-compressed multi-output emulator with
a tree ensemble regressor. This is more suitable than reporting a simple RF
comparison alone because it treats the multi-species output vector as a
correlated object and requires explicit validation before BED conclusions are
claimed.

# Global Sensitivity And Admissible-Bank Association

Global sensitivity screening asks which parameters are associated with
variability across a simulation ensemble rather than only near the nominal
parameter vector. The analysis uses the observable biomarker panel:

FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose.

For biomarker b and parameter vector θ, each simulation
is summarized using the stored-trajectory AUC endpoint:

$$
\mathrm{AUC}b(\theta)
=
\int{t_0}^{t_f}
y_b(t;\theta),dt
$$

This is evaluated numerically using the trapezoidal rule over the stored
sampling days.

## Full-Prior Variance-Based Screening

The unfiltered parameter bank was generated using ordinary independent uniform
Monte Carlo sampling. For each parameter, the screening statistic estimates:

$$
S^{screen}{i,b}
=
\frac{
\mathrm{Var}
\left[
\mathrm{E}
\left(
\mathrm{AUC}b \mid \theta_i
\right)
\right]
}{
\mathrm{Var}
\left(
\mathrm{AUC}b
\right)
}
$$

The conditional mean is approximated with equal-count parameter bins. This is
a variance-based importance screen, but it is not a strict Sobol index:
the simulation bank does not use a Saltelli/Sobol sampling design and therefore
does not support formal Sobol variance decomposition.

## Admissible-Bank Spearman Screening

Spearman rank correlation measures the direction and strength of monotonic
association between a parameter and biomarker AUC:

$$
\rho^S{i,b}
=
\mathrm{Corr}
\left(
\mathrm{rank}(\theta_i),
\mathrm{rank}(\mathrm{AUC}b)
\right)
$$

Positive values indicate that larger parameter values tend to accompany
larger biomarker AUC; negative values indicate the opposite relationship.

## Admissible-Bank PRCC Screening

Partial rank correlation coefficients evaluate the association between
(\theta_i) and (\mathrm{AUC}b) after linearly removing the ranked
effects of the other sampled parameters.

If (r{\theta_i}) and (r{\mathrm{AUC}b}) are residuals from those
rank-based regressions, then:

$$
\mathrm{PRCC}{i,b}
=
\mathrm{Corr}
\left(
r{\theta_i},
r{\mathrm{AUC}_b}
\right)
$$

PRCC values summarize parameter–biomarker AUC associations across biologically
admissible simulations. They indicate direction and strength of monotonic
association, not formal Sobol variance decomposition.

# Uncertainty Propagation

Uncertainty propagation summarizes trajectory variability across the
biologically admissible `+/-0.5%` Monte Carlo bank. For each biomarker $b$
and stored time $t$, the reported ensemble summaries are:

```math
\tilde y_b(t)=Q_{0.50}\{y_b(t;\theta_j)\},
```

```math
y^{low}_b(t)=Q_{0.05}\{y_b(t;\theta_j)\},
\qquad
y^{high}_b(t)=Q_{0.95}\{y_b(t;\theta_j)\}.
```

The figures show the 5th-95th percentile interval as a shaded band, the
ensemble median $\tilde y_b(t)$ as a solid blue line, and the nominal
trajectory $y_b(t;\theta_0)$ as a dashed black line.

Because the ensemble uses a narrow parameter range and biological
admissibility filtering, these bands should be interpreted as local robustness
around the calibrated model rather than full population variability.

# Cross-Method Perturbation Scales

Different perturbation scales are used because each analysis answers a
different question: local sensitivity uses `+1%` one-at-a-time perturbations,
SVD identifiability uses small finite differences, profile likelihood explores
a wider parameter range, global sensitivity uses simulation-bank associations,
uncertainty propagation uses biologically admissible ensembles, and BED uses
prior-based information calculations.

# Appendix

<img width="1169" height="928" alt="image" src="https://github.com/user-attachments/assets/450ae339-29bf-403c-b5f2-a5865417f19d" />


<img width="1326" height="995" alt="image" src="https://github.com/user-attachments/assets/29d46477-d6ec-4d05-baef-975e812575e0" />
