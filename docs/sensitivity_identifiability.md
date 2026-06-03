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

$$
S =
\frac{\Delta O}{\Delta \theta}
$$

where $O$ is a model output summary and $\theta$ is a parameter.

More specifically, the public sensitivity table uses one-at-a-time relative
sensitivity of output AUC values. For parameter $\theta_j$ and output $y_k(t)$:

$$
\mathrm{AUC}_k(\theta) =
\int y_k(t; \theta)\,dt
$$

$$
S_{kj} =
\frac{
  \left(\mathrm{AUC}_k(\theta_j(1+h)) - \mathrm{AUC}_k(\theta_j)\right)
  / \mathrm{AUC}_k(\theta_j)
}{h}
$$

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

$$
S_{:,j} =
\frac{
  Y(\theta_j + h_j) - Y(\theta_j - h_j)
}{2h_j}
$$

$$
h_j =
\max\left(
  \mathrm{absoluteStepMin},
  \mathrm{relativeStep}\cdot |\theta_j|
\right)
$$

$Y(\theta)$ is made by stacking all selected output trajectories over the
simulation time points. In the default identifiability script,
`relative_step = 1e-3`. If a nominal parameter value is zero, the code uses
`relative_step * 1.0` before applying the absolute minimum step.

The matrix is then decomposed as:

$$
S = U \Sigma V^T
$$

where:

- `S` is the stacked sensitivity matrix.
- `U` contains output-space directions.
- $\Sigma$ contains the singular values.
- $V^T$ contains parameter-space directions.

The numerical rank is calculated using a relative threshold:

$$
\tau = c\,\sigma_1
$$

$$
r =
\left|
  \{\sigma_i : \sigma_i > \tau\}
\right|
$$

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

$$
\mathcal{N} =
\{v_i^T : i > r\}
$$

where $r$ is the numerical rank and $\mathcal{N}$ is the nullspace basis.

Nullspace participation for each parameter is summarized as the Euclidean norm
of that parameter's coefficients across all nullspace directions:

$$
n_j =
\sqrt{
  \sum_q \mathcal{N}_{qj}^2
}
$$

where $n_j$ is the nullspace participation score for parameter $j$.

The local SVD ranking score used in the curated result is `rel2_colnorm`:

$$
s_j =
\left\|
  \frac{\theta_j}{\max(|Y|,\epsilon)}
  S_{:,j}
\right\|_2
$$

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

$$
\bar{s}_j =
\frac{s_j-\min(s)}{\max(s)-\min(s)}
$$

$$
\bar{n}_j =
\frac{n_j-\min(n)}{\max(n)-\min(n)}
$$

$$
\begin{aligned}
s_{hi} &= Q_{0.75}(\bar{s}) \\
s_{lo} &= Q_{0.25}(\bar{s}) \\
n_{hi} &= Q_{0.75}(\bar{n}) \\
n_{lo} &= Q_{0.25}(\bar{n})
\end{aligned}
$$

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

$$
L(\theta) =
\frac{1}{2}
\sum_i
\left(
  \frac{m_i(\theta)-z_i}{\sigma_i}
\right)^2
$$

The profile for parameter $\theta_j$ fixes $\theta_j$ on a grid of multipliers
and re-optimizes selected nuisance parameters $\eta$:

$$
P_j(a) =
\min_{\eta}
\left[
  L(\theta_j = a\theta_{j,\mathrm{ref}}, \eta)
  + A(\theta)
\right]
$$

where $A(\theta)$ is the biological admissibility penalty.

The plotted profile is the increase from the best value:

$$
\Delta L_j(a)
= P_j(a) - \min_a P_j(a)
$$

The horizontal cutoff used in the plots is:

$$
c_{\mathrm{profile}} = 1.92
$$

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

## Connection To BED

The identifiability result motivates Bayesian experimental design. Parameters
that are weak, flat, boundary-limited, or compensatory are not simply failures;
they indicate where the current measurement set is not informative enough.

BED addresses the next question:

```text
Given the model and its weakly informed directions,
which sampling times and measured species would add the most information?
```

This makes the workflow constructive: identifiability diagnoses the problem,
and BED proposes how future experiments could improve it.

In practical terms, BED should be aimed at:

- parameters with strong nullspace participation;
- compensation pairs or clusters in the node graph;
- profile-likelihood cases that are weak, flat, or boundary-limited;
- outputs and sampling times likely to distinguish between compensatory
  mechanisms.
