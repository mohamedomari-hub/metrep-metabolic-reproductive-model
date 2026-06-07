# Methodology

This document describes the final public methodology for the MetRep
metabolic-reproductive model portfolio. The analyses are complementary: local
sensitivity identifies influential nominal mechanisms, admissibility filtering
defines physiologically plausible ensembles, global sensitivity links
parameters to biomarkers, identifiability diagnostics classify learnability,
uncertainty propagation identifies informative windows, and Bayesian
experimental design tests whether selected observations reduce parameter
uncertainty.

## 1. Local Sensitivity Analysis

Local sensitivity quantifies the one-at-a-time response of model outputs to a
small perturbation of each parameter near the calibrated parameter vector.
For parameter $\theta_j$, the perturbed value is:

```math
\theta_j^{+}
=
\theta_j(1+\delta),
\qquad
\delta = 0.01 .
```

For biomarker $y_k(t)$, the scalar endpoint is the trajectory AUC:

```math
\mathrm{AUC}_k(\theta)
=
\int_{t_0}^{t_f} y_k(t;\theta)\,dt .
```

The reported local sensitivity score is a normalized AUC response:

```math
S_{j,k}
=
\frac{
  \mathrm{AUC}_k(\theta_j^{+})-\mathrm{AUC}_k(\theta)
}{
  \delta\,\mathrm{AUC}_k(\theta)
}.
```

Local sensitivity is a nominal diagnostic. It identifies parameters that can
move outputs near the calibrated regime, but it does not prove practical
identifiability because other parameters may compensate for the same output
change.

The practical reason for using AUC is that many model outputs are trajectories
rather than single endpoint measurements. A parameter may change timing,
amplitude, or duration of a hormone/metabolite response. Integrating the
trajectory gives one comparable scalar per biomarker while retaining the fact
that the model is dynamic. The local sensitivity table therefore answers:
"If one parameter is nudged slightly, which observable trajectories change the
most near the calibrated solution?"

This analysis is deliberately one-at-a-time. It does not sample all parameters
together, and it does not account for parameter interactions. Its role in the
portfolio is to identify candidate influential mechanisms before moving to
global sensitivity, compensation, and identifiability analyses.

## 2. Biological Admissibility Filtering

Monte Carlo parameter banks can generate unrealistic endocrine or metabolic
trajectories even when parameter values remain close to the nominal model.
Before global sensitivity, uncertainty propagation, and BED ranking, simulated
ODE trajectories are filtered for biological admissibility.

The historical admissibility panel is:

```text
FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose
```

Glucagon was excluded from the biological admissibility filter but retained as
an observable biomarker for downstream uncertainty propagation, global
sensitivity, and Bayesian experimental design.

The final ODE-confirmed admissibility rule uses finite outputs and stored
trajectory diagnostics:

```text
finite outputs
AND penalty <= 0
AND min_correlation >= 0.75
AND max_average_difference <= 0.30
AND max_norm_difference <= 0.30
```

These terms have the following meaning:

- `finite outputs`: every simulated value used by the filter must be a real
  finite number. Simulations with `NaN`, `Inf`, solver failure, or numerical
  blow-up are rejected before biological scoring.
- `penalty`: an aggregate violation score from the historical biological
  filter. A value `<= 0` means that no stored hard biological rule was
  violated. A positive value means that at least one trajectory-level rule was
  outside the allowed biological envelope.
- `min_correlation`: the minimum, across the admissibility biomarkers, of the
  best lag-aligned normalized correlation between a simulated trajectory and
  the reference trajectory shape. The threshold `>= 0.75` requires every
  filtered biomarker to preserve the expected temporal pattern reasonably well.
- `max_average_difference`: the largest average normalized trajectory
  difference across the admissibility biomarkers. The threshold `<= 0.30`
  rejects simulations whose average biomarker levels deviate too far from the
  reference behavior.
- `max_norm_difference`: the largest normalized trajectory-distance score
  across the admissibility biomarkers. The threshold `<= 0.30` rejects
  simulations with excessive whole-trajectory deviation even if pointwise
  averages remain acceptable.

In short, the rule keeps only simulations that are numerically valid, satisfy
hard biological constraints, preserve the reference trajectory shape, and stay
within a normalized amplitude/distance envelope for the historical observable
panel.

Only ODE-confirmed rows are treated as scientific truth. Surrogate or ML models
may rank proposals for enrichment, but they do not define admissibility unless
the ODE trajectory is subsequently confirmed by the exact rule.

This distinction is central. The admissible bank is not simply a random subset
of parameter draws; it is the subset whose ODE solutions remain biologically
reasonable under the historical trajectory-shape and amplitude rules. The
downstream ensemble analyses therefore describe variability and associations
inside the biologically plausible region of the model, not across every
mathematically possible trajectory produced by the prior.

## 3. Global Sensitivity And Admissible-Bank Association

Global sensitivity is evaluated on saved ODE simulation banks. The final
admissible-bank analysis uses the enriched ODE-confirmed ensemble and all
observable biomarkers:

```text
FSH, PGF, P4, E2, INH, IGF1, Insulin, Glucose, Glucagon
```

For parameter $\theta_j$ and biomarker $y_k(t)$, the scalar endpoint is:

```math
\mathrm{AUC}_{i,k}
=
\int_{t_0}^{t_f} y_k(t;\theta_i)\,dt .
```

Spearman association is the Pearson correlation of ranks:

```math
\rho^{S}_{j,k}
=
\mathrm{corr}
\left(
  \mathrm{rank}(\theta_{\cdot j}),
  \mathrm{rank}(\mathrm{AUC}_{\cdot k})
\right).
```

PRCC removes the linear rank effects of all other sampled parameters from both
the target parameter and biomarker endpoint, then correlates the residuals:

```math
\mathrm{PRCC}_{j,k}
=
\mathrm{corr}
\left(
  r_{\theta_j},
  r_{\mathrm{AUC}_k}
\right).
```

The final output is a `98 x 9` parameter-biomarker matrix for PRCC and
Spearman. These values summarize monotonic parameter-biomarker AUC association
across ODE-confirmed admissible simulations. They are not strict Sobol indices:
the banks use ordinary Monte Carlo and SMC-enriched admissible sampling rather
than a Saltelli/Sobol design.

PRCC values should be interpreted as monotonic association measures. They are
useful for identifying parameter-biomarker links inside the admissible
ensemble, but they do not capture all possible nonlinear or non-monotonic
relationships. Strong PRCC/Spearman links are therefore used as interpretable
screening signals, not as formal causal effects or Sobol variance
decompositions.

This analysis answers a different question from local sensitivity. Local
sensitivity asks how outputs respond to a small one-at-a-time perturbation near
the nominal parameter vector. Global PRCC/Spearman analysis asks whether, after
sampling many biologically admissible ODE solutions, larger or smaller values
of a parameter are consistently associated with larger or smaller biomarker AUC
values. A parameter can be locally influential but weak globally if its effect
is compensated by other parameters, or if the relationship changes across the
admissible region.

The PRCC step is especially useful because each parameter is correlated with a
biomarker after removing rank-linear effects of the remaining parameters. In
plain terms, PRCC asks whether a parameter still carries a monotonic signal for
a biomarker after accounting for the rest of the sampled parameter vector.
Spearman is more direct and less adjusted. The two are interpreted together:
agreement strengthens confidence in a robust monotonic link, while
disagreement suggests compensation, confounding, or a more complex
relationship.

## 4. SVD Identifiability Screen

The SVD identifiability screen linearizes selected model outputs with respect
to parameters around the nominal parameter vector. Let $J$ be the sensitivity
matrix:

```math
J_{m,j}
=
\frac{\partial y_m}{\partial \theta_j}.
```

The singular value decomposition is:

```math
J
=
U\Sigma V^T .
```

Large singular values indicate parameter directions that are well expressed in
the measured outputs. Small singular values indicate weak directions. The
right-singular vectors identify combinations of parameters that can compensate
for each other.

SVD is a local structural/practical screen. It helps classify parameters into
estimate, anchor, or fix candidates before nonlinear profile likelihood is
used for confirmation.

The SVD screen is useful because identifiable parameters are not judged one by
one in isolation. The model output may be sensitive to a combination of
parameters while being unable to separate the individual parameters inside that
combination. Small singular values indicate output directions where changes in
parameter combinations produce little observable change. Parameters with strong
participation in those weak right-singular vectors are candidates for weak
identifiability or compensation.

The compensation network is derived from these weak directions. Edges indicate
parameters that repeatedly appear together in poorly informed combinations.
Such pairs or groups should be handled carefully in downstream fitting because
estimating them together may produce unstable or non-unique parameter values
without stronger data or priors.

## 5. Profile Likelihood Confirmation

Profile likelihood fixes one selected parameter over a grid of values and
allows nuisance parameters to compensate. The resulting loss curve provides a
nonlinear practical-identifiability check.

For fixed parameter value $\theta_j=c$, the profiled objective is:

```math
\mathcal{L}_{p,j}(c)
=
\min_{\theta_{-j}}
\mathcal{L}
\left(
  \theta_j=c,\theta_{-j}
\right).
```

The normalized profile is reported relative to the best fit:

```math
\Delta \mathcal{L}_{p,j}(c)
=
\mathcal{L}_{p,j}(c)
-
\min_c \mathcal{L}_{p,j}(c).
```

Representative parameters are classified as practically identifiable,
boundary-limited, weakly identifiable, or flat/non-identifiable based on the
shape of the profile. The fixed 3 x 3 representative set is used consistently
for BED and Bayesian posterior comparisons.

Profile likelihood is the nonlinear follow-up to the SVD screen. For each
profile point, one parameter is fixed and the remaining nuisance parameters are
allowed to re-adjust. This tests whether the data/output scenario truly
requires a specific value of the fixed parameter, or whether other parameters
can compensate. A sharply curved profile indicates that moving the parameter
away from the optimum worsens the fit even after compensation. A flat or broad
profile indicates that the parameter can move without much loss of fit.

The classes used in the final 3 x 3 panel have the following interpretation:

- Practically identifiable: the profile has a clear minimum and rises on both
  sides, so the selected outputs constrain the parameter.
- Boundary-limited: the profile improves or remains acceptable toward a
  boundary, so the parameter is only partly constrained in the explored range.
- Weakly identifiable: the profile has some structure but does not produce a
  strong, well-bounded minimum.
- Flat/non-identifiable: the profile remains shallow, meaning the selected
  outputs do not meaningfully constrain the parameter.

These labels are used as the reference diagnosis when interpreting BED and
Bayesian posterior narrowing. A posterior that remains broad for a weak/flat
profile is not a failure of BED; it is consistent with the identifiability
diagnosis.

## 6. Uncertainty Propagation

Uncertainty propagation summarizes the ODE-confirmed admissible ensemble at
each time point and observable biomarker. For biomarker $y_k$ at time $t$, the
reported curves are:

```math
q_{0.05,k}(t)
=
Q_{0.05}
\left(
  \{y_k(t;\theta_i)\}_{i=1}^{N}
\right),
```

```math
\tilde{y}_k(t)
=
Q_{0.50}
\left(
  \{y_k(t;\theta_i)\}_{i=1}^{N}
\right),
```

```math
q_{0.95,k}(t)
=
Q_{0.95}
\left(
  \{y_k(t;\theta_i)\}_{i=1}^{N}
\right).
```

Figures show the 5th-95th percentile band, the ensemble median, and the nominal
trajectory. Informative windows are times where admissible trajectories remain
biologically plausible but meaningfully separated.

Uncertainty propagation is not used here as a population variability estimate.
It is an ensemble robustness analysis around the calibrated regime, restricted
to ODE-confirmed admissible simulations. The 5th-95th percentile band shows
how much model trajectories can vary while still satisfying the biological
admissibility filter. The nominal dashed trajectory anchors the ensemble to
the calibrated baseline.

The uncertainty windows are important for experimental design. A biomarker can
be globally associated with a parameter, but if the ensemble trajectories are
nearly indistinguishable at a candidate day, that day is unlikely to be useful
for learning. BED therefore benefits from combining parameter-biomarker links
from global sensitivity with time windows where trajectories separate.

## 7. Bayesian Experimental Design

The thesis BED target is mutual information between uncertain parameters and
future measurements. The primary global design quantity is:

```math
I(\Theta;Y)
=
\mathbb{E}
\left[
  \log
  \frac{p(\Theta,Y)}
       {p(\Theta)\,p(Y)}
\right].
```

The Monte Carlo estimator follows the MATLAB thesis workflow:

```math
\widehat{I}(\Theta;Y)
=
\frac{1}{N}
\sum_{i=1}^{N}
\log
\left(
  \frac{
    \widehat{p}(\Theta_i,Y_i)
  }{
    \widehat{p}(\Theta_i)\,\widehat{p}(Y_i)
  }
\right).
```

Posterior updates use the broad `+/-5%` prior bank. For a selected synthetic
observation $z^*$, the normalized posterior for a representative parameter is:

```math
p(\theta_j \mid z^*)
=
\frac{
  p(\theta_j,z^*)
}{
  \int p(\theta_j,z^*)\,d\theta_j
}.
```

The final BED portfolio includes:

1. independent observation scenarios;
2. global cumulative biomarker acquisition from best 1 through best 9
   biomarkers;
3. highest- versus lowest-information day comparisons;
4. parameter-specific GSA + uncertainty + MI guided updates.

The global cumulative biomarker update and the parameter-specific guided update
are intentionally different analyses. The cumulative update uses a global
acquisition order to show how adding biomarkers from 1 to 9 changes the
posterior. The guided update instead chooses a separate biomarker-day scenario
for each representative parameter using GSA, uncertainty windows, and MI
ranking.

For a representative parameter $\theta_j$, candidate biomarkers are selected
from the strongest admissible-bank sensitivity links:

```math
B_j
=
\operatorname{TopK}
\left(
  |\mathrm{PRCC}_{j,k}|
\right).
```

Candidate measurement days are restricted to high-uncertainty windows:

```math
T_k
=
\left\{
  t :
  q_{0.95,k}(t) - q_{0.05,k}(t) > \tau_k
\right\}.
```

Mutual information is then evaluated for candidate biomarker-day pairs:

```math
I_j(t,k)
=
I\left(\theta_j;Y_k(t)\right).
```

The guided observation scenario is selected as:

```math
(t^*,k^*)
=
\arg\max_{\substack{t \in T_k \\ k \in B_j}}
I_j(t,k).
```

This links profile-likelihood class, global sensitivity, uncertainty
propagation, and BED into one parameter-specific observation-selection rule.
The SMC+ML enriched 12,721-row admissible bank is used for stable MI ranking
and robustness checks. It is not used as the plotted prior distribution.

The BED implementation has two levels:

1. Global full-vector BED: estimates information between the full parameter
   vector and all selected observable measurements. This is the faithful thesis
   target because the original PhD workflow asked which measurements are
   informative about the uncertain model as a whole.
2. Parameter-specific BED: estimates information between one representative
   parameter and candidate biomarker-day observations. This is secondary and
   interpretive. It explains why particular biomarkers or days are useful for
   the profile-likelihood representative parameters.

Measurement noise follows the thesis convention: synthetic observations are
generated from ODE outputs with Gaussian noise scaled to approximately 10% of
the mean output level for each biomarker. This prevents the design from
assuming unrealistically perfect measurements.

Posterior plots are generated from the broad `+/-5%` prior bank so the
visualized prior remains the scientific prior. The enriched admissible bank is
used to stabilize ranking, coverage, and MI estimates; it is not treated as the
prior in the posterior plots.

## 8. Bayesian Inference

The final Bayesian inference layer uses saved broad-prior ODE archive rows.
Surrogate models are not used for the final posterior distributions.

The Bayesian inference layer is designed as a consistency check after BED, not
as a replacement for the ODE model. The selected BED scenarios define
synthetic observations such as "P4 at day 81" or parameter-specific
biomarker-day combinations. The inference scripts then ask whether those
observations concentrate probability over the same representative parameters
that the identifiability and BED analyses highlighted.

All final posterior distributions are archive-based. That means the candidate
parameter values are rows from an existing broad `+/-5%` ODE simulation archive.
The scripts evaluate how compatible each archived trajectory is with the
selected synthetic observations and then reweight or filter the archive. No
surrogate-predicted parameter sets are used as scientific posterior samples.

The reduced archive-based posterior uses Gaussian observation weights on
precomputed ODE rows:

```math
\log p(y^{obs}\mid \theta_i)
=
-\frac{1}{2}
\sum_k
\left(
  \frac{y^{obs}_k-y_{i,k}}{\sigma_k}
\right)^2 .
```

The normalized archive weights are:

```math
w_i
=
\frac{p(y^{obs}\mid \theta_i)}
{\sum_j p(y^{obs}\mid \theta_j)} .
```

This is a reduced archive-based posterior over simulated rows. It should not be
described as live ODE MCMC, because no new ODE simulations are proposed or
accepted during the update. The method is useful as a transparent Bayesian
posterior approximation over a trusted ODE simulation archive.

The observation scale $\sigma_k$ controls how strongly a mismatch between the
synthetic observation and archived trajectory is penalized. Smaller values
produce sharper posteriors because only archive rows close to the observation
receive substantial weight. Larger values produce broader posteriors because
more rows remain plausible. In the final workflow this scale is tied to the
same measurement-error logic used by BED, so posterior narrowing reflects the
assumed measurement precision.

The plotted one-dimensional posterior for parameter $\theta_j$ is obtained by
applying the normalized archive weights to the sampled values of
$\theta_{i,j}$. KDE curves are used only for visualization; the underlying
posterior approximation is the weighted set of ODE-confirmed archive rows.

The reduced archive posterior is useful because it is transparent: every
posterior particle corresponds to an actual ODE simulation. Its limitation is
coverage. If the broad archive is too sparse near the true high-probability
region, the posterior can look rough or multimodal because it is limited by
the archive rows available.

Archive-based sequential ABC filtering also uses only precomputed broad-prior
ODE rows. At round $r$, a row is retained when:

```math
\rho(S(y_i),S(y^{obs})) \le \epsilon_r .
```

The tolerance $\epsilon_r$ is reduced across rounds using archive distance
quantiles. This is archive-based sequential ABC filtering, not full adaptive
ABC-SMC: particles are not perturbed and ODEs are not rerun.

Therefore, this analysis should be reported as archive-based sequential ABC
filtering, not full adaptive ABC-SMC. It is a likelihood-free posterior
approximation over precomputed ODE rows. Its strength is computational
efficiency and use of trusted simulations; its limitation is that it cannot
discover posterior regions that were not represented in the original archive.

The ABC distance is computed on selected summary observations rather than a
Gaussian likelihood. At early rounds, a loose tolerance retains many archive
rows. At later rounds, the tolerance is reduced so only rows with smaller
distance to the synthetic observation remain. More rounds can sharpen the
posterior, but only until the archive resolution becomes the limiting factor.
If a parameter remains broad after many tolerance reductions, the likely
interpretation is weak information in the selected observations or insufficient
archive support, not necessarily an implementation error.

The final Bayesian interpretation uses all archive-based methods together:

- Reweighting gives a likelihood-style posterior over existing ODE rows.
- Archive ABC filtering gives a likelihood-free tolerance-based posterior over
  existing ODE rows.
- BED posterior plots show how specific biomarker-day choices reshape the
  broad prior.

Agreement across these methods supports a robust interpretation. Disagreement
is treated as diagnostic evidence that the posterior is sensitive to the
chosen observation scenario, noise scale, distance metric, or archive coverage.

## 9. Cross-Method Perturbation Scales

Different perturbation scales are used because each analysis answers a
different question:

- Local sensitivity uses `+1%` one-at-a-time perturbations.
- SVD identifiability uses small finite differences.
- Profile likelihood explores a wider parameter range.
- Global sensitivity uses simulation-bank associations.
- Uncertainty propagation uses biologically admissible ensembles.
- BED uses prior-based information calculations.

These scales should not be interpreted as contradictory. They reflect distinct
diagnostic questions.

## 10. Final Integrated Interpretation

The final workflow combines classical identifiability diagnostics with
Bayesian experimental design. Local sensitivity and SVD identify influential
and compensatory directions. Profile likelihood classifies representative
parameters as practically identifiable, boundary-limited, weakly identifiable,
or flat. Biological admissibility filtering restricts ensemble analyses to
physiologically plausible ODE trajectories. Global sensitivity and uncertainty
propagation then define candidate biomarkers and time windows for BED.

Bayesian updating confirms the identifiability diagnosis. Practically
identifiable parameters narrow strongly, boundary-limited parameters show
partial or one-sided learning, and weak/flat parameters remain broad even under
guided observations. BED improves learning where information exists, but it
does not magically rescue parameters that are practically non-identifiable
under the available output panel.
