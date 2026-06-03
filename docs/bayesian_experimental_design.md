# Bayesian Experimental Design

The Bayesian experimental design analysis evaluates which sampling times and
measured species provide the most information about model quantities or
parameters.

In this project, BED is used as the experimental-design answer to the
identifiability analysis. Sensitivity and identifiability diagnose which model
directions are weakly informed by existing outputs; BED asks how future
measurements should be chosen to improve those directions.

## Current Repository Status

The original PhD BED implementation is MATLAB code. It is retained for
provenance and methodological traceability.

For a clearer GitHub presentation, the repository also includes a v3 baseline
port:

- `Bayesian_Experimental_Design/matlab_original/BED_1M_ALL.m`
- `Bayesian_Experimental_Design/matlab_original/BovSys_run_v3_baseline.m`

This port uses the published v3 model equations with Dexa PK/PD switched off.
It should be interpreted as BED on the baseline MetRep model, not as a Dexa
simulation.

The BED workflow is still not the recommended lightweight reproducibility path,
because:

- the original script is large and combines several analysis variants
- it uses MATLAB parallel loops
- it depends on Statistics and Machine Learning Toolbox functions
- it contains hardcoded historical cluster paths
- it saves intermediate `.mat` files in the working directory
- exact repeatability is limited because the original script does not define a
  clean public random-seed/output convention
- the existing BED figures used a surrogate comparison that should not be
  reported publicly until surrogate validation, convergence, and biological
  admissibility checks are documented

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

$$
I(W; Z) =
\iint
p(w,z)
\log
\left(
  \frac{p(w,z)}{p(w)p(z)}
\right)
\,dw\,dz
$$

Equivalently:

$$
I(W; Z) = H(W) - H(W \mid Z)
$$

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

$$
\theta_i \sim p(\theta),
\qquad i=1,\ldots,N
$$

In this workflow, the prior was uniform over the selected uncertain
parameters.

2. For each parameter sample, run the model and store the predicted measured
species at the candidate sampling day:

$$
y_i = y(\theta_i)
$$

Here, $y_i$ can contain measured species such as FSH, PGF, P4, E2, INH, IGF1,
insulin, and glucose.

3. Generate one fixed synthetic observation vector from the nominal/reference
simulation. If $\mu$ is the nominal model prediction at that sampling day, then:

$$
z_{\mathrm{obs}}
= \mu + \sigma \odot \varepsilon,
\qquad
\varepsilon \sim \mathcal{N}(0,I)
$$

The observation standard deviation is defined from a relative noise level:

$$
\sigma_j =
\mathrm{relSigma}\,|\mu_j|
$$

with clipping away from zero so that nearly zero outputs do not give a zero
measurement error.

The same fixed observation vector can be used when comparing a full ODE model
and a surrogate calculation, so both workflows are evaluated against the same
synthetic data.

4. Compute the Gaussian likelihood for each prior sample:

$$
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
$$

Equivalently, the log-likelihood is:

$$
\ell_i =
-\frac{1}{2}
\sum_j
\left(
  \frac{z_{\mathrm{obs},j}-y_{i,j}}{\sigma_j}
\right)^2
$$

5. Stabilize and normalize the likelihood weights:

$$
\tilde{w}_i =
\exp(\ell_i - \max_k \ell_k)
$$

$$
w_i =
\frac{\tilde{w}_i}{\sum_k \tilde{w}_k}
$$

The normalized weights satisfy:

$$
\sum_i w_i = 1
$$

These weights measure how compatible each parameter sample is with the
synthetic observation under the assumed Gaussian measurement noise.

6. Estimate the posterior for the parameter of interest, for example
$\theta_p$, using a weighted distribution of the prior samples:

$$
p(\theta_p \mid z_{\mathrm{obs}})
\approx
\sum_i
w_i
K_h(\theta_p - \theta_{i,p})
$$

where $K_h$ is a kernel density estimate with bandwidth $h$. Practically, this
means that samples with higher likelihood contribute more strongly to the
posterior density.

This is Bayes' rule written in an importance-sampling form:

$$
p(\theta \mid z_{\mathrm{obs}})
\propto
p(\theta)\,p(z_{\mathrm{obs}}\mid\theta)
$$

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

$$
(W_i, Z_i)
$$

where $W_i$ is the target quantity for sample $i$ and $Z_i$ is the simulated
candidate measurement for the same sample.

2. Estimate the marginal and joint densities from the Monte Carlo cloud:

$$
p(w), \qquad p(z), \qquad p(w,z)
$$

In the MATLAB workflow these densities are estimated with KDE/Gaussian density
tools such as `ksdensity`, `mvksdensity`, `normpdf`, and `mvnpdf`.

3. Compute the information contribution:

$$
\log
\left(
  \frac{p(W_i,Z_i)}
       {p(W_i)p(Z_i)}
\right)
$$

4. Average this quantity across the Monte Carlo samples:

$$
\widehat{I}(W;Z)
=
\frac{1}{N}
\sum_{i=1}^{N}
\log
\left(
  \frac{p(W_i,Z_i)}
       {p(W_i)p(Z_i)}
\right)
$$

A candidate sampling day/species combination receives a high mutual information
score when the simulated measurement $Z$ is strongly informative about the
target $W$. This is why mutual information is used to rank candidate designs
before collecting new data.

## Posterior Interpretation

The posterior is based on Bayes' rule:

$$
p(w \mid z^*) =
\frac{
  p(z^* \mid w)p(w)
}{p(z^*)}
$$

where $z^*$ is a hypothetical or selected observation. In practical terms, the
BED result asks whether observing $z^*$ would make the distribution of $W$
narrower or more concentrated than the prior distribution.

This is why BED naturally follows identifiability analysis:

```text
identifiability finds weak or compensatory directions
BED asks which new measurements would reduce those weaknesses
```

## Recommended Public Presentation

Present v3 as the canonical model line. Use the single clean MATLAB BED script
`Bayesian_Experimental_Design/matlab_original/BED_1M_ALL.m`, which is
the v3 baseline port.

The full BED result is reported in the PhD thesis. The current public GitHub
repository does not include BED result figures in `results_final` because the
available thesis-style figures rely on an ODE-versus-surrogate comparison that
should be rechecked with a documented public validation workflow.

The Python surrogate workflow is now implemented as a reproducible comparison
pipeline. It is intended to reproduce the main BED story from precomputed ODE
tables rather than port the original MATLAB file line by line.

## Surrogate BED Workflow

The preferred approach is to keep the ODE model as the reference model and use
a surrogate only as an accelerator. The implemented script and detailed
mathematical documentation are in:

`Bayesian_Experimental_Design/surrogate_bed/`

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

## Interpretation In The Project

The intended message is:

```text
Classical identifiability analysis shows the current measurement set does not
equally inform all model parameters.

Bayesian experimental design ranks candidate sampling times and measured
species by expected information gain, suggesting how future experiments can be
made more informative.
```

This connects the identifiability and BED parts into one workflow rather than
two unrelated analyses.

## Link Back To Identifiability

The BED section should be presented as the answer to the specific weaknesses
found by sensitivity and identifiability:

| Identifiability finding | BED interpretation |
|---|---|
| Sensitive and separable parameters | Current output panel is informative enough; these can be estimated and checked by profile likelihood. |
| High-impact compensation pairs | Future designs should target sampling times/species that separate the paired mechanisms. |
| Strong nullspace participation | Add measurements expected to reduce uncertainty in those weak directions. |
| `Fix (irrelevant)` parameters | Do not spend estimation effort on them unless BED suggests a different output/time window can make them informative. |
| Weak, flat, or boundary-limited profiles | Use BED to propose more informative observations before claiming precise estimates. |

This framing keeps BED from looking like an unrelated extra analysis. The
classical analysis says where the model is under-informed; BED says how a
future experiment could improve that information.
