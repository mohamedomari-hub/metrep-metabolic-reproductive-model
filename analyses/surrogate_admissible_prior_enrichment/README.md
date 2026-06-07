# Surrogate-Assisted Admissible Prior Enrichment

This workflow uses the existing 50,000-simulation MetRep `+/-5%` ODE bank to
train:

1. an admissibility classifier; and
2. a biomarker-trajectory/output regression surrogate.

The source bank contains 2,957 biologically admissible simulations, an
acceptance rate of 5.914%. Precision is prioritized because false admissible
candidates are scientifically risky.

## Leakage Control

A single stratified row-index split is shared by classification and regression:
70% train, 15% validation, and 15% test. Scaling, PCA, fitting, and model
selection use training/validation rows only. Test rows are reserved for final
reporting. Split indices are saved to `outputs/split_indices.csv`.

## Scientific Status

The classifier enriches the prior by selecting parameter sets likely to produce
biologically admissible trajectories. The regression surrogate predicts stored
biomarker outputs. The surrogate-window filter applies the original metric
thresholds to days 54-89, but it is not equivalent to exact full-trajectory ODE
admissibility.

A real ODE audit is required before treating the enriched bank as scientifically
reliable. Only ODE-confirmed audit samples are true ODE-admissible. Until audit
precision is acceptable, the downstream bank is exploratory and must be called
a **surrogate-selected enriched candidate bank**.

The initial classifier attempt performed poorly and must not be used for
downstream analysis. The improved classifier runner uses the strict biological
admissibility truth label, explicitly prioritizes precision, tests rare-class
training strategies, and generates candidates only when validation/test results
show useful enrichment. Parameter-only classification may remain intrinsically
difficult when admissibility depends on nonlinear trajectory features.

The preferred workflow is now `run_trajectory_based_enrichment.py`. It predicts
biomarker trajectories, extracts P4 and multi-biomarker biological features,
and applies trajectory rules plus a trajectory-feature classifier. The
parameter-only classifiers must not gate downstream analyses.

The most PhD-faithful refinement is
`run_p4_phd_admissibility_enrichment.py`, which focuses on P4 and applies the
original biological admissibility metrics from
`MetRep_Python/model_definition/admissibility.py`: shifted normalized
cross-correlation, normalized average absolute difference, and normalized
squared-norm difference. Candidate banks from this workflow are still
surrogate-selected and exploratory until exact ODE audit confirms high
precision.

The corrected enrichment framing is
`run_envelope_guided_surrogate_enrichment.py`. It trains a PCA + ExtraTrees
trajectory surrogate, extracts biological trajectory features, builds accepted
feature envelopes from ODE-confirmed admissible training simulations, and ranks
ordinary Monte Carlo candidates by distance from those envelopes. This avoids
treating raw parameter classification or P4-only metrics as final admissibility
truth. Exact ODE audit remains mandatory.

The metric-faithful workflow is
`run_metric_faithful_surrogate_enrichment.py`. Debugging showed the stored
2,957 labels are exactly reproduced by the aggregate admissibility columns:
`penalty`, `min_correlation`, `max_average_difference`, and
`max_norm_difference`. This workflow predicts those metrics directly, applies
the same admissibility rule with conservative margins/top-k ranking, and saves
only ODE-audit candidates for scientific confirmation.

The PGF/E2-weighted workflow is `run_pgf_e2_weighted_ranking.py`. Limiting
biomarker diagnostics showed that PGF and E2 dominate admissibility failures,
while P4-only screening is insufficient for the full bank. PGF/E2 are therefore
used as soft ranking signals, not removed from the biological rule. Final truth
remains exact ODE simulation plus the aggregate admissibility rule.

## Run

```bash
python analyses/surrogate_admissible_prior_enrichment/run_surrogate_admissible_enrichment.py \
  --bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs \
  --n-candidates 1000000 \
  --seed 42
```

Run the exact ODE audit subset separately:

```bash
python analyses/surrogate_admissible_prior_enrichment/run_ode_audit.py \
  --audit-parameters analyses/surrogate_admissible_prior_enrichment/outputs/ode_audit_parameter_samples.csv \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs
```

The first command does not run ODE simulations. The second command runs only
the prepared audit subset.

Improved strict-label classifier:

```bash
python analyses/surrogate_admissible_prior_enrichment/run_improved_admissibility_classifier.py \
  --bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs \
  --n-candidates 1000000 \
  --seed 42
```

Envelope-guided trajectory-feature enrichment:

```bash
python analyses/surrogate_admissible_prior_enrichment/run_envelope_guided_surrogate_enrichment.py \
  --bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs/envelope_guided \
  --n-candidates 1000000 \
  --max-save-selected 100000 \
  --minimum-useful-precision 0.30 \
  --seed 42
```

Metric-faithful aggregate admissibility enrichment:

```bash
python analyses/surrogate_admissible_prior_enrichment/run_metric_faithful_surrogate_enrichment.py \
  --bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs/metric_faithful \
  --n-candidates 1000000 \
  --max-save-selected 100000 \
  --minimum-useful-precision 0.50 \
  --seed 42
```

PGF/E2-weighted ranking for ODE audit:

```bash
python analyses/surrogate_admissible_prior_enrichment/run_pgf_e2_weighted_ranking.py \
  --bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs/pgf_e2_weighted_ranking \
  --n-candidates 1000000 \
  --seed 42
```

Tuned margin ranking for ODE audit:

```bash
python analyses/surrogate_admissible_prior_enrichment/run_tuned_margin_ranking.py \
  --bank-dir analyses/bayesian_experimental_design/surrogate_bed/phd_bed_bank_5pct_50k_glucagon \
  --output-dir analyses/surrogate_admissible_prior_enrichment/outputs/tuned_margin_ranking \
  --n-candidates 1000000 \
  --n-iter 40 \
  --seed 42
```
