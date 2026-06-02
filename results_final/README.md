# Curated Results

This folder contains selected result files for public reporting. The full
working output directory is `../MetRep_Python/results/` and may include
exploratory, smoke-test, or intermediate files.

## Tables

- `validation_report.csv`: metadata validation against MATLAB reference values
- `metrep_non_lactating_standard_local_sensitivity_allparams.csv`: local
  sensitivity table
- `sensitivity_top_parameters.csv`: curated top-parameter sensitivity summary
- `structid_50d_measurable_holistic_table.csv`: SVD-based Estimate / Fix table
- `structid_50d_measurable_singular_values.csv`: singular values from the SVD
  identifiability screen
- `structid_50d_measurable_ranking_and_participation.csv`: SVD ranking,
  nullspace participation, and identifiable-score table
- `structid_50d_measurable_compensation_edges.csv`: parameter-pair
  compensation relationships inferred from nullspace directions
- `identifiability_class_counts.csv`: count of Estimate / Fix classes
- `identifiability_fixed_parameters.csv`: parameters fixed by the SVD screen
- `profile_50d_balanced_relaxed_summary.csv`: practical identifiability
  classification from profile likelihood
- `profile_50d_balanced_relaxed_selected_parameters.csv`: profiled and nuisance
  parameter list
- `profile_likelihood_class_counts.csv`: compact profile-class counts

## Figures

- `baseline_selected_states.png`: baseline Python simulation check
- `sensitivity_top_parameters.png`: top local sensitivity parameters
- `identifiability_singular_values.png`: singular-value spectrum of the local
  sensitivity matrix
- `identifiability_svd_ranking.png`: top parameters by local SVD ranking score
- `identifiability_nullspace_participation.png`: parameters most involved in
  weak/non-identifiable directions
- `identifiability_compensation_edges.png`: strongest parameter-pair
  compensation relationships
- `identifiability_compensation_network.png`: graph view of compensation
  relationships, with nodes colored by SVD class
- `identifiability_class_counts.png`: Estimate / Fix class counts
- `identifiability_decision_map.png`: sensitivity vs nullspace decision map
- `profile_50d_balanced_relaxed_combined_profiles.png`: combined profile
  likelihood result
- `profile_50d_balanced_synthetic_outputs.png`: synthetic measurable outputs
  used for profile likelihood
- `bed_MI.png`: BED mutual information by candidate sampling day
- `bed_MI_Individual_Species.png`: BED per-species information
- `bed_Posteriors.png`: posterior comparison
- `bed_Posteriors_Combinations.png`: posterior combinations by species set
- `bed_Time_Comparison_3000_ite.png`: BED timing comparison

## Scope

The published MetRep scenario simulations, Dexa perturbation simulation, and
full BED analysis are primarily reported in the associated paper/thesis
materials. This curated folder focuses on GitHub-facing summaries for the
Python sensitivity and identifiability workflow, plus selected BED figures for
context.
