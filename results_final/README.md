# Curated Results

This folder contains selected result files for public reporting.

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
- `combined_parameter_summary.csv`: representative parameters compared across
  local sensitivity, identifiability class, PRCC, and Spearman association
- `global_sensitivity_spearman_auc.csv`: admissible-bank rank associations
- `global_sensitivity_prcc_auc.csv`: admissible-bank partial rank correlations
- `top_parameters_by_biomarker_auc.csv`: strongest global associations per
  observable biomarker
- `uncertainty_summary.csv`: compact uncertainty-band summary
- `trajectory_quantiles.csv`: median and 5th-95th percentile trajectories

## Figures

- `baseline_selected_states.png`: baseline Python simulation check for the
  main metabolic and reproductive states
- `baseline_all_states.png`: baseline Python simulation check for all 22 model
  states
- `sensitivity_top_parameters.png`: top local sensitivity parameters
- `identifiability_singular_values.png`: singular-value spectrum of the local
  sensitivity matrix
- `identifiability_svd_ranking.png`: top parameters by local SVD ranking score
- `identifiability_nullspace_participation.png`: parameters most involved in
  weak/non-identifiable directions
- `identifiability_nullspace_participation_all_parameters.png`: all analyzed
  parameters ranked by weak/nullspace participation
- `identifiability_sensitivity_vs_nullspace_all_parameters.png`: all-parameter
  map comparing sensitivity against nullspace participation
- `identifiability_compensation_edges.png`: strongest parameter-pair
  compensation relationships
- `identifiability_compensation_network_core.png`: readable core network of
  the strongest compensation relationships; node color shows recommendation
  class, node size shows sensitivity, and edge width shows compensation
  strength
- `identifiability_parameter_scenario_map.png`: all-parameter diagnostic map
  showing sensitivity versus nullspace/compensation involvement
- `identifiability_sensitivity_ranked_by_class.png`: all analyzed parameters
  ranked by sensitivity and colored by estimate/fix recommendation class
- `identifiability_compensation_network.png`: older graph view of
  compensation relationships, with nodes colored by SVD class
- `identifiability_compensation_network_sensitivity.png`: compensation network
  with node size and outline indicating local sensitivity
- `identifiability_compensation_network_all_parameters.png`: all analyzed
  parameters shown in a zoned compensation network; strongest compensation
  edges are central and no-edge parameters are placed in peripheral groups
- `identifiability_class_counts.png`: Estimate / Fix class counts
- `identifiability_decision_map.png`: sensitivity vs nullspace decision map
- `profile_likelihood_representative_3x3.png`: representative profile
  likelihood examples for practically identifiable, boundary-limited, weakly
  identifiable, and flat/non-identifiable classes
- `profile_50d_balanced_synthetic_outputs.png`: synthetic measurable outputs
  used for profile likelihood
- `combined_parameter_diagnostics.png`: combined local sensitivity, global
  association, and identifiability diagnostic summary
- `global_sensitivity_representative_identifiability_parameters.png`:
  representative-parameter PRCC heatmap
- `global_sensitivity_heatmap_top_parameters.png`: supplementary strongest
  global associations
- `uncertainty_readme_summary.png`: compact uncertainty propagation summary
- `uncertainty_reproductive.png`: reproductive-biomarker uncertainty bands
- `uncertainty_metabolic.png`: metabolic-biomarker uncertainty bands

## Scope

BED scripts and generated outputs remain under
`analyses/bayesian_experimental_design/` while the expanded prior-based result
is being regenerated. This folder contains selected stable GitHub-facing
outputs only.
