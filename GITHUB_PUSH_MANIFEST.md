# GitHub Push Manifest

This file separates the curated repository content from local/generated files.

## Push These

- `README.md`
- `GITHUB_PUSH_MANIFEST.md`
- `CITATION.cff`
- `pyproject.toml`
- `requirements.txt`
- `.gitignore`
- `MetRep_Python/`
- `MetRep_Matlab/`
- `analyses/`
- `docs/`
- `results_final/`
- `MetRep_model.pdf` if redistribution is allowed

## Keep Local Or Ignored

- `../MetRep_Model_local_archive/`
- `_local_archive/`, if you create one again later
- `MetRep_Python/results/figures/*`
- `MetRep_Python/results/tables/*`
- `MetRep_Python/results/simulations/*`
- `MetRep_Python/.matplotlib/`
- `__pycache__/`
- `*.pyc`

## Result Policy

Use `results_final/` for curated figures and tables that should appear on
GitHub. Use `MetRep_Python/results/` only as a local working-output folder when
running scripts.

## BED Policy

The BED MATLAB reference belongs in
`analyses/bayesian_experimental_design/matlab_original/`.

Use `BED_1M_ALL.m` as the GitHub-facing BED script. It calls
`BovSys_run_v3_baseline.m`, which uses the published v3 model equations with
Dexa disabled. The historical `n2` BED files are kept outside the repository
folder in `../MetRep_Model_local_archive/bed_n2_historical/`.

Do not place BED scripts directly in `MetRep_Matlab/`, because that folder is
reserved for the published MetRep/Dexa v3 MATLAB reference files.
