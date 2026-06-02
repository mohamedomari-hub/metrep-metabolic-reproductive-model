# Release Checklist

Before publishing this repository on GitHub:

- Initialize the repository and review `git status`.
- Confirm the intended license for code, data, PDFs, and figures.
- Update `CITATION.cff` with the final GitHub URL and publication details.
- Decide whether `MetRep_model.pdf` can be redistributed publicly.
- Keep `MetRep_Python/results/` out of the public commit unless a result is
  intentionally curated into `results_final/`.
- Do not commit `__pycache__`, `.pyc`, `.DS_Store`, or Matplotlib cache files.
- Make Dexa status clear: MATLAB reference completed; Python Dexa extension in
  progress.
- Make BED status clear: use `BED_1M_ALL.m` as the single GitHub-facing MATLAB
  baseline port; keep historical `n2` files out of the repository folder;
  compact Python BED reproduction is future work.
