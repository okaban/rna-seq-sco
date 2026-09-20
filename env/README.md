# Environment pins for the analysis (epi-rna)

Captured 2026-09-20 from the environment that produced every table and figure in
the R4 manuscript. Two equivalent forms:

- `epi-rna.requirements.lock.txt` — exact `pip freeze` of the live environment (authoritative).
- `environment.yml` — the same pins arranged as a conda spec (compiled packages via conda-forge/bioconda, the rest via pip).
- `python_version.txt` — interpreter and platform.

Recreate: `conda env create -f env/environment.yml` or, inside any Python of the
same minor version, `pip install -r env/epi-rna.requirements.lock.txt`.

Re-capture after any package change: `python -m pip freeze > env/epi-rna.requirements.lock.txt`.
(A hash-level cross-platform lock via `pixi` would be the stronger option if the
analysis ever needs to run on a machine other than this Mac; it was not adopted
because it adds a tool the author does not use and the analysis is single-host.)
