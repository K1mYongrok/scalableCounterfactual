# API and output stability policy

## Scope

The 1.x public API consists of the functions exported in `NAMESPACE`, their
documented arguments, the registered methods for the `cfdecomp`,
`summary.cfdecomp`, and `cf_simulation_validation` S3 classes, and the files
written by `write_cf_outputs()`. Additive arguments and columns may be
introduced in a minor release. Removing or renaming an exported function,
argument, effect, registered method, required output file, or required output
column requires a major release.

The output schema version is stored as `output_schema_version` in
`run_metadata.csv` and the fitted object's metadata. Version `1.1` requires:

- `decomposition.csv` with `model`, `solver`, `quantile`, `effect`,
  `estimate`, `identified`, `std_error`, pointwise and uniform confidence
  limits, and bootstrap replication counts;
- `distribution_diagnostics.csv`, `point_resources.csv`,
  `marginalization_diagnostics.csv`, `run_metadata.csv`, and `fit.rds`;
- conditional files such as bootstrap resources, failures, functional tests,
  and crossing diagnostics only when the corresponding calculation exists.

Adding or reordering columns increments the schema minor version and requires
an explicit contract update. A change in the meaning, sign, normalization, or
units of an existing column increments the schema major version.

The ordered exported-function signatures and output fields are recorded in
`inst/schema` and enforced by `tests/api_contract.R`. `public_api_1.0.csv`
preserves the initial interface, while `public_api_1.1.csv` appends the
explicit distribution-regression noncrossing policy and the optional
`fit_weighted_qr()` row-frequency argument. `s3_methods_1.1.csv` carries the
initial decomposition methods forward and records the already registered
simulation-validation print method that was absent from the 1.0 inventory.
Output schema 1.1 keeps the decomposition table and file set unchanged while
adding documented columns to `marginalization_diagnostics.csv`. An intentional
interface change requires an explicit edit to the contract, its versioned
schema files, documentation, and test in the same change.

## Statistical defaults

Group effects are always comparison group 1 minus reference group 0. Sampling
weights are normalized within each evaluated distribution, making the default
weighted-quantile estimand invariant to common weight rescaling. The legacy
scale-sensitive rule remains opt-in and is not the 1.x default estimand.
Version 1.1 preserves the version-1.0 `legacy_qr_shift = TRUE` prediction
default and the running-maximum DR CDF correction, now named
`dr_noncrossing = "cummax"`. Rearrangement and isotonic projection are additive
opt-in policies rather than silent changes to existing 1.x results.

Exact QR solvers and approximate solvers remain explicitly distinguished in
metadata. `onestep` and `cuda_admm` are not promoted to exact estimators by the
1.x stability promise. CUDA may change the computation device but not the
documented estimand.

## Numerical compatibility

Patch releases preserve the estimator and output contract, not bitwise
floating-point identity across operating systems, BLAS libraries, CPUs, GPUs,
or solver implementations. Parity is evaluated with documented coefficient,
effect, objective, confidence-interval, and convergence tolerances.

Version 1.1 evaluates the normalized weighted type-7 rank rule directly in
matrix, chunked, and CUDA marginalization. It agrees with the prior Hmisc path
for ordinary weight ranges and avoids `approx()` collapsing numerically tied
cumulative ranks when weights span an extreme dynamic range. This is treated
as a numerical correction to the same rank definition, not a new estimand.

Integer bootstrap frequency weights use the expanded resample size for type-7
rank positions, location and CQR selection quantiles, and sample-size gates.
This makes compressed empirical and counterfactual draws equivalent to
explicit row duplication while leaving point-estimate survey-weight
normalization unchanged.

For CUDA work, the runtime fingerprint records the Python, NumPy, CuPy, CUDA
and device versions plus the canonical module file and its SHA-256. The module
loader also verifies Python's `__file__`. These execution details are metadata,
not part of the statistical estimand, but a fingerprint change invalidates
checkpoint reuse.

Every public model formula currently requires an intercept; non-Cox models may
be intercept-only subject to model identification, while Cox also requires a
non-intercept covariate. Public output writers replace only their declared
managed file set and commit it transactionally. Unrelated destination files
and checkpoint directories are outside that contract and remain untouched.

## Deprecation

An argument scheduled for removal is retained with a warning for at least one
minor 1.x release. Checkpoint identities include the package version and active
backend fingerprints, so checkpoints from an incompatible implementation are
not silently reused.
