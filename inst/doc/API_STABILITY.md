# API and output stability policy

## Scope

The 1.x public API consists of the functions exported in `NAMESPACE`, their
documented arguments, the `cfdecomp` S3 class, and the files written by
`write_cf_outputs()`. Additive arguments and columns may be introduced in a
minor release. Removing or renaming an exported function, argument, effect,
required output file, or required output column requires a major release.

The output schema version is stored as `output_schema_version` in
`run_metadata.csv` and the fitted object's metadata. Version `1.0` requires:

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
`inst/schema` and enforced by `tests/api_contract.R`. This is the stable 1.0
contract. An intentional interface change therefore requires an explicit edit to the contract,
its versioned schema files, documentation, and test in the same change.

## Statistical defaults

Group effects are always comparison group 1 minus reference group 0. Sampling
weights are normalized within each evaluated distribution, making the default
weighted-quantile estimand invariant to common weight rescaling. The legacy
scale-sensitive rule remains opt-in and is not the 1.x default estimand.

Exact QR solvers and approximate solvers remain explicitly distinguished in
metadata. `onestep` and `cuda_admm` are not promoted to exact estimators by the
1.x stability promise. CUDA may change the computation device but not the
documented estimand.

## Numerical compatibility

Patch releases preserve the estimator and output contract, not bitwise
floating-point identity across operating systems, BLAS libraries, CPUs, GPUs,
or solver implementations. Parity is evaluated with documented coefficient,
effect, objective, confidence-interval, and convergence tolerances.

## Deprecation

An argument scheduled for removal is retained with a warning for at least one
minor 1.x release. Checkpoint identities include the package version and active
backend fingerprints, so checkpoints from an incompatible implementation are
not silently reused.
