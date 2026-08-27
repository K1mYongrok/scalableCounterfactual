# Architecture

The package deliberately separates the statistical estimator from execution.

## Naming and comment conventions

Public names describe statistical roles rather than a particular empirical
application. `outcome`, `group`, and `weights` denote the response, binary
comparison indicator, and optional sampling weights. Group 0 is the reference
group and group 1 is the comparison group throughout. A `model` determines the
conditional statistical model; a QR `solver` changes the numerical algorithm
for a fixed QR objective; a `backend` changes an implementation used within a
model. These terms are not interchangeable.

Plural names are used for vectors or collections (`weights`, `taus`,
`workers`, `quantiles`); singular names identify one value or one column.
`nreg` is retained as the established public name for the conditional process
grid size. Public arguments are not abbreviated further. Comments explain a
statistical or numerical reason and avoid project-specific variable names,
informal status notes, and claims about quality that are not tied to a test or
reference.

## 1. Data and common design

`prepare_cf_data()` evaluates the formula once, applies one complete-case mask,
normalizes survey weights, splits reference group 0 and comparison group 1,
and enforces a common full-rank design. Removed columns are recorded in run
metadata.

## 2. Conditional-distribution estimation

`fit_conditional_model()` dispatches to the following model families:

- conditional quantile regression (`qr`) and multi-step censored quantile
  regression (`cqr`);
- location or location-scale regression (`loc`, `locsca`);
- Cox duration/transformation regression (`cox`);
- distribution regression (`logit`, `probit`, `cloglog`, `lpm`).

The non-QR layer separates the statistical model from its computational
backend. Linear models can reuse QR or Cholesky WLS factorizations or use the
optional RcppEigen-backed `fastglm` implementation. Binary distribution
regression can use base-R GLM, `fastglm`, or `speedglm`. Logit, probit, and
cloglog fits
use an invertible weighted-design transformation and return coefficients in
the original design units.

Sequential DR fits may warm-start neighboring thresholds. Alternatively,
thresholds can be distributed across a PSOCK cluster after initializing the
design once per worker. Threshold, group, and bootstrap parallelism are kept
as mutually exclusive layers to avoid nested oversubscription.

The QR backend is solver-neutral. `fit_weighted_qr()` applies positive case
weights before calling the selected public `quantreg` fitting routine. This is
important for the multi-quantile QFNB and PFNB routines, whose low-level APIs do
not accept a separate weights argument.

The QR design is preconditioned by default with an invertible transformation to
reduce numerical problems caused by mixed scales and correlated columns such as
age and age squared. Coefficients are transformed back to their original units
before marginalization. The option is explicit in `cf_control()` and recorded
in metadata.

CQR adds a weighted logit uncensoring model followed by repeated selected-
sample QR fits. Censoring points are carried through group splitting,
subsampling, bootstrap resampling, and checkpoint fingerprints. Cox removes
the formula intercept, fits weighted partial likelihood with Breslow ties, and
constructs a matching weighted Breslow baseline hazard. Cox warnings,
iteration counts, iteration-limit failures, and possible infinite coefficients
are retained in standardized convergence diagnostics.

The fast path uses a Cholesky preconditioner. An ill-conditioned Gram matrix
triggers a pivoted-QR fallback. The method, reciprocal condition estimate, and
fallback reason are retained in each fit and benchmark record.

## 3. Marginalization and decomposition

Each group-specific conditional model is integrated over the observed
covariate distribution. With group 0 as the reference and group 1 as the
comparison, the package constructs

- fitted reference distribution: `F_0|0` integrated over `X_0`;
- counterfactual distribution: group-0 conditional structure integrated over
  `X_1`;
- fitted comparison distribution: `F_1|1` integrated over `X_1`.

At every reported quantile, total equals structure plus composition. The
identity residual is checked numerically.

Draw-based models have two marginalization implementations. `matrix` preserves
the original direct materialization. `chunked` first obtains the prediction
range, builds a survey-weighted histogram in bounded row chunks, and retains
only values in bins containing the low/high weighted order statistics required
by the Hmisc-compatible interpolation rule. The final selection is exact within
those bins. `auto` resolves between the implementations from the estimated
matrix size, and every resolution is stored in diagnostics.

QR and CQR point estimates optionally make one chunked diagnostic pass over each of
the three relevant design/structure combinations. Adjacent conditional
quantiles differing by less than minus the configured tolerance are summarized
without altering estimates. If row-wise rearrangement is requested, diagnostics
are emitted for both the raw and corrected predictions. Sorting within rows
preserves the draw multiset used for marginalization. Cox marginalization uses
binary search over its baseline-hazard grid and never constructs the full
observation-by-event-time CDF matrix.

Observed Cox reference distributions use survey-weighted Kaplan--Meier
quantiles. Requested probabilities above the identified Kaplan--Meier CDF are
handled by the explicit `cox_boundary` policy: `na` (default), `error`, or the
legacy-compatible `cap`. Unidentified point quantiles do not enter bootstrap or
functional inference.

## 4. Execution and inference

Point-estimation workers and bootstrap workers are separate settings. Point
estimation can fit the two groups concurrently; bootstrap workers distribute
independent replications. Each bootstrap replication is written atomically to
a checkpoint file and can be reused after interruption when the run signature
matches.

When `bootstrap_workers = 1`, the standard refit bootstrap may reuse
`point_workers = 2` to fit the two CPU group models concurrently within a draw.
Their fitted processes return to the parent process before CUDA prediction and
marginalization, so only one process accesses the GPU. With multiple bootstrap
workers, within-draw point fitting resolves to one worker to avoid nested CPU
parallelism.

Device planning is also stage-specific. The package validates model/backend
combinations before reducing worker counts and records the requested backend,
resolved fit backend, fit device, prediction device, and marginalization device
separately. Thus CPU group fits may run concurrently before serial CUDA
marginalization, whereas CUDA conditional fitting and CUDA bootstrap reject
process-level device oversubscription.

The bootstrap sampling law is independent of the fitting engine. The
`counterfactual` scheme uses exponential-probability within-group resampling,
`empirical` uses ordinary within-group resampling, and `multiplier` uses
exponential weights. A QR draw can then be fitted with a standard refit, the
exact XY preprocessing path where supported, or the saved-Jacobian one-step
path. Thus timing comparisons do not inadvertently compare different
resampling laws.

On Windows PSOCK clusters, the prepared matrices are initialized once per
bootstrap worker. Replication tasks contain only integer replication IDs, so the
full design is not serialized for every replication. Checkpoint signatures hash
the complete prepared data, package version, sampling law, fitting engine, and
estimation configuration but exclude the target replication count, allowing an
existing run to be extended without accepting stale checkpoints.

Each requested replication has a stable checkpoint identity. Failed numerical
draws receive deterministic replacement seeds and are retried within the same
identity. Attempt errors are persisted before a terminal failure, while a
successful replacement checkpoint records the attempt, seed, and prior errors.
Parallel replications are dispatched in worker-sized batches to permit progress
and ETA updates without repeatedly serializing the prepared matrices.

## 5. Output and audit trail

Every run writes estimates, distribution diagnostics, phase-level timing and R
heap measurements, metadata, and a serialized fit object. Metadata include the
formula, model, solver, design columns, removed columns, sample sizes, worker
settings, seeds, and package versions.

The point-estimation seed is applied to sequential fits and parallel RNG
streams. Bootstrap replications use deterministic replication-specific seeds.

Weighted quantiles normalize weights within every evaluated distribution. This
ensures invariance to weight scale and intentionally differs from the
small-sample behavior of `Counterfactual` 1.2 with the default Hmisc settings.
An explicit audit-only switch restores the original scale-sensitive call and
forces matrix marginalization; the preserved CRAN source validation harness
uses it to distinguish source parity from the recommended estimator default.

Bootstrap effect draws also support joint studentized KS and Cramer--von Mises
tests of zero, constant, quantile-invariant, and one-sided effect curves. These
tests operate on structure, composition, and total effects separately.

The scaling benchmark runs every solver/repetition in a fresh R process and
polls its operating-system resident set size. Timing and numerical comparisons
therefore remain separated from memory retained by earlier solver calls.

## Optional GPU computation layer

`gpu_backend.R` owns optional runtime discovery and R/Python conversion.
`inst/python/scalablecf_cuda.py` contains CuPy kernels. Statistical model
selection remains in `conditional_model.R`; CUDA is only a computation device.

- `marginal_distribution.R` dispatches dense prediction and DR probability
  aggregation to CUDA when requested.
- `model_dr.R` dispatches an entire threshold process to batched CUDA IRLS and
  uses CPU GLM only for thresholds that do not converge.
- `solver_qr.R` exposes `cuda_admm` as an explicit experimental solver; it is
  not part of automatic solver selection.
- `decomposition_engine.R` prevents multiple R processes from oversubscribing
  one GPU. Existing bootstrap checkpoints and deterministic seeds are retained.
