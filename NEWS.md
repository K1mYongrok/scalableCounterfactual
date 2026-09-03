# scalableCounterfactual 1.1.0

- Corrected bootstrap inference when successful draws have partial missing
  support, and skipped automatic functional tests when fewer than two usable
  quantiles are available.
- Corrected the legacy weighted-bootstrap evaluation denominator and hardened
  deterministic checkpoint identities against runtime, extension, and
  simulation-scenario changes.
- Hardened censored quantile regression against failed selection fits,
  inconsistent censoring definitions, non-finite designs, and ambiguous
  boundary solutions. The selected first-stage backend now falls back to
  base-R GLM after nonconvergence or a boundary solution, records initial and
  final backend diagnostics, and stops if the fallback is still invalid.
- Made design reduction scale aware, aligned matrix, vector, CUDA, and
  bounded-memory evaluation on a stable weighted type-7 quantile rule, and
  added explicit distribution-regression noncrossing policies and diagnostics.
  The direct type-7 implementation agrees with
  `Hmisc::wtd.quantile(normwt = TRUE)` for ordinary weights when cumulative
  weighted ranks are numerically distinct, while avoiding `approx()` tie
  collapse for extreme weight ranges.
- Preserved the 1.0 numerical default `legacy_qr_shift = TRUE`; setting it to
  `FALSE` removes the legacy scalar prediction shift. The shift cancels from
  decomposition effects but changes fitted-distribution levels and diagnostics.
- Exposed the 1.0 running-maximum DR rule as the explicit, backward-compatible
  default `dr_noncrossing = "cummax"`. New opt-in alternatives are increasing
  `"rearrange"` and equal-weight `"isotonic"`; `"none"` disables monotonicity
  correction but still enforces `[0, 1]` probability bounds.
  Diagnostics now distinguish raw and corrected crossings from the count and
  maximum size of probability-bound corrections.
- Advanced the output schema to version 1.1 for the additive DR CDF,
  CQR-selection, and GPU-runtime diagnostics; archived 1.0 contracts remain
  available as prior-release records. The 1.1 S3 inventory also records the
  existing simulation-validation print method omitted from the initial 1.0
  machine-readable inventory.
- Reduced CPU distribution-regression marginalization memory use and improved
  benchmark, simulation, plotting, worker, and all-failed-run diagnostics.
- Strengthened optional CUDA capability checks, module identity, precision
  handling, command-line controls, provenance, and continuous-integration
  isolation.
- Preserved full midpoint conditional grids for CQR, location, and
  location-scale models while limiting tail trimming to QR, and fixed the
  intercept-only LPM matrix shape.
- Added scale-aware complete and quasi-complete separation diagnostics shared
  by CPU and CUDA binary fits, and preserved the process dimension for a
  single regular CUDA threshold.
- Recorded requested and effective point, threshold, and simulation worker
  counts; point fitting is capped at the two available group models and nested
  process parallelism remains disallowed.
- Made decomposition, simulation, benchmark, and scaling output sets
  transactional; added early CLI path-collision guards and structural
  validation for bootstrap, simulation, and scaling checkpoints.
- Made resumed scaling schedules reproduce uninterrupted repetition order,
  retained missing RSS as `NA`, and added runtime fingerprints to scaling and
  simulation identities.
- Corrected the stable type-7 bootstrap rank normalization so integer-frequency
  compressed empirical and counterfactual resamples match explicit row
  duplication, including fitted residual distributions, CQR selection cutoffs,
  effective sample-size gates, matrix/chunked marginalization, and CUDA
  marginal quantiles. `fit_weighted_qr()` now exposes the same optional
  positive-integer row-frequency representation without changing its 1.0
  ordered argument prefix.

# scalableCounterfactual 1.0.0

- Froze the public API and output schema at version 1.0.
- Replaced placeholder package authorship and maintainer metadata with the
  package author's information and added a package citation.
- Created the first stable release and added the public source repository and
  issue-tracker URLs.

# scalableCounterfactual 0.15.0

- Added a versioned 1.0 output contract and regression tests for required
  exports, public arguments, CSV columns, metadata, and RDS classes.
- Added cross-platform CPU and self-hosted CUDA continuous-integration
  workflows plus a release preflight script.
- Recorded the validated CuPy CUDA environment in install and exact lock
  files, and excluded Python bytecode caches from source tarballs.
- Added an API stability policy and a 1.0.0 release-evidence checklist.

# scalableCounterfactual 0.14.1

- Standard bootstrap replications can now fit the two CPU group models with
  `point_workers = 2` when `bootstrap_workers = 1`. CUDA prediction and
  marginalization remain in the parent process, so a single device is not
  oversubscribed.
- Bootstrap checkpoint signatures and metadata record the resolved number of
  point-fit workers used inside each replication.

# scalableCounterfactual 0.14.0

- Cox and weighted Kaplan-Meier quantiles now use an explicit `cox_boundary`
  policy. The safe default returns unidentified censored upper quantiles as
  `NA`; `error` and legacy-compatible `cap` policies are also available.
- Cox observed-distribution diagnostics use survey-weighted Kaplan-Meier
  quantiles, and unidentified point quantiles are excluded from bootstrap and
  functional inference.
- CUDA runtime paths and imported modules are cached without repeatedly
  extending `PATH` or `PYTHONPATH`.
- Backend requests are validated by model, while metadata separately records
  fit, prediction, and marginalization devices and requested/resolved
  backends. CPU PFNB fitting can use two group workers before single-device
  CUDA marginalization.
- Experimental `cuda_admm` now reports its actual internal preconditioning and
  final-iterate convergence flags. CUDA DR records the backend and fallback
  reason for every threshold.
- Cox convergence warnings, iteration limits, and possible infinite
  coefficients are retained as standardized diagnostics.

# scalableCounterfactual 0.13.0

- QR and CQR CUDA marginalization now keeps conditional draws, sorting,
  cumulative survey weights, rank selection, and interpolation on the GPU.
  Only the reported marginal quantiles are returned to R.
- The CUDA weighted-quantile rule matches `Hmisc::wtd.quantile(...,
  normwt = TRUE)` and retains float64 cumulative weights under either GPU
  prediction precision.

# scalableCounterfactual 0.12.0

- Added an optional CuPy CUDA backend for dense prediction, counterfactual
  marginalization, and batched logit/probit/cloglog distribution regression.
- Added float64/float32 controls, project-local Python runtime discovery,
  device status reporting, GPU column batching, and bootstrap/checkpoint
  identity fields.
- Added the experimental `cuda_admm` QR solver with explicit nonconvergence
  safeguards. It is not selected automatically and does not replace PFNB.
- Added CPU/GPU parity tests and documented boundary-fit, precision, worker,
  and sample-size break-even diagnostics.

# scalableCounterfactual 0.11.0

- Added multi-step left- and right-censored linear quantile regression with
  scalar or observation-specific censoring points and selectable exact QR
  solvers.
- Added weighted Cox duration regression with optional event indicators,
  Breslow baseline hazards, and bounded-memory marginal-CDF inversion.
- Added complementary log-log distribution regression across the existing
  GLM backends, warm starts, threshold parallelism, and bootstrap paths.
- Added opt-in row-wise quantile rearrangement for QR and CQR, with raw and
  corrected crossing diagnostics. Rearrangement is distribution preserving,
  so marginal decomposition estimates remain unchanged.
- Added CQR/Cox bootstrap tests, Cox coefficient validation against
  `survival::coxph`, cloglog smoke tests, CLI options, metadata, and docs.

# scalableCounterfactual 0.10.1

- Made custom extension functions self-contained by capturing direct lexical
  dependencies, added explicit dependency and implementation-version fields,
  and included captured closure state in function fingerprints.
- Restricted checkpoint extension fingerprints to active solvers/backends and
  added the previously omitted legacy weighted-quantile setting to bootstrap
  checkpoint identities.
- Added replication-level simulation checkpoints, deterministic retries,
  best-effort task time limits, parallel progress reporting, overflow-safe
  seeds, and a stable `fn` default.
- Separated curve-level uniform coverage from quantile-level pointwise
  coverage, strengthened custom distribution-backend validation, and removed
  stale optional CSV files when output directories are reused.
- Added regression tests for closure fingerprints, global dependencies on
  PSOCK workers, checkpoint invalidation and resume, large seeds, backend
  diagnostics, and stale-output cleanup.

# scalableCounterfactual 0.10.0

- Added public registration APIs for custom QR solvers and linear or
  distribution-regression computation backends. Registered functions are
  propagated to PSOCK point-estimation, bootstrap, and isolated-benchmark
  workers and are included in checkpoint identities.
- Added `summary()` and base-R `plot()` methods for decomposition objects,
  including pointwise and uniform confidence-band selection.
- Added an analytically identified Gaussian simulation-validation suite for
  location, composition, and combined scale changes, with optional nested
  bootstrap coverage assessment and a command-line runner.
- Added integration tests covering custom extensions in parallel bootstrap,
  S3 output, and serial/parallel simulation execution.

# scalableCounterfactual 0.9.0

- Added bootstrap KS and Cramer--von Mises tests of zero, constant,
  quantile-invariant, nonnegative, and nonpositive structure, composition, and
  total-effect curves.
- Added fresh-process QR scaling benchmarks with sampled operating-system RSS,
  randomized solver order, repeated timing, numerical parity, and a generic
  command-line runner.
- Added an explicit audit-only weighted-quantile mode that reproduces the
  preserved `Counterfactual` 1.2 QR decomposition while retaining the
  scale-invariant correction as the default for new analyses.
- Added a frozen-source parity harness, project-specific scaling adapter,
  executable README example, and expanded release-validation documentation.
- Fixed atomic checkpoint replacement so a successful removal of the previous
  RDS is no longer interpreted as an error; added resumable scaling checkpoints.

# scalableCounterfactual 0.8.0

- Added bounded-memory marginal quantiles. The chunked implementation uses a
  weighted histogram to isolate target bins and performs exact weighted
  selection within those bins without materializing the full observation by
  conditional-draw matrix. `auto` selects it above a configurable matrix-size
  threshold.
- Added point-estimate QR crossing diagnostics on the reference design,
  counterfactual covariate design, and comparison design.
- Added deterministic bootstrap replacement draws, configurable retry limits,
  failure logs, per-replication attempt/seed metadata, batched parallel progress,
  and ETA reporting.
- Added marginalization and crossing CSV outputs, CLI controls, regression
  tests against the legacy matrix calculation, and retry-path tests.

# scalableCounterfactual 0.7.1

- Added strict public-API validation for `cf_control()` objects and rejected
  QR solver arguments supplied to non-QR models instead of silently ignoring
  them.
- Made constant-outcome distribution-regression thresholds invariant to any
  invertible design preconditioner, including the pivoted-QR fallback.
- Synchronized the documented `dr_backend = "auto"` behavior with the actual
  `fastglm`-when-available dispatcher.
- Hardened integer validation and expanded run metadata with effective DR
  execution settings, boundary counts, optional backend versions, and the
  resolved checkpoint run directory.

# scalableCounterfactual 0.7.0

- Added objective-preserving computation backends for location,
  location-scale, LPM distribution-regression, logit, and probit models.
- Added reusable QR and Cholesky WLS factorizations, optional RcppEigen-backed
  `fastglm`, and optional `speedglm` binary-response fitting.
- Added invertible DR design preconditioning, neighboring-threshold warm
  starts, threshold-level parallelism, configurable IRLS iteration controls,
  and explicit convergence and backend metadata.
- Added `conditional_backend_registry()` and
  `benchmark_conditional_backends()` for reproducible non-QR backend parity
  and performance comparisons.
- Kept regularized `glmnet` models outside the backend interface because their
  penalized objective changes the statistical estimator.

# scalableCounterfactual 0.6.0

- Separated the bootstrap sampling law (`counterfactual`, `empirical`, or
  `multiplier`) from the QR fitting engine (`standard`, `xy_preprocess`, or
  `onestep`) so engines can be compared under the same resampling draws.
- Added strict binary-group validation, feasible stratified benchmark
  allocation for imbalanced groups, and explicit distribution-regression
  convergence and separation diagnostics.
- Standardized the default QR solver as `auto` across public entry points.
- Added package-version-aware checkpoint signatures and refreshed the generic
  and project command-line adapters, tests, and reproducible benchmarks.

# scalableCounterfactual 0.5.0

- Removed the urban-rural wage-project adapter from the package API. The core
  package now contains no application-specific variable names, outcomes, group
  definitions, or default data paths.
- Replaced the project CLI with generic CSV/RDS interfaces accepting an
  arbitrary formula, binary group column, and optional weight column.
- Added group and weight definitions to run metadata.
- Moved project-specific adapters and historical wage-analysis artifacts back
  to the surrounding research project.

# scalableCounterfactual 0.4.0

- Promoted the source-traceable `qrprocess` implementation to the sole public
  `onestep` solver and renamed its specialized bootstrap engine to `onestep`.
- Removed the earlier generic approximation from all execution paths after
  replacing it with the source-traceable implementation.
- Consolidated one-step controls, CLI options, metadata, documentation, and
  tests under one stable naming scheme.

# scalableCounterfactual 0.3.0

- Added `onestep_stata`, a native R transcription of the point-estimation core
  in `qrprocess` 1.1.3 with source commit and line-level provenance.
- Added Stata-compatible type-2 weighted quantiles, weighted residual scale,
  Hall-Sheather and Bofinger bandwidths, quantile traversal, and exact-fit
  fallback behavior.
- Added the `stata_onestep` bootstrap engine. It uses exponential multipliers
  or empirical bootstrap counts and reuses point-estimate inverse Jacobians as
  in `rq_boot_1step()`.
- Added automatic specialized-bootstrap selection for `onestep_stata` and new
  solver, bandwidth, fallback, source-version, and source-commit metadata.
- Added an optional cross-language parity harness enabled by `STATA_EXE`.

# scalableCounterfactual 0.2.0

- Added exact quantile-process preprocessing solvers `proqreg` and `profn`.
- Added the approximate CFM-style `onestep` solver with explicit fallback
  diagnostics and exactness metadata.
- Added the documented `auto` solver selector and `qr_solver_registry()`.
- Added the exact QR bootstrap residual-sign preprocessing engine for
  `proqreg` and `profn`.
- Added resolved-solver, exactness, implementation, preprocessing, and
  bootstrap-engine metadata to results and benchmarks.
- Extended the CLI and benchmark runner to expose all solver and bootstrap
  options.

# scalableCounterfactual 0.1.0

- Added modular counterfactual distribution and quantile decomposition.
- Added QR, location, location-scale, logit, probit, and LPM conditional models.
- Added BR, FN, PFN, QFNB, and PFNB QR backends.
- Added weighted, checkpointed bootstrap inference with separate worker control.
- Added frozen-design QR solver benchmarks and numerical-parity diagnostics.
- Preserved the legacy project implementation and its provenance hashes.
- Made weighted quantiles invariant to common weight rescaling.
- Initialized bootstrap data once per PSOCK worker and enabled checkpoint
  reuse when extending the replication count.
- Added complete-data checkpoint fingerprints, solver warnings, convergence
  availability, condition diagnostics, and a pivoted-QR preconditioning fallback.
- Added repeated randomized-order benchmarks with warmup support.
