# scalableCounterfactual

`scalableCounterfactual` is a research implementation of counterfactual
distribution and quantile decompositions designed for large repeated
cross-sections.

## Start here

- Public repository: <https://github.com/K1mYongrok/scalableCounterfactual>
- Self-contained example: `inst/examples/quick_start.R`
- Architecture: `inst/doc/ARCHITECTURE.md`
- Current release validation: `inst/doc/RELEASE_VALIDATION_1.1.0.md`
- Prior release validation: `inst/doc/RELEASE_VALIDATION_1.0.0.md` (historical
  evidence for version 1.0; it is not a validation claim for 1.1)

The R package is application-neutral. Project-specific wage-worker cleaning,
urban-rural definitions, and empirical specifications remain outside this
folder.

```text
R/          estimation, decomposition, bootstrap, and backend source
man/        function documentation
tests/      API, numerical, GPU, and Stata-parity tests
inst/doc/   architecture, validation, and technical notes
inst/examples/  runnable examples installed with the package
inst/python/    optional CUDA backend
```

Install the public source with:

```r
install.packages("remotes")
remotes::install_github("K1mYongrok/scalableCounterfactual")
```

The package separates three layers:

1. conditional-distribution estimation;
2. counterfactual marginalization and decomposition;
3. execution, bootstrap scheduling, checkpointing, and benchmarking.

For large designs, marginalization is a separate scalability layer. The
default `marginal_method = "auto"` materializes the conditional-draw matrix
only when its estimated size is below `marginal_matrix_max_mb`. Larger runs use
a chunked weighted-histogram selection algorithm. It scans bounded row chunks,
isolates only the bins containing the requested weighted order statistics, and
performs exact weighted selection within those bins. The point estimates retain
the same stable weighted type-7 definition as the matrix implementation.
The resolved method, pass count, histogram size, candidate count, and avoided
matrix size are recorded for each fitted/counterfactual distribution.

## Supported conditional-distribution models

- `qr`: linear conditional quantile regression;
- `cqr`: multi-step censored linear conditional quantile regression;
- `loc`: linear location-shift model;
- `locsca`: linear location-scale model;
- `cox`: weighted Cox duration/transformation regression with optional
  right-censoring status;
- `logit`: logit distribution regression;
- `probit`: probit distribution regression;
- `cloglog`: complementary log-log distribution regression;
- `lpm`: linear-probability distribution regression.

The QR model supports single-quantile, preprocessing, and quantile-process
backends:

- `br`, `fn`: exact single-quantile fits;
- `pfn`: exact single-quantile preprocessing fit;
- `qfnb`, `pfnb`: exact multi-quantile fits;
- `proqreg`, `profn`: exact quantile-process preprocessing fits using
  `quantreg::rq.fit.ppro()` with `br` or `fn` as the reduced-problem solver;
- `onestep`: source-traceable native R transcription of the point algorithm in
  Stata `qrprocess` 1.1.3. It reproduces the Mata traversal,
  weighted scale, bandwidth, Jacobian, score, and fallback sequence and always
  disables design preconditioning;
- `auto`: a documented size- and grid-aware exact-solver selector (`br`/`pfn`
  for short grids and `qfnb`/`pfnb` for quantile processes). Approximate
  `onestep` estimation is never selected implicitly.

Inspect the machine-readable registry with `qr_solver_registry()`. The
`solver_requested`, resolved group-specific solver, exact/approximate status,
implementation, preprocessing size, and one-step fallbacks are retained in the
result metadata. If `rq.fit.ppro()` reports a singular reduced problem or
numerical fixups, `proqreg`/`profn` transparently refit the process with their
corresponding full exact solver and record the reason and fallback quantiles.
No failed solver is silently replaced without metadata.

`model` and `solver` are independent options. A QR solver is used by `qr` and
by each selected-sample refit in `cqr`. CQR accepts `br`, `fn`, `pfn`, `qfnb`,
`pfnb`, or `auto`; process-only approximations are deliberately excluded.
`trimming` changes only the QR conditional grid. CQR, `loc`, and `locsca`
always use the complete midpoint grid `(j - 0.5) / nreg` because dropping its
tails would change those model definitions.

Non-QR computation is also backend-neutral. `loc`, `locsca`, and `lpm` accept
`linear_backend = "qr"`, `"chol"`, or `"fastglm"`; the first two reuse one
weighted factorization across responses. `logit`, `probit`, and `cloglog` accept
`dr_backend = "glm"`, `"fastglm"`, `"speedglm"`, or the optional `"cuda"`
process backend. `fastglm` is the optional RcppEigen-backed path. `auto` uses
`fastglm` for distribution regression when installed and otherwise base-R
GLM; the conservative linear default is QR.
Inspect these mappings with `conditional_backend_registry()`.

Sequential distribution regression can warm-start neighboring thresholds.
Set `dr_workers > 1` to fit thresholds in parallel; warm starts are then
disabled. Threshold parallelism cannot be nested inside group-level or
bootstrap-replication parallelism. An invertible weighted-design
preconditioner is enabled by default for logit/probit/cloglog and coefficients are
returned in their original units.

For logit, probit, and cloglog, built-in CPU and CUDA fits supplement backend
convergence flags with a scale-aware signed-margin diagnostic for complete and
quasi-complete separation. A separated fit is reported as a boundary fit; an
isolated extreme fitted probability is not sufficient by itself. CQR retries
an invalid first-stage backend with base-R GLM and stops if that fit also
remains nonconverged or separated.

Separately fitted distribution-regression thresholds can yield a nonmonotone
marginal CDF in finite samples. The 1.x-compatible default
`dr_noncrossing = "cummax"` preserves the running-maximum rule used in version
1.0. Set it to `"rearrange"` for the increasing rearrangement used by
`Counterfactual` 1.2, `"isotonic"` for an equal-weight isotonic projection, or
`"none"` to disable monotonicity correction. Every policy still enforces the
probability bounds `[0, 1]` before inversion. Marginalization diagnostics record the raw crossing
count and largest violation, corrected crossing count, number of out-of-bounds
values (`dr_out_of_bounds_values`), largest bound correction
(`dr_max_bound_adjustment`), and largest total correction.

Point-estimate QR and CQR runs also diagnose adjacent conditional-quantile crossing on
the reference design, the reference structure evaluated on comparison
covariates, and the comparison design. Set `crossing_diagnostics = FALSE` to
skip the additional passes. Diagnostics report row and adjacent-pair shares,
survey-weighted row shares, and crossing magnitudes. The default
`quantile_noncrossing = "none"` never changes the fitted process. The explicit
`"rearrange"` option sorts conditional quantiles within each covariate row and
records both raw and corrected diagnostics. Rearrangement preserves each row's
conditional draw distribution, so counterfactual marginal quantiles are
unchanged by construction.

CQR follows the multi-step selection/refitting estimator: a weighted logit
model first estimates the uncensoring probability, and exact weighted QR fits
are iteratively restricted to observations whose predicted latent quantiles
clear the censoring point. `censoring` may be a column, vector, or scalar;
`cqr_right = TRUE` handles right censoring by sign reversal. The first-stage
logit uses the selected `dr_backend`; a nonconverged or boundary fit is retried
with base-R `glm`. Estimation stops if that fallback remains nonconverged or on
the probability boundary. Group-specific metadata retain the initial and final
backend, fallback indicator, convergence and boundary status, iteration count,
and accumulated warnings.

Cox fits use Breslow ties and optional `event` status. Marginal Cox quantiles
are inverted with a bounded-memory binary search over the estimated baseline
hazard rather than materializing an observation-by-event-time matrix. Observed
reference quantiles use the same survey-weighted Kaplan--Meier definition. If a
requested quantile lies above the identified Kaplan--Meier CDF, the safe default
`cox_boundary = "na"` returns `NA`; `"error"` stops the run and `"cap"` retains
the former last-event-time behavior. Unidentified quantiles are excluded from
bootstrap and functional inference rather than treated as estimates.
Every model formula must include an intercept. All non-Cox models accept an
intercept-only formula, subject to the CQR overlap and identified-quantile
conditions; this estimates an unconditional two-group comparison with a zero
composition component. Cox regression requires at least one non-intercept
covariate.

These backends preserve the corresponding unpenalized model objective.
Regularized `glmnet` fits are deliberately not labelled as solver backends,
because a penalty changes the estimator rather than only its computation.

Before estimation, both groups are forced to use one common, full-rank design
matrix. A scale-aware maximum common independent-column calculation retains as
many formula columns as possible while requiring full rank in each group.
Columns that cannot be jointly identified are removed from both groups and
recorded as `dropped_design_columns` in the run metadata. Solver failures are
retained as benchmark results rather than silently replaced by another solver.

By default, QR fits also apply an invertible design preconditioner internally
and transform coefficients back to their original units. This changes numerical
conditioning, not the model or estimand. Set `qr_precondition = FALSE` in
`cf_control()` for a strict unscaled legacy comparison.

All empirical and simulated weighted quantiles are normalized within the
distribution being evaluated. This makes results invariant to a common
rescaling of sampling weights. It is an intentional numerical correction to
the default `Hmisc::wtd.quantile` behavior used by `Counterfactual` 1.2; strict
source-code output parity is therefore not expected in small samples even when
the conditional BR coefficients are identical. The default rank calculation
agrees with `Hmisc::wtd.quantile(normwt = TRUE)` when cumulative weighted ranks
are numerically distinct. It evaluates those ranks directly so that weights
far below machine precision do not trigger `approx()` tie collapsing; in such
extreme dynamic-range cases, 1.1 results can differ from the numerically
unstable 1.0/Hmisc output while preserving the intended weighted type-7 rule.
For empirical and counterfactual bootstrap schemes, repeated sampled rows are
stored once with integer frequency weights. The type-7 effective sample size
and the location/CQR selection quantiles and sample-size gates retain the
original resample size, so this compressed representation is numerically
equivalent to explicitly duplicating the sampled rows. Multiplier bootstrap
weights keep the ordinary fixed-row normalization.

## Install

From the cloned repository root:

```powershell
R CMD INSTALL .
```

The optional CUDA dependencies are recorded in
`inst/python/requirements-cuda.txt`. Install them in an isolated project-local
virtual environment rather than into a shared user site or `--target`
directory:

```powershell
py -3.12 -m venv tmp/cuda-venv
$env:PYTHONNOUSERSITE = "1"
tmp/cuda-venv/Scripts/python.exe -m pip install `
  -r inst/python/requirements-cuda.txt
```

The exact Windows/Python 3.12 package set used for validation is retained in
`inst/python/requirements-cuda-windows-py312.lock`. Pass the virtual
environment's Python executable as `gpu_python` and its `Lib/site-packages`
directory as `gpu_python_path`. Keep `PYTHONNOUSERSITE=1` set before R
initializes Python. CUDA is optional; the CPU package does not import Python or
CuPy.

## Executable quick start

This example is self-contained and runs after installation:

```r
library(scalableCounterfactual)

set.seed(42)
n <- 600
example_data <- data.frame(
  outcome = rnorm(n),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.4),
  group = rep(0:1, each = n / 2),
  sampling_weight = runif(n, 0.5, 2)
)

fit <- counterfactual_decompose(
  outcome ~ x1 + x2,
  data = example_data,
  group = "group",
  weights = "sampling_weight",
  model = "qr",
  solver = "pfnb",
  control = cf_control(
    nreg = 9,
    reported_quantiles = c(0.25, 0.5, 0.75),
    bootstrap_progress = FALSE
  ),
  bootstrap_reps = 3,
  seed = 42
)

as.data.frame(fit)
fit$functional_tests
```

This quick start is deliberately application-neutral. Empirical sample rules,
variable construction, and the paper specification belong in a separate
analysis repository or script and are not defaults of the package.

## Generic R API

```r
fit <- counterfactual_decompose(
  log_outcome ~ age + I(age^2) + education + occupation,
  data = analysis_data,
  group = "group_indicator",
  weights = "sampling_weight",
  model = "qr",
  solver = "profn",
  control = cf_control(
    bootstrap_scheme = "counterfactual",
    qr_bootstrap_engine = "xy_preprocess"
  ),
  bootstrap_reps = 100,
  point_workers = 2,
  bootstrap_workers = 5,
  checkpoint_dir = "output/checkpoints"
)
```

Group 0 is the reference population and group 1 is the comparison population.
All reported effects are group 1 minus group 0.

## Generic command-line interface

The included CLI accepts an arbitrary CSV or RDS data frame, model formula,
binary group column, and optional sampling-weight column. Dataset preparation
and variable construction remain the caller's responsibility.

```powershell
Rscript inst/scripts/cfdecomp.R `
  --data data/analysis_sample.csv `
  --formula "log_outcome ~ age + I(age^2) + education + occupation" `
  --group group_indicator `
  --weights sampling_weight `
  --model qr `
  --solver profn `
  --reps 100 `
  --point-workers 2 `
  --bootstrap-workers 5 `
  --bootstrap-scheme counterfactual `
  --qr-precondition true `
  --qr-bootstrap-engine xy_preprocess `
  --output output/counterfactual
```

For a distribution-regression run, for example:

```powershell
Rscript inst/scripts/cfdecomp.R `
  --data data/analysis_sample.csv `
  --formula "log_outcome ~ age + I(age^2) + education + occupation" `
  --group group_indicator `
  --weights sampling_weight `
  --model logit `
  --dr-backend fastglm `
  --dr-workers 1 `
  --dr-warm-start true `
  --dr-precondition true `
  --dr-noncrossing rearrange `
  --output output/counterfactual_logit
```

The resampling law and fitting engine are separate controls. Set
`bootstrap_scheme` to `"counterfactual"` for the exponential-probability
resampling used by `Counterfactual` 1.2, `"empirical"` for ordinary within-group
resampling, or `"multiplier"` for exponential multiplier weights. Then choose
the computational engine with `qr_bootstrap_engine`: `"standard"` refits the
requested solver, `"xy_preprocess"` uses an exact residual-sign preprocessing
refit for `profn` and `proqreg`, and `"onestep"` reuses the point-estimate
inverse Jacobians. `"auto"` selects a compatible specialized engine. The
`xy_preprocess` engine does not implement multiplier draws. Point-estimation
and bootstrap worker counts remain separate options.
The requested point count is retained in metadata, while the effective count
is capped at two because there are exactly two group-specific conditional
fits. With multiple bootstrap workers, the effective within-draw point count
is one to prevent nested process parallelism.

Example of the Stata-compatible path:

```powershell
Rscript inst/scripts/cfdecomp.R `
  --data data/analysis_sample.csv `
  --formula "log_outcome ~ age + I(age^2) + education + occupation" `
  --group group_indicator `
  --weights sampling_weight `
  --model qr `
  --solver onestep `
  --onestep-first-solver auto `
  --onestep-bandwidth hall_sheather `
  --qr-bootstrap-engine onestep `
  --bootstrap-scheme multiplier `
  --reps 100 `
  --bootstrap-workers 5 `
  --output output/counterfactual_onestep
```

`onestep` is algorithm-compatible, not a claim of bitwise identity with
Mata: its initial exact QR fit and matrix inversion are supplied by R and
`quantreg`. Set the `STATA_EXE` environment variable before running the package
tests to enable the optional cross-language coefficient parity test.

## Solver benchmark

The benchmark command uses one frozen design matrix, quantile grid, weight
vector, and worker count for every solver. It records elapsed time, peak R heap
per process, convergence flags, and numerical differences from a selected
reference solver.

```powershell
Rscript inst/scripts/cfbenchmark.R `
  --data data/analysis_sample.csv `
  --formula "log_outcome ~ age + I(age^2) + education + occupation" `
  --group group_indicator `
  --weights sampling_weight `
  --sample-n 5000 `
  --solvers br,fn,pfn,qfnb,pfnb,proqreg,profn,onestep,auto `
  --reference br `
  --repetitions 5 `
  --warmup 1 `
  --output output/benchmark/qr_solvers.csv
```

The benchmark command writes repeat-level results to the requested path and a
second `_summary.csv` file containing median times and relative speed. The two
files are staged and committed as one managed output set, so a failed write
does not leave a new raw file paired with an old summary.

Use `benchmark_conditional_backends()` for frozen-design comparisons of
linear and distribution-regression backends. It records resolved backends,
errors, warnings, time, R heap, and maximum effect differences from a chosen
reference backend. When groups have different numbers of distinct outcome
thresholds, group-specific grid sizes and effective threshold-worker counts
are reported; the single combined field is `NA` rather than implying that the
two groups used the same process.

The ordinary benchmark's peak-memory field is the peak R heap reported by
`gc()`. `benchmark_qr_scaling()` additionally samples fresh-process RSS from
the operating system.

Bootstrap data are copied to each PSOCK worker once at initialization. Only a
replication identifier is sent for subsequent tasks. Checkpoint identities use
a complete hash of the prepared outcomes, weights, and design matrices and do
not include the requested replication count, so a run can be extended without
discarding completed replications. The signature also contains the package
version, sampling scheme, fitting engine, solver controls, and random seed, so
checkpoints from a changed implementation are not silently reused.

Bootstrap replications are fault-tolerant. A failed draw is retried with a
deterministic replacement seed up to `bootstrap_max_retries`; every failed
attempt is written to `bootstrap_failures.csv`, and successful checkpoints
retain their attempt number and actual seed. Parallel work is submitted in
worker-sized batches so `bootstrap_progress = TRUE` can report completed draws,
elapsed time, and ETA. If the retry budget is exhausted, completed checkpoints
remain reusable and the error points to the failure log.
Retries protect long runs from isolated numerical failures; a material failure
rate is still a diagnostic problem and should be reported rather than treated
as harmless missing computation.

`write_cf_outputs()` additionally writes
`marginalization_diagnostics.csv`, `quantile_crossing_diagnostics.csv`, and,
when any retry occurred, `bootstrap_failures.csv`. Runs with at least two
bootstrap replications also write `functional_effect_tests.csv`. Use
`functional_effect_tests(fit, constants = c(0.05, 0.10))` to add specific
constant-effect nulls. The reported studentized KS and Cramer--von Mises
p-values jointly evaluate the selected effect curve rather than testing one
quantile at a time.

For publication-quality scaling measurements, use the isolated runner:

```powershell
Rscript inst/scripts/cfscaling.R `
  --data data/analysis_sample.csv `
  --formula "log_outcome ~ age + I(age^2) + education + occupation" `
  --group group_indicator `
  --weights sampling_weight `
  --sample-sizes 5000,20000,100000 `
  --solvers fn,qfnb,pfnb,proqreg,profn `
  --reference pfnb `
  --repetitions 3 `
  --warmup 1 `
  --checkpoint output/benchmark/qr_scaling_checkpoint.rds `
  --output output/benchmark/qr_scaling.csv
```

Each solver/repetition is launched in a fresh R process. The output therefore
adds sampled process RSS to the ordinary in-session R-heap measurement and
avoids memory retained by a previously run solver. With `point_workers > 1`,
RSS refers to the parent worker process only; use one point worker for a clean
solver-memory comparison. The checkpoint is atomically updated after every
recorded solver run, so an interrupted command resumes completed conditions in
the same repetition order. Its identity includes the R, package, BLAS/LAPACK,
extension, and applicable GPU runtime; stored rows and fit shapes are validated
before reuse.

Strict source auditing is separate from the recommended default. Retaining the
1.x-compatible `legacy_qr_shift = TRUE` and setting
`legacy_weighted_quantile = TRUE`,
`qr_precondition = FALSE`, and `marginal_method = "matrix"` reproduces the
scale-sensitive QR path in `Counterfactual` 1.2. The default
`legacy_weighted_quantile = FALSE` is weight-scale invariant and should remain
the choice for new empirical work. Setting `legacy_qr_shift = FALSE` removes
the legacy scalar level shift; it does not change decomposition differences.

## Custom solvers and model backends

The computation layer is extensible without editing package internals. A QR
solver receives the original design, outcome, normalized case weights,
quantile vector, and control list; it returns a coefficient matrix in original
design units. Registered functions are serialized to PSOCK workers and their
fingerprints are included in bootstrap and scaling-checkpoint identities.

```r
register_qr_solver(
  "my_solver",
  fit = function(X, y, weights, taus, control) {
    coef <- vapply(taus, function(tau) {
      quantreg::rq.fit.fnb(X * weights, y * weights, tau)$coefficients
    }, numeric(ncol(X)))
    list(coefficients = coef, flag = rep(0L, length(taus)))
  },
  exact = TRUE,
  version = "1.0"
)

fit <- counterfactual_decompose(
  outcome ~ x1 + x2,
  data = analysis_data, group = "group", weights = "sampling_weight",
  solver = "my_solver"
)
unregister_qr_solver("my_solver")
```

Direct lexical dependencies, including values referenced from the global
environment, are captured when an extension is registered. Supply
`dependencies = list(name = value)` for objects reached through dynamic lookup
and change `version` whenever external implementation details change. Use
`register_conditional_backend()` for custom linear or binary-response engines.
The complete function contracts are documented in the corresponding help
pages.

## Summary, plotting, and simulation validation

```r
summary(fit, effects = c("structure", "composition"))
plot(fit, interval = "uniform")

validation <- simulate_counterfactual_validation(
  replications = 100,
  n_per_group = 2000,
  solver = "fn",
  workers = 4,
  checkpoint_dir = "output/simulation_validation/checkpoints",
  seed = 20260809
)
write_simulation_validation(validation, "output/simulation_validation")
```

The validation suite uses Gaussian designs whose structure, composition, and
total quantile effects are available analytically. Optional nested bootstrap
replications add pointwise and simultaneous-band coverage diagnostics.
Replication-level checkpoints and deterministic retries preserve completed
work across interrupted runs. `task_timeout_seconds` is best-effort because
compiled solvers may not check R interrupts. The same workflow is available
through `inst/scripts/cfsimulate.R`.

Simulation metadata retain both requested and effective worker counts; the
effective count cannot exceed the number of tasks. Decomposition and
simulation output writers stage and validate their required managed files
before replacing an existing output set. Unrelated files and checkpoint
subdirectories in the destination are left untouched. The command-line
adapters reject input, output, summary, and checkpoint path collisions before
estimation.

## API and output stability

The fitted object and `run_metadata.csv` record an
`output_schema_version`. Schema version 1.1 retains the group direction as
group 1 minus group 0 and the required 1.0 decomposition columns, and adds DR
CDF-bound/noncrossing diagnostics and additive run metadata for CQR and GPU
execution. The 1.0 contracts remain archived as the prior baseline. Renaming a
required output or changing an existing column's meaning requires a major
release. See `inst/doc/API_STABILITY.md` for the complete policy.

Approximate `onestep` and experimental `cuda_admm` computation remain clearly
labelled and are not covered by a claim of exact-solver equivalence. The 1.x
stability promise concerns documented estimands and interfaces rather than
bitwise identity across hardware and numerical libraries.

## Attribution and scope

The counterfactual estimands, QR marginalization, and weighted-bootstrap scheme
follow the CRAN `Counterfactual` package and its underlying methodology. The
quantile-process preprocessing and one-step interfaces are informed by the CFM
algorithms and Blaise Melly's Stata `qrprocess` implementation. See
`inst/provenance/METHODS.md` for the component-level paper and source map, and
`inst/NOTICE` for licensing information. The map explicitly distinguishes
methods implemented from papers, delegated R-package fits, and translated
source code.

This is an independent modular implementation, not a fork or replacement of
either package. `proqreg` and `profn` call the public `quantreg::rq.fit.ppro()`
implementation. `onestep` is approximate and is labelled as such throughout
the output. Censored QR and Cox conditional models are separate extensions;
they do not change the definitions of the reported structure, composition,
and total effects. Functional inference remains limited to the implemented
KS/CMS-style interfaces documented by the package.

## Optional CUDA backend

The CPU path remains the default and has no Python dependency. A CUDA device
can optionally accelerate dense prediction and counterfactual marginalization:

```r
cuda_control <- cf_control(
  gpu_backend = "cuda",
  gpu_precision = "float64",
  gpu_python = "C:/path/to/python.exe",
  gpu_python_path = "C:/path/to/project-local-python-packages"
)
```

For logit, probit, and complementary-log-log distribution regression,
`dr_backend = "cuda"` fits thresholds in GPU batches with the unpenalized
binomial objective. Boundary thresholds should be compared through fitted
probabilities and objective values because separated binary regressions need
not have unique finite coefficients.

`solver = "cuda_admm"` is an explicitly experimental QR solver. It is not an
exact replacement for PFNB and stops on nonconvergence by default. Keep
`gpu_qr_allow_nonconvergence = FALSE` for substantive estimation and use an
exact CPU solver with `gpu_backend = "cuda"` when only marginalization should
be accelerated.

For QR and CQR, CUDA keeps conditional draws, sorting, cumulative survey
weights, and marginal-quantile selection on the device when
`legacy_weighted_quantile = FALSE`. The returned quantiles use the same rank
and interpolation rule as the CPU implementation. Set
`marginal_method = "chunked"` explicitly when the GPU sorting workspace would
not fit in device memory.

CUDA bootstrap runs currently use one R process so that the device is not
oversubscribed. Every replication uses the selected GPU kernels, while the
existing checkpoint, retry, and RNG logic remains unchanged. CUDA runtime paths
and imported Python modules are cached by normalized Python/runtime paths plus
the CUDA module's path and content hash; repeated calls do not repeatedly
extend `PATH` or `PYTHONPATH`.

Backend requests are validated against the model before fitting. For example,
CUDA DR fitting is available only for logit, probit, and cloglog, and Cox does
not advertise CUDA execution. Metadata records `fit_device`,
`prediction_device`, `marginalization_device`, `requested_backend`, and the
resolved group-specific fit backends separately. Consequently an exact CPU
PFNB fit may use two group workers and then pass its fitted process to a
single-device CUDA marginalization step without being mislabeled as a CUDA fit.
