#' Configure counterfactual estimation
#'
#' @param nreg Number of conditional quantile regressions or distribution
#'   regression thresholds.
#' @param trimming QR tail trimming fraction.
#' @param reported_quantiles Quantiles of the marginal distributions to report.
#' @param alpha Inference significance level.
#' @param weighted_bootstrap Deprecated logical compatibility argument. `TRUE`
#'   maps to `bootstrap_scheme = "counterfactual"`; `FALSE` maps to
#'   `"empirical"`.
#' @param bootstrap_scheme Bootstrap sampling law. `counterfactual` uses the
#'   exponential-probability resampling in `Counterfactual` 1.2; `empirical`
#'   uses ordinary within-group resampling; `multiplier` applies independent
#'   exponential multipliers directly to observation weights.
#' @param legacy_qr_shift Reproduce the scalar trimming shift used in
#'   `Counterfactual` 1.2 QR predictions.
#' @param legacy_weighted_quantile Reproduce the scale-sensitive
#'   `Hmisc::wtd.quantile()` call in `Counterfactual` 1.2. This audit option
#'   forces matrix marginalization and is not recommended for new analysis.
#' @param qr_precondition Apply an invertible QR design preconditioner internally
#'   for numerical stability, then return coefficients in their original units.
#' @param onestep_first_solver Exact R solver used for the first
#'   (nearest-median) quantile in the one-step process. `auto` follows
#'   the large-sample choice in `qrprocess`: `br` below 100,000 observations
#'   and `pfn` otherwise.
#' @param onestep_bandwidth Bandwidth rule used by the one-step process.
#'   Hall-Sheather is the `qrprocess` default.
#' @param qr_bootstrap_engine QR bootstrap fitting engine. `standard` refits the
#'   selected solver on every draw; `xy_preprocess` uses
#'   [quantreg::boot.rq.pxy()] with the point estimates; `onestep` reuses
#'   the point-estimate inverse Jacobians like `qrprocess`; `auto` selects the
#'   matching specialized engine for supported solvers.
#' @param robust_se Use the bootstrap interquartile-range standard error.
#' @param bootstrap_max_retries Maximum deterministic replacement draws after
#'   a failed bootstrap replication. A value of two allows three attempts per
#'   requested replication.
#' @param bootstrap_progress Print completed replications and an estimated
#'   remaining time while the bootstrap runs.
#' @param marginal_method Marginal-quantile implementation. `matrix` materializes
#'   every conditional draw, `chunked` uses a bounded-memory weighted histogram
#'   selection, and `auto` selects from the estimated draw-matrix size.
#' @param marginal_chunk_rows Maximum design rows predicted in one chunk.
#' @param marginal_matrix_max_mb Matrix-memory threshold used by `auto`.
#' @param marginal_histogram_bins Initial number of bins used by chunked
#'   weighted selection.
#' @param marginal_candidate_max Maximum candidate draws retained in the final
#'   exact-within-bin selection step; the histogram is refined when needed.
#' @param crossing_diagnostics Evaluate QR quantile crossing on the reference,
#'   counterfactual-covariate, and comparison designs for the point estimate.
#' @param crossing_tolerance Minimum negative adjacent-quantile difference
#'   counted as a crossing.
#' @param quantile_noncrossing Conditional-quantile monotonicity correction.
#'   `none` retains the fitted process. `rearrange` sorts fitted conditional
#'   quantiles within each covariate row. Rearrangement removes crossing while
#'   preserving the conditional draw distribution and hence the reported
#'   counterfactual marginal quantiles.
#' @param cqr_right Whether `cqr` uses right rather than left censoring.
#' @param cqr_nsteps Number of censored-QR selection/refitting steps; at least
#'   three.
#' @param cqr_first_cut Weighted fraction trimmed after the first-stage
#'   uncensoring-probability model.
#' @param cqr_later_cut Weighted positive-margin fraction trimmed during the
#'   third and later CQR steps.
#' @param cox_boundary Policy for requested Cox/Kaplan-Meier quantiles above
#'   the largest identified CDF value. `na` returns unidentified quantiles as
#'   `NA`, `error` stops, and `cap` reproduces the legacy boundary behavior of
#'   returning the final event time.
#' @param linear_backend Weighted least-squares backend for `loc`, `locsca`,
#'   and `lpm`: robust base-R QR, Cholesky, RcppEigen-backed `fastglm`, or
#'   `auto` (currently QR for backward-compatible numerical robustness).
#' @param dr_backend Binary-response backend for logit, probit, and cloglog distribution
#'   regression: base-R `glm`, `fastglm`, `speedglm`, or `auto` (currently
#'   `fastglm` when installed and base-R GLM otherwise).
#' @param dr_workers Workers used across distribution-regression thresholds.
#'   This execution layer cannot be combined with group-level or
#'   replication-level parallelism.
#' @param dr_warm_start Reuse neighboring-threshold coefficients as starting
#'   values during sequential logit/probit/cloglog distribution regression.
#' @param dr_maxit Maximum IRLS iterations at each distribution threshold.
#' @param dr_tolerance IRLS convergence tolerance.
#' @param dr_precondition Apply an invertible weighted-design transformation
#'   before logit/probit/cloglog fitting and return coefficients in original units.
#' @param gpu_backend Prediction and marginalization device. `cuda` uses the
#'   optional CuPy runtime; `cpu` preserves the dependency-free default.
#' @param gpu_precision CUDA arithmetic precision. `float64` is the auditable
#'   default. `float32` is experimental.
#' @param gpu_python Optional Python executable containing CuPy.
#' @param gpu_python_path Optional project-local Python package directory.
#' @param gpu_module_path Optional development override for the bundled module.
#' @param gpu_block_columns Quantile or threshold columns per CUDA block.
#' @param gpu_qr_rho ADMM penalty for the experimental `cuda_admm` QR solver.
#' @param gpu_qr_maxit Maximum ADMM iterations.
#' @param gpu_qr_tolerance ADMM primal and dual convergence tolerance.
#' @param gpu_qr_allow_nonconvergence Allow an experimental `cuda_admm` result
#'   that reached its iteration limit to be returned with diagnostics. The
#'   default stops before decomposition or inference.
#' @return A validated control list.
#' @export
cf_control <- function(
    nreg = 100L,
    trimming = 0.005,
    reported_quantiles = seq(0.1, 0.9, by = 0.1),
    alpha = 0.05,
    weighted_bootstrap = NULL,
    bootstrap_scheme = c("counterfactual", "empirical", "multiplier"),
    legacy_qr_shift = TRUE,
    legacy_weighted_quantile = FALSE,
    qr_precondition = TRUE,
    onestep_first_solver = c("auto", "br", "fn", "pfn"),
    onestep_bandwidth = c("hall_sheather", "bofinger"),
    qr_bootstrap_engine = c(
      "standard", "xy_preprocess", "onestep", "auto"
    ),
    robust_se = FALSE,
    bootstrap_max_retries = 2L,
    bootstrap_progress = interactive(),
    marginal_method = c("auto", "matrix", "chunked"),
    marginal_chunk_rows = 50000L,
    marginal_matrix_max_mb = 512,
    marginal_histogram_bins = 262144L,
    marginal_candidate_max = 2000000L,
    crossing_diagnostics = TRUE,
    crossing_tolerance = 1e-8,
    quantile_noncrossing = c("none", "rearrange"),
    cqr_right = FALSE,
    cqr_nsteps = 3L,
    cqr_first_cut = 0.1,
    cqr_later_cut = 0.05,
    cox_boundary = c("na", "error", "cap"),
    linear_backend = "auto",
    dr_backend = "auto",
    dr_workers = 1L,
    dr_warm_start = TRUE,
    dr_maxit = 100L,
    dr_tolerance = 1e-8,
    dr_precondition = TRUE,
    gpu_backend = c("cpu", "cuda"),
    gpu_precision = c("float64", "float32"),
    gpu_python = NULL,
    gpu_python_path = NULL,
    gpu_module_path = NULL,
    gpu_block_columns = 16L,
    gpu_qr_rho = 1,
    gpu_qr_maxit = 5000L,
    gpu_qr_tolerance = 1e-6,
    gpu_qr_allow_nonconvergence = FALSE) {
  nreg <- assert_scalar_integer(nreg, "nreg", 3L)
  trimming <- assert_probability(trimming, "trimming", open = FALSE)
  if (trimming >= 0.5) stop("trimming must be less than 0.5", call. = FALSE)
  reported_quantiles <- sort(unique(as.numeric(reported_quantiles)))
  if (!length(reported_quantiles) || any(!is.finite(reported_quantiles)) ||
      any(reported_quantiles <= 0 | reported_quantiles >= 1)) {
    stop("reported_quantiles must lie strictly between 0 and 1", call. = FALSE)
  }
  alpha <- assert_probability(alpha, "alpha")
  scheme_supplied <- !missing(bootstrap_scheme)
  bootstrap_scheme <- match.arg(bootstrap_scheme)
  if (!is.null(weighted_bootstrap)) {
    weighted_bootstrap <- assert_scalar_logical(
      weighted_bootstrap, "weighted_bootstrap"
    )
    mapped_scheme <- if (weighted_bootstrap) "counterfactual" else "empirical"
    if (scheme_supplied && bootstrap_scheme != mapped_scheme) {
      stop(
        "weighted_bootstrap conflicts with bootstrap_scheme",
        call. = FALSE
      )
    }
    bootstrap_scheme <- mapped_scheme
    warning(
      "weighted_bootstrap is deprecated; use bootstrap_scheme instead",
      call. = FALSE
    )
  }
  legacy_qr_shift <- assert_scalar_logical(
    legacy_qr_shift, "legacy_qr_shift"
  )
  legacy_weighted_quantile <- assert_scalar_logical(
    legacy_weighted_quantile, "legacy_weighted_quantile"
  )
  qr_precondition <- assert_scalar_logical(
    qr_precondition, "qr_precondition"
  )
  robust_se <- assert_scalar_logical(robust_se, "robust_se")
  bootstrap_max_retries <- assert_scalar_integer(
    bootstrap_max_retries, "bootstrap_max_retries", 0L
  )
  bootstrap_progress <- assert_scalar_logical(
    bootstrap_progress, "bootstrap_progress"
  )
  marginal_method <- match.arg(marginal_method)
  marginal_chunk_rows <- assert_scalar_integer(
    marginal_chunk_rows, "marginal_chunk_rows", 1L
  )
  if (length(marginal_matrix_max_mb) != 1L ||
      !is.finite(marginal_matrix_max_mb) || marginal_matrix_max_mb <= 0) {
    stop("marginal_matrix_max_mb must be a positive finite scalar", call. = FALSE)
  }
  marginal_histogram_bins <- assert_scalar_integer(
    marginal_histogram_bins, "marginal_histogram_bins", 256L
  )
  marginal_candidate_max <- assert_scalar_integer(
    marginal_candidate_max, "marginal_candidate_max", 1L
  )
  crossing_diagnostics <- assert_scalar_logical(
    crossing_diagnostics, "crossing_diagnostics"
  )
  if (length(crossing_tolerance) != 1L ||
      !is.finite(crossing_tolerance) || crossing_tolerance < 0) {
    stop("crossing_tolerance must be a nonnegative finite scalar", call. = FALSE)
  }
  quantile_noncrossing <- match.arg(quantile_noncrossing)
  cqr_right <- assert_scalar_logical(cqr_right, "cqr_right")
  cqr_nsteps <- assert_scalar_integer(cqr_nsteps, "cqr_nsteps", 3L)
  cqr_first_cut <- assert_probability(cqr_first_cut, "cqr_first_cut")
  cqr_later_cut <- assert_probability(cqr_later_cut, "cqr_later_cut")
  cox_boundary <- match.arg(cox_boundary)
  linear_backend <- match.arg(linear_backend, supported_linear_backends())
  dr_backend <- match.arg(dr_backend, supported_dr_backends())
  dr_workers <- assert_scalar_integer(dr_workers, "dr_workers", 1L)
  dr_warm_start <- assert_scalar_logical(dr_warm_start, "dr_warm_start")
  dr_precondition <- assert_scalar_logical(
    dr_precondition, "dr_precondition"
  )
  dr_maxit <- assert_scalar_integer(dr_maxit, "dr_maxit", 1L)
  if (length(dr_tolerance) != 1L || !is.finite(dr_tolerance) ||
      dr_tolerance <= 0) {
    stop("dr_tolerance must be a positive finite scalar", call. = FALSE)
  }
  gpu_backend <- match.arg(gpu_backend)
  gpu_precision <- match.arg(gpu_precision)
  validate_gpu_path <- function(path, name, directory = FALSE) {
    if (is.null(path)) return(NULL)
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path) || !file.exists(path) || (directory && !dir.exists(path))) {
      stop(name, " must be NULL or an existing path", call. = FALSE)
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  gpu_python <- validate_gpu_path(gpu_python, "gpu_python")
  gpu_python_path <- validate_gpu_path(
    gpu_python_path, "gpu_python_path", directory = TRUE
  )
  gpu_module_path <- validate_gpu_path(gpu_module_path, "gpu_module_path")
  gpu_block_columns <- assert_scalar_integer(
    gpu_block_columns, "gpu_block_columns", 1L
  )
  if (length(gpu_qr_rho) != 1L || !is.finite(gpu_qr_rho) || gpu_qr_rho <= 0) {
    stop("gpu_qr_rho must be a positive finite scalar", call. = FALSE)
  }
  gpu_qr_maxit <- assert_scalar_integer(gpu_qr_maxit, "gpu_qr_maxit", 1L)
  if (length(gpu_qr_tolerance) != 1L || !is.finite(gpu_qr_tolerance) ||
      gpu_qr_tolerance <= 0) {
    stop("gpu_qr_tolerance must be a positive finite scalar", call. = FALSE)
  }
  gpu_qr_allow_nonconvergence <- assert_scalar_logical(
    gpu_qr_allow_nonconvergence, "gpu_qr_allow_nonconvergence"
  )
  onestep_first_solver <- match.arg(onestep_first_solver)
  onestep_bandwidth <- match.arg(onestep_bandwidth)
  qr_bootstrap_engine <- match.arg(qr_bootstrap_engine)
  conditional_quantiles <- (seq_len(nreg) - 0.5) / nreg
  conditional_quantiles <- conditional_quantiles[
    conditional_quantiles >= trimming &
      conditional_quantiles <= 1 - trimming
  ]
  if (!length(conditional_quantiles)) {
    stop("trimming removed every conditional quantile", call. = FALSE)
  }
  structure(list(
    nreg = nreg,
    trimming = trimming,
    conditional_quantiles = conditional_quantiles,
    reported_quantiles = reported_quantiles,
    alpha = alpha,
    bootstrap_scheme = bootstrap_scheme,
    weighted_bootstrap = if (bootstrap_scheme == "counterfactual") {
      TRUE
    } else if (bootstrap_scheme == "empirical") {
      FALSE
    } else {
      NA
    },
    legacy_qr_shift = legacy_qr_shift,
    legacy_weighted_quantile = legacy_weighted_quantile,
    qr_precondition = qr_precondition,
    onestep_first_solver = onestep_first_solver,
    onestep_bandwidth = onestep_bandwidth,
    qr_bootstrap_engine = qr_bootstrap_engine,
    robust_se = robust_se,
    bootstrap_max_retries = bootstrap_max_retries,
    bootstrap_progress = bootstrap_progress,
    marginal_method = marginal_method,
    marginal_chunk_rows = marginal_chunk_rows,
    marginal_matrix_max_mb = as.numeric(marginal_matrix_max_mb),
    marginal_histogram_bins = marginal_histogram_bins,
    marginal_candidate_max = marginal_candidate_max,
    crossing_diagnostics = crossing_diagnostics,
    crossing_tolerance = as.numeric(crossing_tolerance),
    quantile_noncrossing = quantile_noncrossing,
    cqr_right = cqr_right,
    cqr_nsteps = cqr_nsteps,
    cqr_first_cut = cqr_first_cut,
    cqr_later_cut = cqr_later_cut,
    cox_boundary = cox_boundary,
    linear_backend = linear_backend,
    dr_backend = dr_backend,
    dr_workers = dr_workers,
    dr_warm_start = dr_warm_start,
    dr_maxit = dr_maxit,
    dr_tolerance = as.numeric(dr_tolerance),
    dr_precondition = dr_precondition,
    gpu_backend = gpu_backend,
    gpu_precision = gpu_precision,
    gpu_python = gpu_python,
    gpu_python_path = gpu_python_path,
    gpu_module_path = gpu_module_path,
    gpu_block_columns = gpu_block_columns,
    gpu_qr_rho = as.numeric(gpu_qr_rho),
    gpu_qr_maxit = gpu_qr_maxit,
    gpu_qr_tolerance = as.numeric(gpu_qr_tolerance),
    gpu_qr_allow_nonconvergence = gpu_qr_allow_nonconvergence
  ), class = "cf_control")
}

validate_cf_control <- function(control) {
  if (!inherits(control, "cf_control")) {
    stop("control must be created by cf_control()", call. = FALSE)
  }
  required <- c(
    "nreg", "trimming", "conditional_quantiles", "reported_quantiles",
    "alpha", "bootstrap_scheme", "legacy_qr_shift",
    "legacy_weighted_quantile", "qr_precondition",
    "onestep_first_solver", "onestep_bandwidth", "qr_bootstrap_engine",
    "robust_se", "bootstrap_max_retries", "bootstrap_progress",
    "marginal_method", "marginal_chunk_rows", "marginal_matrix_max_mb",
    "marginal_histogram_bins", "marginal_candidate_max",
    "crossing_diagnostics", "crossing_tolerance", "quantile_noncrossing",
    "cqr_right", "cqr_nsteps", "cqr_first_cut", "cqr_later_cut",
    "cox_boundary", "linear_backend",
    "dr_backend", "dr_workers",
    "dr_warm_start", "dr_maxit", "dr_tolerance", "dr_precondition",
    "gpu_backend", "gpu_precision", "gpu_python", "gpu_python_path",
    "gpu_module_path", "gpu_block_columns", "gpu_qr_rho",
    "gpu_qr_maxit", "gpu_qr_tolerance", "gpu_qr_allow_nonconvergence"
  )
  missing_fields <- setdiff(required, names(control))
  if (length(missing_fields)) {
    stop(
      "control is missing field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(control)
}
