# Counterfactual estimands and functional inference follow Chernozhukov,
# Fernandez-Val, and Melly (2013), doi:10.3982/ECTA10582. See
# inst/provenance/METHODS.md for the distinction between method and code source.
validate_model_solver <- function(model, solver) {
  model <- match.arg(model, supported_cf_models())
  if (model == "qr") {
    if (is.null(solver)) solver <- "auto"
    solver <- match.arg(solver, supported_qr_solvers())
  } else if (model == "cqr") {
    if (is.null(solver)) solver <- "auto"
    solver <- validate_cqr_solver(solver)
  } else {
    internal_na <- length(solver) == 1L && is.na(solver)
    if (!is.null(solver) && !internal_na) {
      stop("solver is only valid when model = 'qr'", call. = FALSE)
    }
    solver <- NA_character_
  }
  list(model = model, solver = solver)
}

execution_backend_plan <- function(model, solver, control) {
  dr_models <- c("logit", "probit", "cloglog")
  if (identical(control$dr_backend, "cuda") && !model %in% dr_models) {
    stop(
      "dr_backend='cuda' is only supported for logit, probit, and cloglog; ",
      "model='", model, "' does not use the CUDA DR fitter",
      call. = FALSE
    )
  }
  if (identical(control$gpu_backend, "cuda") && model == "cox") {
    stop(
      "gpu_backend='cuda' is not implemented for model='cox'; use ",
      "gpu_backend='cpu' so metadata matches the executed Cox path",
      call. = FALSE
    )
  }
  fit_device <- if (identical(solver, "cuda_admm") ||
                    (model %in% dr_models &&
                     identical(control$dr_backend, "cuda"))) {
    "cuda"
  } else {
    "cpu"
  }
  prediction_device <- if (
    identical(control$gpu_backend, "cuda") && model != "cox"
  ) "cuda" else "cpu"
  marginalization_device <- if (prediction_device == "cpu") {
    "cpu"
  } else if (model %in% c("qr", "cqr") &&
             !isTRUE(control$legacy_weighted_quantile) &&
             !identical(control$marginal_method, "chunked")) {
    "cuda"
  } else if (model %in% c("logit", "probit", "cloglog", "lpm")) {
    "cuda"
  } else {
    "cuda+cpu"
  }
  requested_backend <- if (model %in% c("qr", "cqr")) {
    solver
  } else if (model %in% c("loc", "locsca", "lpm")) {
    control$linear_backend
  } else if (model == "cox") {
    "survival::coxph.fit(breslow)"
  } else {
    control$dr_backend
  }
  list(
    fit_device = fit_device,
    prediction_device = prediction_device,
    marginalization_device = marginalization_device,
    requested_backend = requested_backend,
    any_cuda = fit_device == "cuda" || prediction_device == "cuda"
  )
}

fit_backend_name <- function(fit) {
  fit$solver %||% fit$backend %||% "unknown"
}

resolved_fit_device <- function(point, planned) {
  threshold_backends <- unique(unlist(lapply(point$fits, function(fit) {
    fit$threshold_backend %||% character()
  })))
  if (length(threshold_backends)) {
    uses_cuda <- "cuda" %in% threshold_backends
    uses_cpu <- any(!threshold_backends %in% c("cuda", "analytic"))
    if (uses_cuda && uses_cpu) return("cuda+cpu_fallback")
    if (uses_cuda) return("cuda")
    return("cpu")
  }
  planned
}

fit_two_group_models <- function(
    prepared, model, solver, control, point_workers, point_seed = NULL) {
  point_workers <- min(assert_scalar_integer(point_workers, "point_workers", 1L), 2L)
  tasks <- list(
    list(
      X = prepared$X0, y = prepared$y0, weights = prepared$w0,
      censoring = prepared$censoring0, event = prepared$event0,
      model = model, solver = solver, control = control
    ),
    list(
      X = prepared$X1, y = prepared$y1, weights = prepared$w1,
      censoring = prepared$censoring1, event = prepared$event1,
      model = model, solver = solver, control = control
    )
  )
  if (point_workers == 1L) {
    if (!is.null(point_seed)) set.seed(point_seed)
    return(lapply(tasks, fit_group_task))
  }
  cluster <- parallel::makeCluster(point_workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  package_worker_init(cluster)
  if (!is.null(point_seed)) {
    parallel::clusterSetRNGStream(cluster, iseed = point_seed)
  }
  parallel::parLapply(cluster, tasks, function(task) {
    fit_task <- get("fit_group_task", envir = asNamespace("scalableCounterfactual"))
    fit_task(task)
  })
}

evaluate_decomposition <- function(prepared, fit0, fit1, control) {
  marginal_resource <- measure_resources(function() {
    marginal_w0 <- if (isTRUE(control$legacy_weighted_quantile)) {
      prepared$w0 / prepared$n
    } else {
      prepared$w0
    }
    marginal_w1 <- if (isTRUE(control$legacy_weighted_quantile)) {
      prepared$w1 / prepared$n
    } else {
      prepared$w1
    }
    fitted0 <- marginal_quantiles(
      fit0, prepared$X0, marginal_w0,
      control$reported_quantiles, control
    )
    counterfactual01 <- marginal_quantiles(
      fit0, prepared$X1, marginal_w1,
      control$reported_quantiles, control
    )
    fitted1 <- marginal_quantiles(
      fit1, prepared$X1, marginal_w1,
      control$reported_quantiles, control
    )
    marginal_diagnostics <- do.call(rbind, Map(
      function(distribution, values) {
        diagnostics <- attr(values, "marginal_diagnostics")
        data.frame(
          distribution = distribution,
          method = diagnostics$method,
          passes = diagnostics$passes,
          histogram_bins = diagnostics$histogram_bins,
          candidate_draws = diagnostics$candidate_draws,
          estimated_matrix_mb = diagnostics$estimated_matrix_mb,
          boundary_quantiles = diagnostics$boundary_quantiles %||% 0L,
          identified_cdf_max = diagnostics$identified_cdf_max %||% NA_real_,
          stringsAsFactors = FALSE
        )
      },
      c("reference_model_X0", "reference_model_X1_counterfactual",
        "comparison_model_X1"),
      list(fitted0, counterfactual01, fitted1)
    ))
    fitted0 <- as.numeric(fitted0)
    counterfactual01 <- as.numeric(counterfactual01)
    fitted1 <- as.numeric(fitted1)
    if (inherits(fit0, "cf_cox_fit")) {
      observed0 <- weighted_km_quantiles(
        prepared$y0, prepared$event0, marginal_w0,
        control$reported_quantiles, control$cox_boundary
      )
      observed1 <- weighted_km_quantiles(
        prepared$y1, prepared$event1, marginal_w1,
        control$reported_quantiles, control$cox_boundary
      )
    } else {
      observed0 <- weighted_quantile(
        prepared$y0, marginal_w0, control$reported_quantiles,
        legacy = control$legacy_weighted_quantile
      )
      observed1 <- weighted_quantile(
        prepared$y1, marginal_w1, control$reported_quantiles,
        legacy = control$legacy_weighted_quantile
      )
    }
    structure_effect <- fitted1 - counterfactual01
    composition_effect <- counterfactual01 - fitted0
    total_effect <- fitted1 - fitted0
    list(
      effects = rbind(
        structure = structure_effect,
        composition = composition_effect,
        total = total_effect
      ),
      diagnostics = rbind(
        reference_observed = observed0,
        reference_model = fitted0,
        counterfactual_reference_structure_comparison_X = counterfactual01,
        comparison_observed = observed1,
        comparison_model = fitted1
      ),
      marginal_diagnostics = marginal_diagnostics
    )
  })
  effects <- marginal_resource$value$effects
  residual <- effects["total", ] -
    effects["structure", ] - effects["composition", ]
  finite_residual <- is.finite(residual)
  if (any(finite_residual) &&
      max(abs(residual[finite_residual])) > 1e-10) {
    stop("decomposition identity failed", call. = FALSE)
  }
  list(
    effects = effects,
    diagnostics = marginal_resource$value$diagnostics,
    marginal_diagnostics = marginal_resource$value$marginal_diagnostics,
    identity_residual = residual,
    resource = marginal_resource
  )
}

estimate_point_prepared <- function(
    prepared,
    model,
    solver,
    control,
    point_workers = 1L,
    point_seed = NULL,
    keep_fits = TRUE) {
  model_solver <- validate_model_solver(model, solver)
  model <- model_solver$model
  solver <- model_solver$solver
  started <- proc.time()[["elapsed"]]
  group_results <- fit_two_group_models(
    prepared, model, solver, control, point_workers, point_seed
  )
  fit0 <- group_results[[1L]]$fit
  fit1 <- group_results[[2L]]$fit
  evaluated <- evaluate_decomposition(prepared, fit0, fit1, control)
  marginal_resource <- evaluated$resource
  effects <- evaluated$effects
  residual <- evaluated$identity_residual
  elapsed <- proc.time()[["elapsed"]] - started
  resources <- data.frame(
    phase = c("group0_fit", "group1_fit", "marginalization"),
    elapsed_seconds = c(
      group_results[[1L]]$elapsed_seconds,
      group_results[[2L]]$elapsed_seconds,
      marginal_resource$elapsed_seconds
    ),
    peak_r_heap_mb = c(
      group_results[[1L]]$peak_r_heap_mb,
      group_results[[2L]]$peak_r_heap_mb,
      marginal_resource$peak_r_heap_mb
    ),
    stringsAsFactors = FALSE
  )
  result <- list(
    model = model,
    solver = solver,
    quantiles = control$reported_quantiles,
    effects = effects,
    diagnostics = evaluated$diagnostics,
    marginal_diagnostics = evaluated$marginal_diagnostics,
    identity_residual = residual,
    elapsed_seconds = unname(elapsed),
    resources = resources,
    warnings = unique(unlist(lapply(
      list(fit0, fit1),
      function(fit) fit$warnings
    )))
  )
  if (keep_fits) result$fits <- list(group0 = fit0, group1 = fit1)
  result
}

#' Estimate a counterfactual quantile decomposition
#'
#' @param formula Outcome and covariate formula.
#' @param data Data frame containing analysis variables.
#' @param group Binary group column name or vector. Zero is the reference group
#'   and one is the comparison group.
#' @param weights Optional positive sampling-weight column name or vector.
#' @param model Conditional-distribution model: `qr`, `cqr`, `loc`, `locsca`,
#'   `cox`, `logit`, `probit`, `cloglog`, or `lpm`.
#' @param solver QR solver: `br`, `fn`, `pfn`, `qfnb`, `pfnb`, `proqreg`,
#'   `profn`, `onestep`, or `auto`. CQR accepts the exact solvers `br`, `fn`,
#'   `pfn`, `qfnb`, `pfnb`, and `auto`. Leave `NULL` for other models.
#' @param censoring Censoring-point column, vector, or scalar for `model =
#'   "cqr"`.
#' @param event Optional 0/1 event-status column or vector for `model = "cox"`.
#'   If omitted, every outcome is treated as an observed event.
#' @param control A [cf_control()] object.
#' @param bootstrap_reps Number of bootstrap replications.
#' @param point_workers Workers used for the two group-specific point fits. If
#'   `bootstrap_workers = 1`, the standard bootstrap also reuses up to two of
#'   these workers for the two group fits within each replication.
#' @param bootstrap_workers Independent bootstrap-replication workers.
#' @param checkpoint_dir Optional persistent checkpoint directory.
#' @param seed Random seed.
#' @return An object of class `cfdecomp`.
#' @examples
#' set.seed(11)
#' example_data <- data.frame(
#'   outcome = rnorm(120),
#'   x = rnorm(120),
#'   group = rep(0:1, each = 60),
#'   weight = runif(120, 0.5, 1.5)
#' )
#' example_fit <- counterfactual_decompose(
#'   outcome ~ x,
#'   data = example_data,
#'   group = "group",
#'   weights = "weight",
#'   model = "qr",
#'   solver = "fn",
#'   control = cf_control(
#'     nreg = 5,
#'     reported_quantiles = c(0.25, 0.5, 0.75),
#'     crossing_diagnostics = FALSE
#'   ),
#'   seed = 11
#' )
#' as.data.frame(example_fit)
#' @export
counterfactual_decompose <- function(
    formula,
    data,
    group,
    weights = NULL,
    model = "qr",
    solver = NULL,
    control = cf_control(),
    bootstrap_reps = 0L,
    point_workers = 1L,
    bootstrap_workers = 1L,
    checkpoint_dir = NULL,
    seed = 20260719L,
    censoring = NULL,
    event = NULL) {
  validate_cf_control(control)
  model_solver <- validate_model_solver(model, solver)
  group_definition <- if (is.character(group) && length(group) == 1L) {
    group
  } else {
    "<supplied vector>"
  }
  weights_definition <- if (is.null(weights)) {
    "<equal weights>"
  } else if (is.character(weights) && length(weights) == 1L) {
    weights
  } else {
    "<supplied vector>"
  }
  describe_vector_argument <- function(value, absent) {
    if (is.null(value)) return(absent)
    if (is.character(value) && length(value) == 1L) return(value)
    if (length(value) == 1L) return(as.character(value))
    "<supplied vector>"
  }
  censoring_definition <- describe_vector_argument(censoring, "<not used>")
  event_definition <- describe_vector_argument(event, if (
    model_solver$model == "cox"
  ) "<all events>" else "<not used>")
  bootstrap_reps <- assert_scalar_integer(
    bootstrap_reps, "bootstrap_reps", 0L
  )
  point_workers <- assert_scalar_integer(point_workers, "point_workers", 1L)
  bootstrap_workers <- assert_scalar_integer(
    bootstrap_workers, "bootstrap_workers", 1L
  )
  backend_plan <- execution_backend_plan(
    model_solver$model, model_solver$solver, control
  )
  if (backend_plan$fit_device == "cuda" && point_workers > 1L) {
    stop(
      "CUDA conditional fitting requires point_workers=1; use GPU column ",
      "batching instead of multiple CUDA fitting processes",
      call. = FALSE
    )
  }
  if (backend_plan$any_cuda && bootstrap_reps > 0L &&
      bootstrap_workers > 1L) {
    stop(
      "CUDA bootstrap currently requires bootstrap_workers=1 to avoid ",
      "multiple R processes oversubscribing one device",
      call. = FALSE
    )
  }
  if (control$dr_workers > 1L &&
      (point_workers > 1L || (bootstrap_reps > 0L && bootstrap_workers > 1L))) {
    stop(
      "dr_workers > 1 cannot be combined with point_workers > 1 or ",
      "bootstrap_workers > 1; choose one parallel execution layer",
      call. = FALSE
    )
  }
  seed <- assert_scalar_integer(seed, "seed", 1L)
  prepared <- prepare_cf_data(
    formula, data, group, weights,
    model = model_solver$model, censoring = censoring, event = event
  )
  if (length(prepared$dropped_design_columns)) {
    warning(
      "Dropped columns that were unidentified in at least one group: ",
      paste(prepared$dropped_design_columns, collapse = ", "),
      call. = FALSE
    )
  }
  point <- estimate_point_prepared(
    prepared,
    model_solver$model,
    model_solver$solver,
    control,
    point_workers = point_workers,
    point_seed = seed,
    keep_fits = TRUE
  )
  if (is_quantile_process_model(model_solver$model) &&
      isTRUE(control$crossing_diagnostics)) {
    crossing_resource <- measure_resources(function() {
      evaluate_qr_crossing_diagnostics(
        prepared, point$fits$group0, point$fits$group1, control
      )
    })
    point$crossing_diagnostics <- crossing_resource$value
    point$elapsed_seconds <- point$elapsed_seconds +
      crossing_resource$elapsed_seconds
    point$resources <- rbind(point$resources, data.frame(
      phase = "crossing_diagnostics",
      elapsed_seconds = crossing_resource$elapsed_seconds,
      peak_r_heap_mb = crossing_resource$peak_r_heap_mb,
      stringsAsFactors = FALSE
    ))
  } else {
    point$crossing_diagnostics <- data.frame()
  }
  if (length(point$warnings)) {
    warning(
      "Conditional-model warnings: ",
      paste(point$warnings, collapse = " | "),
      call. = FALSE
    )
  }
  bootstrap <- NULL
  inference <- point_effect_table(point, bootstrap_reps, control$alpha)
  if (bootstrap_reps > 0L) {
    if (is.null(checkpoint_dir)) {
      checkpoint_dir <- file.path(tempdir(), "scalableCounterfactual_checkpoints")
    }
    bootstrap <- run_cf_bootstrap(
      prepared = prepared,
      point = point,
      model = model_solver$model,
      solver = model_solver$solver,
      control = control,
      reps = bootstrap_reps,
      workers = bootstrap_workers,
      draw_point_workers = point_workers,
      checkpoint_dir = checkpoint_dir,
      seed = seed
    )
    if (nrow(bootstrap$failures)) {
      warning(
        "Bootstrap used ", nrow(bootstrap$failures),
        " failed attempt(s) before obtaining ", bootstrap_reps,
        " successful replications; inspect bootstrap_failures.csv",
        call. = FALSE
      )
    }
    inference <- bootstrap_inference(point, bootstrap, control)
  }
  metadata <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    package = "scalableCounterfactual",
    package_version = as.character(utils::packageVersion("scalableCounterfactual")),
    output_schema_version = .cf_output_schema_version,
    extension_registry_fingerprint = active_extension_fingerprint(
      model_solver$model, model_solver$solver, control, point
    ),
    formula = paste(deparse(formula), collapse = " "),
    group = group_definition,
    weights = weights_definition,
    censoring = censoring_definition,
    event = event_definition,
    model = model_solver$model,
    solver = model_solver$solver,
    fit_device = resolved_fit_device(point, backend_plan$fit_device),
    prediction_device = backend_plan$prediction_device,
    marginalization_device = backend_plan$marginalization_device,
    requested_backend = backend_plan$requested_backend,
    resolved_backend_group0 = fit_backend_name(point$fits$group0),
    resolved_backend_group1 = fit_backend_name(point$fits$group1),
    nreg = control$nreg,
    trimming = control$trimming,
    legacy_qr_shift = control$legacy_qr_shift,
    legacy_weighted_quantile = control$legacy_weighted_quantile,
    qr_precondition = control$qr_precondition,
    marginal_method_requested = control$marginal_method,
    marginal_methods_resolved = paste(
      unique(point$marginal_diagnostics$method), collapse = ","
    ),
    marginal_chunk_rows = control$marginal_chunk_rows,
    marginal_matrix_max_mb = control$marginal_matrix_max_mb,
    marginal_histogram_bins = control$marginal_histogram_bins,
    marginal_candidate_max = control$marginal_candidate_max,
    crossing_diagnostics = control$crossing_diagnostics,
    crossing_tolerance = control$crossing_tolerance,
    quantile_noncrossing = control$quantile_noncrossing,
    cqr_right = control$cqr_right,
    cqr_nsteps = control$cqr_nsteps,
    cqr_first_cut = control$cqr_first_cut,
    cqr_later_cut = control$cqr_later_cut,
    cox_boundary = control$cox_boundary,
    crossing_rows_max_share = if (nrow(point$crossing_diagnostics)) {
      max(point$crossing_diagnostics$crossing_row_share)
    } else {
      NA_real_
    },
    onestep_first_solver = control$onestep_first_solver,
    onestep_bandwidth = control$onestep_bandwidth,
    linear_backend_requested = control$linear_backend,
    dr_backend_requested = control$dr_backend,
    dr_workers = control$dr_workers,
    dr_warm_start = control$dr_warm_start,
    dr_maxit = control$dr_maxit,
    dr_tolerance = control$dr_tolerance,
    dr_precondition = control$dr_precondition,
    gpu_backend = control$gpu_backend,
    gpu_backend_requested = control$gpu_backend,
    gpu_precision = control$gpu_precision,
    gpu_block_columns = control$gpu_block_columns,
    gpu_qr_rho = control$gpu_qr_rho,
    gpu_qr_maxit = control$gpu_qr_maxit,
    gpu_qr_tolerance = control$gpu_qr_tolerance,
    gpu_qr_allow_nonconvergence = control$gpu_qr_allow_nonconvergence,
    dr_warm_start_effective_group0 = if (
      model_solver$model %in% c("logit", "probit", "cloglog")
    ) point$fits$group0$warm_start else NA,
    dr_warm_start_effective_group1 = if (
      model_solver$model %in% c("logit", "probit", "cloglog")
    ) point$fits$group1$warm_start else NA,
    dr_threshold_workers_group0 = if (
      is_distribution_regression_model(model_solver$model)
    ) point$fits$group0$threshold_workers else NA_integer_,
    dr_threshold_workers_group1 = if (
      is_distribution_regression_model(model_solver$model)
    ) point$fits$group1$threshold_workers else NA_integer_,
    dr_boundary_thresholds_group0 = if (
      is_distribution_regression_model(model_solver$model)
    ) sum(point$fits$group0$convergence_flag == 2L) else NA_integer_,
    dr_boundary_thresholds_group1 = if (
      is_distribution_regression_model(model_solver$model)
    ) sum(point$fits$group1$convergence_flag == 2L) else NA_integer_,
    dr_threshold_backends_group0 = if (
      is_distribution_regression_model(model_solver$model)
    ) paste(point$fits$group0$threshold_backend, collapse = ",") else NA_character_,
    dr_threshold_backends_group1 = if (
      is_distribution_regression_model(model_solver$model)
    ) paste(point$fits$group1$threshold_backend, collapse = ",") else NA_character_,
    dr_cuda_fallbacks_group0 = if (
      is_distribution_regression_model(model_solver$model)
    ) sum(!is.na(point$fits$group0$threshold_fallback_reason)) else 0L,
    dr_cuda_fallbacks_group1 = if (
      is_distribution_regression_model(model_solver$model)
    ) sum(!is.na(point$fits$group1$threshold_fallback_reason)) else 0L,
    cox_boundary_quantiles = if (model_solver$model == "cox") {
      sum(point$marginal_diagnostics$boundary_quantiles)
    } else 0L,
    conditional_backend_group0 = if (is_quantile_process_model(model_solver$model)) {
      point$fits$group0$solver
    } else {
      point$fits$group0$backend
    },
    conditional_backend_group1 = if (is_quantile_process_model(model_solver$model)) {
      point$fits$group1$solver
    } else {
      point$fits$group1$backend
    },
    solver_group0 = if (is_quantile_process_model(model_solver$model)) {
      point$fits$group0$solver
    } else {
      NA_character_
    },
    solver_group1 = if (is_quantile_process_model(model_solver$model)) {
      point$fits$group1$solver
    } else {
      NA_character_
    },
    solver_group0_exact = if (is_quantile_process_model(model_solver$model)) {
      point$fits$group0$solver_exact
    } else {
      NA
    },
    solver_group1_exact = if (is_quantile_process_model(model_solver$model)) {
      point$fits$group1$solver_exact
    } else {
      NA
    },
    solver_group0_process_fallbacks = if (model_solver$model == "qr") {
      length(point$fits$group0$process_fallback_taus)
    } else {
      0L
    },
    solver_group1_process_fallbacks = if (model_solver$model == "qr") {
      length(point$fits$group1$process_fallback_taus)
    } else {
      0L
    },
    solver_group0_process_fallback_reason = if (model_solver$model == "qr") {
      point$fits$group0$process_fallback_reason
    } else {
      NA_character_
    },
    solver_group1_process_fallback_reason = if (model_solver$model == "qr") {
      point$fits$group1$process_fallback_reason
    } else {
      NA_character_
    },
    stata_source_version = if (
      model_solver$model == "qr" &&
        all(vapply(point$fits, function(x) x$solver == "onestep", logical(1L)))
    ) {
      point$fits$group0$stata_source_version
    } else {
      NA_character_
    },
    stata_source_commit = if (
      model_solver$model == "qr" &&
        all(vapply(point$fits, function(x) x$solver == "onestep", logical(1L)))
    ) {
      point$fits$group0$stata_source_commit
    } else {
      NA_character_
    },
    onestep_implementation_version = if (
      model_solver$model == "qr" && point$fits$group0$solver == "onestep"
    ) {
      point$fits$group0$onestep_implementation_version
    } else {
      NA_character_
    },
    group0_onestep_fallbacks = if (
      model_solver$model == "qr" && point$fits$group0$solver == "onestep"
    ) {
      length(point$fits$group0$onestep_fallback_taus)
    } else {
      0L
    },
    group1_onestep_fallbacks = if (
      model_solver$model == "qr" && point$fits$group1$solver == "onestep"
    ) {
      length(point$fits$group1$onestep_fallback_taus)
    } else {
      0L
    },
    bootstrap_reps = bootstrap_reps,
    bootstrap_max_retries = control$bootstrap_max_retries,
    bootstrap_progress = control$bootstrap_progress,
    bootstrap_retries_used = if (is.null(bootstrap)) {
      0L
    } else {
      bootstrap$retries
    },
    bootstrap_failed_attempts = if (is.null(bootstrap)) {
      0L
    } else {
      nrow(bootstrap$failures)
    },
    bootstrap_scheme = control$bootstrap_scheme,
    weighted_bootstrap = control$weighted_bootstrap,
    qr_bootstrap_engine_requested = control$qr_bootstrap_engine,
    qr_bootstrap_engine = if (is.null(bootstrap)) {
      NA_character_
    } else {
      bootstrap$engine
    },
    point_workers = min(point_workers, 2L),
    bootstrap_workers = bootstrap_workers,
    bootstrap_point_workers = if (is.null(bootstrap)) {
      NA_integer_
    } else {
      bootstrap$draw_point_workers
    },
    seed = seed,
    analysis_rows = prepared$n,
    group0_rows = prepared$n0,
    group1_rows = prepared$n1,
    omitted_rows = prepared$omitted_rows,
    point_elapsed_seconds = point$elapsed_seconds,
    conditional_model_warnings = paste(point$warnings, collapse = " | "),
    checkpoint_root = if (bootstrap_reps > 0L) normalizePath(
      checkpoint_dir, winslash = "/", mustWork = FALSE
    ) else NA_character_,
    checkpoint_dir = if (bootstrap_reps > 0L) normalizePath(
      bootstrap$checkpoint_dir, winslash = "/", mustWork = FALSE
    ) else NA_character_,
    bootstrap_signature = if (is.null(bootstrap)) {
      NA_character_
    } else {
      bootstrap$signature
    },
    bootstrap_data_fingerprint = if (is.null(bootstrap)) {
      NA_character_
    } else {
      bootstrap$data_fingerprint
    },
    bootstrap_warning_count = if (is.null(bootstrap)) {
      0L
    } else {
      length(bootstrap$warnings)
    },
    bootstrap_warnings = if (is.null(bootstrap)) {
      NA_character_
    } else {
      paste(bootstrap$warnings, collapse = " | ")
    },
    bootstrap_dropped_design_columns = if (is.null(bootstrap)) {
      NA_character_
    } else {
      paste(bootstrap$dropped_design_columns, collapse = ",")
    },
    R_version = R.version.string,
    quantreg_version = as.character(utils::packageVersion("quantreg")),
    fastglm_version = if (requireNamespace("fastglm", quietly = TRUE)) {
      as.character(utils::packageVersion("fastglm"))
    } else {
      NA_character_
    },
    speedglm_version = if (requireNamespace("speedglm", quietly = TRUE)) {
      as.character(utils::packageVersion("speedglm"))
    } else {
      NA_character_
    },
    design_columns = paste(prepared$design_columns, collapse = ","),
    dropped_design_columns = paste(
      prepared$dropped_design_columns,
      collapse = ","
    )
  )
  object <- structure(list(
    call = match.call(),
    metadata = metadata,
    control = control,
    point = point,
    bootstrap = bootstrap,
    results = inference
  ), class = "cfdecomp")
  object$functional_tests <- if (bootstrap_reps > 1L) {
    functional_effect_tests(object)
  } else {
    data.frame()
  }
  object
}
