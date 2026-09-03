# Bootstrap estimands follow doi:10.3982/ECTA10582. Solver-specific source
# relationships are recorded in inst/provenance/METHODS.md.
.cf_bootstrap_worker_state <- new.env(parent = emptyenv())

bootstrap_data_fingerprint <- function(prepared) {
  object_md5(list(
    schema_version = 2L,
    design_columns = prepared$design_columns,
    dropped_design_columns = prepared$dropped_design_columns,
    X0 = prepared$X0,
    y0 = prepared$y0,
    w0 = prepared$w0,
    censoring0 = prepared$censoring0,
    event0 = prepared$event0,
    X1 = prepared$X1,
    y1 = prepared$y1,
    w1 = prepared$w1,
    censoring1 = prepared$censoring1,
    event1 = prepared$event1
  ))
}

bootstrap_runtime_package_version <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

bootstrap_gpu_runtime_identity <- function(control, model, solver) {
  gpu_active <- identical(control$gpu_backend, "cuda") ||
    identical(control$dr_backend, "cuda") || identical(solver, "cuda_admm")
  if (!gpu_active) return(NULL)
  module_path <- control$gpu_module_path
  if (is.null(module_path)) {
    module_path <- tryCatch(gpu_python_module_path(), error = function(error) NULL)
  }
  module_hash <- if (!is.null(module_path) && file.exists(module_path)) {
    object_md5(readBin(module_path, what = "raw", n = file.info(module_path)$size))
  } else {
    NA_character_
  }
  if (exists("gpu_runtime_metadata", mode = "function", inherits = TRUE)) {
    metadata <- tryCatch(
      gpu_runtime_metadata(control = control),
      error = function(error) NULL
    )
    if (!is.null(metadata)) {
      metadata$module_file <- metadata$module_file %||% module_path
      metadata$module_content_md5 <- module_hash
      return(metadata)
    }
  }
  list(
    backend = control$gpu_backend,
    model = model,
    solver = solver,
    python = control$gpu_python %||% NA_character_,
    python_path = control$gpu_python_path %||% NA_character_,
    module_file = module_path %||% NA_character_,
    module_hash = module_hash,
    reticulate_version = bootstrap_runtime_package_version("reticulate")
  )
}

bootstrap_runtime_identity <- function(model, solver, control, point = NULL) {
  requested_backends <- c(
    solver %||% character(), control$linear_backend, control$dr_backend
  )
  resolved_backends <- if (!is.null(point) && length(point$fits)) {
    unique(unlist(lapply(point$fits, function(fit) {
      c(
        fit$solver %||% character(), fit$backend %||% character(),
        fit$selection_backend %||% character()
      )
    }), use.names = FALSE))
  } else {
    character()
  }
  fitted_backends <- unique(c(requested_backends, resolved_backends))
  versions <- list()
  if (isTRUE(control$legacy_weighted_quantile)) {
    versions$Hmisc <- bootstrap_runtime_package_version("Hmisc")
  }
  if (is_quantile_process_model(model)) {
    versions$quantreg <- bootstrap_runtime_package_version("quantreg")
  }
  if (identical(model, "cox")) {
    versions$survival <- bootstrap_runtime_package_version("survival")
  }
  uses_auto_dr <- identical(control$dr_backend, "auto") &&
    model %in% c("cqr", "logit", "probit", "cloglog")
  needs_fastglm_identity <- uses_auto_dr ||
    any(grepl("fastglm", fitted_backends, fixed = TRUE))
  fastglm_available <- if (needs_fastglm_identity) {
    requireNamespace("fastglm", quietly = TRUE)
  } else {
    NA
  }
  if (needs_fastglm_identity) {
    versions$fastglm <- bootstrap_runtime_package_version("fastglm")
  }
  if (any(grepl("speedglm", fitted_backends, fixed = TRUE))) {
    versions$speedglm <- bootstrap_runtime_package_version("speedglm")
  }
  list(
    R = R.version.string,
    platform = R.version$platform,
    os_type = .Platform$OS.type,
    blas = unname(extSoftVersion()[["BLAS"]] %||% NA_character_),
    lapack = as.character(base::La_version()),
    packages = versions,
    backend_availability = if (uses_auto_dr) {
      list(fastglm = fastglm_available)
    } else {
      list()
    },
    gpu = bootstrap_gpu_runtime_identity(control, model, solver)
  )
}

bootstrap_signature <- function(
    prepared, model, solver, control, seed, bootstrap_engine,
    draw_point_workers = 1L,
    data_fingerprint = NULL, point = NULL) {
  if (is.null(data_fingerprint)) {
    data_fingerprint <- bootstrap_data_fingerprint(prepared)
  }
  object_md5(list(
    schema_version = 14L,
    package_version = as.character(
      utils::packageVersion("scalableCounterfactual")
    ),
    data_fingerprint = data_fingerprint,
    extension_registry = active_extension_fingerprint(
      model, solver, control, point
    ),
    runtime = bootstrap_runtime_identity(model, solver, control, point),
    model = model,
    solver = solver,
    nreg = control$nreg,
    trimming = control$trimming,
    quantiles = control$reported_quantiles,
    conditional_quantiles = control$conditional_quantiles,
    full_conditional_quantiles = control$full_conditional_quantiles,
    bootstrap_scheme = control$bootstrap_scheme,
    legacy_qr_shift = control$legacy_qr_shift,
    legacy_weighted_quantile = control$legacy_weighted_quantile,
    qr_precondition = control$qr_precondition,
    onestep_first_solver = control$onestep_first_solver,
    onestep_bandwidth = control$onestep_bandwidth,
    qr_bootstrap_engine = bootstrap_engine,
    linear_backend = control$linear_backend,
    dr_backend = control$dr_backend,
    dr_workers = control$dr_workers,
    dr_warm_start = control$dr_warm_start,
    dr_maxit = control$dr_maxit,
    dr_tolerance = control$dr_tolerance,
    dr_precondition = control$dr_precondition,
    bootstrap_max_retries = control$bootstrap_max_retries,
    marginal_method = control$marginal_method,
    marginal_chunk_rows = control$marginal_chunk_rows,
    marginal_matrix_max_mb = control$marginal_matrix_max_mb,
    marginal_histogram_bins = control$marginal_histogram_bins,
    marginal_candidate_max = control$marginal_candidate_max,
    quantile_noncrossing = control$quantile_noncrossing,
    dr_noncrossing = control$dr_noncrossing,
    cqr_right = control$cqr_right,
    cqr_nsteps = control$cqr_nsteps,
    cqr_first_cut = control$cqr_first_cut,
    cqr_later_cut = control$cqr_later_cut,
    cox_boundary = control$cox_boundary,
    gpu_backend = control$gpu_backend,
    gpu_precision = control$gpu_precision,
    gpu_block_columns = control$gpu_block_columns,
    gpu_qr_rho = control$gpu_qr_rho,
    gpu_qr_maxit = control$gpu_qr_maxit,
    gpu_qr_tolerance = control$gpu_qr_tolerance,
    gpu_qr_allow_nonconvergence = control$gpu_qr_allow_nonconvergence,
    draw_point_workers = draw_point_workers,
    seed = seed
  ))
}

bootstrap_multipliers <- function(
    n, scheme = c("counterfactual", "empirical", "multiplier")) {
  scheme <- match.arg(scheme)
  if (scheme == "multiplier") return(stats::rexp(n, rate = 1))
  probabilities <- if (scheme == "counterfactual") {
    stats::rexp(n, rate = 1)
  } else {
    NULL
  }
  tabulate(
    sample.int(n, n, replace = TRUE, prob = probabilities),
    nbins = n
  )
}

resolve_qr_bootstrap_engine <- function(
    requested, model, point, bootstrap_scheme) {
  if (model != "qr") {
    if (requested %in% c("xy_preprocess", "onestep")) {
      stop(requested, " is available only for QR models", call. = FALSE)
    }
    return("standard")
  }
  supported <- c("profn", "proqreg")
  fitted_solvers <- vapply(
    point$fits,
    function(fit) fit$solver,
    character(1L)
  )
  if (requested == "auto") {
    return(if (all(fitted_solvers == "onestep")) {
      "onestep"
    } else if (all(fitted_solvers %in% supported) &&
               bootstrap_scheme != "multiplier") {
      "xy_preprocess"
    } else {
      "standard"
    })
  }
  if (requested == "xy_preprocess" &&
      !all(fitted_solvers %in% supported)) {
    stop(
      "xy_preprocess requires profn or proqreg in both group fits; resolved solvers were ",
      paste(fitted_solvers, collapse = ", "),
      call. = FALSE
    )
  }
  if (requested == "xy_preprocess" && bootstrap_scheme == "multiplier") {
    stop(
      "xy_preprocess does not support the multiplier bootstrap scheme",
      call. = FALSE
    )
  }
  if (requested == "onestep" &&
      !all(fitted_solvers == "onestep")) {
    stop(
      "onestep bootstrap requires onestep in both group fits; resolved solvers were ",
      paste(fitted_solvers, collapse = ", "),
      call. = FALSE
    )
  }
  requested
}

# Native R transcription of qrprocess.ado 1.1.3, rq_boot_1step(), lines
# 2362-2455. In particular, the inverse Jacobians saved during point
# estimation are reused instead of being re-estimated in each bootstrap draw.
fit_bootstrap_onestep <- function(
    X, y, weights, multipliers, point_fit) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  weights <- normalize_weights(weights)
  multipliers <- as.numeric(multipliers)
  if (length(multipliers) != nrow(X) || any(!is.finite(multipliers)) ||
      any(multipliers < 0)) {
    stop("invalid onestep bootstrap multipliers", call. = FALSE)
  }
  bootstrap_weights <- weights * multipliers
  active <- bootstrap_weights > 0
  active_design <- X[active, , drop = FALSE]
  if (sum(multipliers) <= ncol(X) || sum(active) < ncol(X) ||
      qr(active_design, LAPACK = FALSE)$rank < ncol(X)) {
    stop("Stata-compatible bootstrap draw has an unidentified design",
         call. = FALSE)
  }
  taus <- point_fit$taus
  p <- ncol(X)
  coefficients <- matrix(NA_real_, p, length(taus))
  traversal <- point_fit$onestep_traversal
  start <- traversal[[1L]]
  directions <- c(
    rep(1L, length(taus) - start + 1L),
    rep(-1L, start - 1L)
  )
  first_solver <- point_fit$onestep_first_solver_resolved
  fallback_taus <- numeric()
  for (position in seq_along(traversal)) {
    current <- traversal[[position]]
    tau <- taus[[current]]
    candidate <- NULL
    if (position > 1L) {
      previous <- current - directions[[position]]
      previous_beta <- coefficients[, previous]
      fitted <- drop(X %*% previous_beta)
      inverse_jacobian <- point_fit$onestep_inverse_jacobians[[current]]
      if (!is.null(inverse_jacobian)) {
        score <- colSums(
          X * as.numeric(
            bootstrap_weights * (tau - (y <= fitted))
          )
        ) / sum(bootstrap_weights)
        update <- previous_beta + drop(inverse_jacobian %*% score)
        if (all(is.finite(update))) candidate <- update
      }
    }
    if (is.null(candidate)) {
      exact <- exact_initial_qr(
        X[active, , drop = FALSE],
        y[active],
        bootstrap_weights[active],
        tau,
        first_solver
      )
      candidate <- exact$coefficients
      if (position > 1L) fallback_taus <- c(fallback_taus, tau)
    }
    coefficients[, current] <- candidate
  }
  structure(list(
    model = "qr",
    solver = "onestep",
    solver_requested = point_fit$solver_requested,
    solver_implementation = paste0(
      "qrprocess 1.1.3 rq_boot_1step-compatible saved-Jacobian bootstrap"
    ),
    solver_exact = FALSE,
    solver_process_aware = TRUE,
    taus = taus,
    coefficients = coefficients,
    conditioned_coefficients = coefficients,
    convergence_flag = rep(0L, length(taus)),
    convergence_diagnostics_available = TRUE,
    iterations = NULL,
    warnings = if (length(fallback_taus)) {
      paste0(
        "onestep bootstrap used exact ", first_solver,
        " fallback at tau=", paste(signif(fallback_taus, 5), collapse = ",")
      )
    } else {
      character()
    },
    preprocessing_sample_size = NA_integer_,
    process_fallback_taus = numeric(),
    onestep_inverse_jacobians = point_fit$onestep_inverse_jacobians,
    onestep_first_solver_resolved = first_solver,
    onestep_traversal = traversal,
    onestep_implementation_version = point_fit$onestep_implementation_version,
    stata_source_version = point_fit$stata_source_version,
    stata_source_commit = point_fit$stata_source_commit,
    precondition_requested = FALSE,
    preconditioned = FALSE,
    preconditioning_method = "none",
    preconditioning_gram_rcond = NA_real_,
    preconditioning_condition_estimate = NA_real_,
    preconditioning_fallback_reason = NA_character_,
    preconditioning_matrix = diag(p)
  ), class = c("cf_qr_fit", "cf_conditional_fit"))
}

bootstrap_qr_xy_single <- function(
    X, y, weights, multipliers, tau, point_coefficient,
    method = c("fn", "br"), Mm_factor = 3) {
  method <- match.arg(method)
  X <- as.matrix(X)
  y <- as.numeric(y)
  weights <- normalize_weights(weights)
  n <- nrow(X)
  p <- ncol(X)
  counts <- as.numeric(multipliers)
  if (length(counts) != n || any(!is.finite(counts)) || any(counts < 0)) {
    stop("invalid xy_preprocess bootstrap multipliers", call. = FALSE)
  }
  weighted_X <- X * weights
  weighted_y <- y * weights
  point_residual <- drop(weighted_y - weighted_X %*% point_coefficient)
  inverse_upper <- tryCatch(
    backsolve(chol(crossprod(weighted_X)), diag(p)),
    error = function(e) NULL
  )
  if (is.null(inverse_upper)) {
    stop("xy_preprocess point design is singular", call. = FALSE)
  }
  band <- sqrt(rowSums((weighted_X %*% inverse_upper)^2))
  band <- pmax(band, sqrt(.Machine$double.eps))
  standardized <- point_residual / band
  ordered <- sort(standardized)
  selected_data <- counts > 0L
  xb_data <- weighted_X[selected_data, , drop = FALSE]
  yb_data <- weighted_y[selected_data]
  rb <- standardized[selected_data]
  wb_data <- counts[selected_data]
  mm <- max(p + 1L, as.integer(round(sqrt(p * n))))
  outer_iteration <- 0L
  repeat {
    outer_iteration <- outer_iteration + 1L
    if (outer_iteration > 25L) {
      stop("xy_preprocess failed to find a valid residual-sign partition", call. = FALSE)
    }
    M <- Mm_factor * mm
    lower_position <- max(1L, as.integer(ceiling(n * tau - M / 2)))
    upper_position <- min(n, as.integer(floor(n * tau + M / 2)))
    kappa <- ordered[c(lower_position, upper_position)]
    lower <- c(FALSE, FALSE, rb < kappa[[1L]])
    upper <- c(FALSE, FALSE, rb > kappa[[2L]])
    xb <- rbind(matrix(0, 2L, p), xb_data)
    yb <- c(0, 0, yb_data)
    wb <- c(1, 1, wb_data)
    restart_outer <- FALSE
    repeat {
      if (any(lower)) {
        xb[1L, ] <- colSums(xb[lower, , drop = FALSE] * wb[lower])
        yb[[1L]] <- sum(yb[lower] * wb[lower])
      } else {
        lower[[1L]] <- TRUE
      }
      if (any(upper)) {
        xb[2L, ] <- colSums(xb[upper, , drop = FALSE] * wb[upper])
        yb[[2L]] <- sum(yb[upper] * wb[upper])
      } else {
        upper[[2L]] <- TRUE
      }
      keep <- !upper & !lower
      fit <- tryCatch(
        quantreg::rq.wfit(
          xb[keep, , drop = FALSE],
          yb[keep],
          weights = wb[keep],
          tau = tau,
          method = method
        ),
        error = function(e) NULL
      )
      if (is.null(fit) || any(!is.finite(fit$coefficients))) {
        mm <- min(n, 2L * mm)
        restart_outer <- TRUE
        break
      }
      coefficient <- as.numeric(fit$coefficients)
      bootstrap_residual <- drop(yb_data - xb_data %*% coefficient)
      upper_bad <- c(FALSE, FALSE, bootstrap_residual < 0) & upper
      lower_bad <- c(FALSE, FALSE, bootstrap_residual > 0) & lower
      bad_signs <- sum(upper_bad | lower_bad)
      if (bad_signs == 0L) return(coefficient)
      if (bad_signs > 0.1 * M) {
        mm <- min(n, 2L * mm)
        restart_outer <- TRUE
        break
      }
      upper <- upper & !upper_bad
      lower <- lower & !lower_bad
    }
    if (!restart_outer) break
  }
  stop("xy_preprocess terminated unexpectedly", call. = FALSE)
}

fit_bootstrap_xy_preprocess <- function(
    X, y, weights, multipliers, point_fit) {
  solver <- point_fit$solver
  method <- if (solver == "proqreg") "br" else "fn"
  transform <- point_fit$preconditioning_matrix
  conditioned_X <- as.matrix(X) %*% transform
  conditioned_coefficients <- point_fit$conditioned_coefficients
  coefficients <- vapply(seq_along(point_fit$taus), function(j) {
    bootstrap_qr_xy_single(
      conditioned_X,
      y,
      weights,
      multipliers,
      point_fit$taus[[j]],
      conditioned_coefficients[, j],
      method = method
    )
  }, numeric(ncol(X)))
  coefficients <- normalize_coefficient_matrix(
    coefficients, ncol(X), point_fit$taus
  )
  coefficients <- transform %*% coefficients
  structure(list(
    model = "qr",
    solver = solver,
    solver_requested = point_fit$solver_requested,
    solver_implementation = paste0(
      "quantreg::boot.rq.pxy-compatible preprocessing; method=", method
    ),
    solver_exact = TRUE,
    solver_process_aware = TRUE,
    taus = point_fit$taus,
    coefficients = coefficients,
    convergence_flag = rep(0L, length(point_fit$taus)),
    convergence_diagnostics_available = TRUE,
    iterations = NULL,
    warnings = character(),
    preprocessing_sample_size = as.integer(round(sqrt(nrow(X) * ncol(X)))),
    preconditioned = point_fit$preconditioned,
    preconditioning_method = point_fit$preconditioning_method,
    preconditioning_gram_rcond = point_fit$preconditioning_gram_rcond,
    preconditioning_condition_estimate = point_fit$preconditioning_condition_estimate,
    preconditioning_fallback_reason = point_fit$preconditioning_fallback_reason,
    preconditioning_matrix = transform
  ), class = c("cf_qr_fit", "cf_conditional_fit"))
}

bootstrap_checkpoint_path <- function(run_dir, rep_id) {
  file.path(run_dir, sprintf("bootstrap_rep_%04d.rds", rep_id))
}

valid_bootstrap_checkpoint <- function(path, common, rep_id) {
  if (!file.exists(path)) return(NULL)
  existing <- tryCatch(readRDS(path), error = function(e) NULL)
  required <- c(
    "signature", "data_fingerprint", "replication", "attempt", "seed",
    "effects", "elapsed_seconds", "peak_r_heap_mb", "warnings",
    "dropped_design_columns", "group0_active_rows", "group1_active_rows",
    "active_rows", "resources", "bootstrap_engine", "bootstrap_scheme",
    "draw_point_workers", "attempt_failures"
  )
  expected_effects <- c("structure", "composition", "total")
  effects_valid <- function(effects) {
    values <- as.numeric(effects)
    if (identical(common$model, "cox")) {
      return(all(is.finite(values) | (is.na(values) & !is.nan(values))))
    }
    all(is.finite(values))
  }
  nonnegative_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 0
  }
  active_counts_valid <- function(object) {
    values <- object[c(
      "group0_active_rows", "group1_active_rows", "active_rows",
      "draw_point_workers"
    )]
    typed <- all(vapply(values, function(value) {
      is.integer(value) && length(value) == 1L && !is.na(value)
    }, logical(1L)))
    typed && object$group0_active_rows >= 1L &&
      object$group0_active_rows <= common$prepared$n0 &&
      object$group1_active_rows >= 1L &&
      object$group1_active_rows <= common$prepared$n1 &&
      object$active_rows == object$group0_active_rows + object$group1_active_rows &&
      object$draw_point_workers >= 1L
  }
  resources_valid <- function(resources) {
    is.data.frame(resources) && nrow(resources) >= 1L &&
      all(c("phase", "elapsed_seconds", "peak_r_heap_mb") %in% names(resources)) &&
      is.character(resources$phase) &&
      is.numeric(resources$elapsed_seconds) &&
      is.numeric(resources$peak_r_heap_mb) &&
      all(is.finite(resources$elapsed_seconds) & resources$elapsed_seconds >= 0) &&
      all(is.finite(resources$peak_r_heap_mb) & resources$peak_r_heap_mb >= 0)
  }
  failures_valid <- function(failures) {
    if (!is.data.frame(failures)) return(FALSE)
    if (!nrow(failures)) return(TRUE)
    required_failures <- c("replication", "attempt", "seed", "error")
    all(required_failures %in% names(failures)) &&
      is.numeric(failures$replication) &&
      is.numeric(failures$attempt) &&
      is.numeric(failures$seed) &&
      is.character(failures$error) &&
      all(is.finite(failures$replication)) &&
      all(is.finite(failures$attempt)) && all(failures$attempt >= 0) &&
      all(is.finite(failures$seed)) &&
      all(nzchar(failures$error))
  }
  valid <- tryCatch({
    is.list(existing) && all(required %in% names(existing)) &&
    identical(existing$signature, common$signature) &&
    identical(existing$data_fingerprint, common$data_fingerprint) &&
    identical(existing$replication, as.integer(rep_id)) &&
    is.integer(existing$attempt) && length(existing$attempt) == 1L &&
    !is.na(existing$attempt) && existing$attempt >= 0L &&
    is.integer(existing$seed) && length(existing$seed) == 1L &&
    !is.na(existing$seed) &&
    is.matrix(existing$effects) && is.numeric(existing$effects) &&
    identical(dim(existing$effects), c(
      3L, as.integer(length(common$control$reported_quantiles))
    )) && identical(rownames(existing$effects), expected_effects) &&
    effects_valid(existing$effects) &&
    nonnegative_scalar(existing$elapsed_seconds) &&
    nonnegative_scalar(existing$peak_r_heap_mb) &&
    is.character(existing$warnings) &&
    is.character(existing$dropped_design_columns) &&
    active_counts_valid(existing) &&
    resources_valid(existing$resources) &&
    identical(existing$bootstrap_engine, common$bootstrap_engine) &&
    identical(existing$bootstrap_scheme, common$control$bootstrap_scheme) &&
    identical(existing$draw_point_workers, common$draw_point_workers) &&
    failures_valid(existing$attempt_failures)
  }, error = function(error) FALSE)
  if (!isTRUE(valid)) {
    return(NULL)
  }
  existing
}

bootstrap_attempt_seed <- function(seed, rep_id, attempt = 0L, stream = 0L) {
  modulus <- .Machine$integer.max - 1
  value <- as.numeric(seed) + as.numeric(rep_id) +
    as.numeric(attempt) * 1000003 + as.numeric(stream) * 10000019
  as.integer((value - 1) %% modulus + 1)
}

bootstrap_replication_attempt <- function(rep_id, attempt, common) {
  rep_id <- as.integer(rep_id)
  attempt <- as.integer(attempt)
  seed_used <- bootstrap_attempt_seed(common$seed, rep_id, attempt)
  set.seed(seed_used)
  prepared <- common$prepared
  multiplier0 <- bootstrap_multipliers(
    prepared$n0, common$control$bootstrap_scheme
  )
  multiplier1 <- bootstrap_multipliers(
    prepared$n1, common$control$bootstrap_scheme
  )
  active0 <- multiplier0 > 0
  active1 <- multiplier1 > 0
  boot_data <- prepared
  boot_data$X0 <- prepared$X0[active0, , drop = FALSE]
  boot_data$y0 <- prepared$y0[active0]
  boot_data$w0 <- prepared$w0[active0] * multiplier0[active0]
  boot_data$quantile_frequency0 <- if (
    common$control$bootstrap_scheme == "multiplier"
  ) rep(1, sum(active0)) else as.numeric(multiplier0[active0])
  if (!is.null(prepared$censoring0)) {
    boot_data$censoring0 <- prepared$censoring0[active0]
  }
  if (!is.null(prepared$event0)) {
    boot_data$event0 <- prepared$event0[active0]
  }
  boot_data$X1 <- prepared$X1[active1, , drop = FALSE]
  boot_data$y1 <- prepared$y1[active1]
  boot_data$w1 <- prepared$w1[active1] * multiplier1[active1]
  boot_data$quantile_frequency1 <- if (
    common$control$bootstrap_scheme == "multiplier"
  ) rep(1, sum(active1)) else as.numeric(multiplier1[active1])
  if (!is.null(prepared$censoring1)) {
    boot_data$censoring1 <- prepared$censoring1[active1]
  }
  if (!is.null(prepared$event1)) {
    boot_data$event1 <- prepared$event1[active1]
  }
  boot_data$n0 <- sum(active0)
  boot_data$n1 <- sum(active1)
  boot_data$n <- boot_data$n0 + boot_data$n1
  estimation_data <- boot_data
  if (isTRUE(common$control$legacy_weighted_quantile)) {
    # Counterfactual 1.2 divided resampled weights by the original resample
    # size. Keep n0/n1/n on boot_data as physical retained-row counts while
    # using that legacy denominator only for decomposition evaluation.
    estimation_data$n <- prepared$n
  }
  if (common$bootstrap_engine == "onestep") {
    started <- proc.time()[["elapsed"]]
    group0 <- measure_resources(function() {
      fit_bootstrap_onestep(
        prepared$X0, prepared$y0, prepared$w0, multiplier0,
        common$point$fits$group0
      )
    })
    group1 <- measure_resources(function() {
      fit_bootstrap_onestep(
        prepared$X1, prepared$y1, prepared$w1, multiplier1,
        common$point$fits$group1
      )
    })
    evaluated <- evaluate_decomposition(
      estimation_data, group0$value, group1$value, common$control
    )
    point <- list(
      effects = evaluated$effects,
      elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
      warnings = unique(c(group0$value$warnings, group1$value$warnings)),
      resources = data.frame(
        phase = c("group0_fit", "group1_fit", "marginalization"),
        elapsed_seconds = c(
          group0$elapsed_seconds,
          group1$elapsed_seconds,
          evaluated$resource$elapsed_seconds
        ),
        peak_r_heap_mb = c(
          group0$peak_r_heap_mb,
          group1$peak_r_heap_mb,
          evaluated$resource$peak_r_heap_mb
        )
      )
    )
  } else if (common$bootstrap_engine == "xy_preprocess") {
    started <- proc.time()[["elapsed"]]
    group0 <- measure_resources(function() {
      fit_bootstrap_xy_preprocess(
        prepared$X0, prepared$y0, prepared$w0, multiplier0,
        common$point$fits$group0
      )
    })
    group1 <- measure_resources(function() {
      fit_bootstrap_xy_preprocess(
        prepared$X1, prepared$y1, prepared$w1, multiplier1,
        common$point$fits$group1
      )
    })
    evaluated <- evaluate_decomposition(
      estimation_data, group0$value, group1$value, common$control
    )
    point <- list(
      effects = evaluated$effects,
      elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
      warnings = character(),
      resources = data.frame(
        phase = c("group0_fit", "group1_fit", "marginalization"),
        elapsed_seconds = c(
          group0$elapsed_seconds,
          group1$elapsed_seconds,
          evaluated$resource$elapsed_seconds
        ),
        peak_r_heap_mb = c(
          group0$peak_r_heap_mb,
          group1$peak_r_heap_mb,
          evaluated$resource$peak_r_heap_mb
        )
      )
    )
  } else {
    estimation_data <- reduce_prepared_design(estimation_data)
    point <- estimate_point_prepared(
      estimation_data,
      common$model,
      common$solver,
      common$control,
      point_workers = common$draw_point_workers,
      point_seed = bootstrap_attempt_seed(
        common$seed, rep_id, attempt, stream = 1L
      ),
      keep_fits = FALSE
    )
  }
  checkpoint <- list(
    signature = common$signature,
    data_fingerprint = common$data_fingerprint,
    replication = rep_id,
    attempt = attempt,
    seed = seed_used,
    effects = point$effects,
    elapsed_seconds = point$elapsed_seconds,
    peak_r_heap_mb = max(point$resources$peak_r_heap_mb),
    warnings = point$warnings,
    dropped_design_columns = estimation_data$dropped_design_columns,
    group0_active_rows = as.integer(boot_data$n0),
    group1_active_rows = as.integer(boot_data$n1),
    active_rows = as.integer(boot_data$n),
    resources = point$resources,
    bootstrap_engine = common$bootstrap_engine,
    bootstrap_scheme = common$control$bootstrap_scheme,
    draw_point_workers = common$draw_point_workers,
    attempt_failures = data.frame()
  )
  checkpoint
}

bootstrap_replication <- function(rep_id, common) {
  rep_id <- as.integer(rep_id)
  path <- bootstrap_checkpoint_path(common$run_dir, rep_id)
  existing <- valid_bootstrap_checkpoint(path, common, rep_id)
  if (!is.null(existing)) {
    return(list(
      status = "ok", path = path, replication = rep_id,
      attempt = existing$attempt %||% 0L, cached = TRUE,
      failures = existing$attempt_failures %||% data.frame()
    ))
  }

  failures <- list()
  for (attempt in 0:common$control$bootstrap_max_retries) {
    checkpoint <- tryCatch(
      bootstrap_replication_attempt(rep_id, attempt, common),
      error = identity
    )
    if (!inherits(checkpoint, "error")) {
      failure_table <- if (length(failures)) {
        do.call(rbind, failures)
      } else {
        data.frame()
      }
      checkpoint$attempt_failures <- failure_table
      atomic_save_rds(checkpoint, path, compress = FALSE)
      return(list(
        status = "ok", path = path, replication = rep_id,
        attempt = as.integer(attempt), cached = FALSE,
        failures = failure_table
      ))
    }
    failures[[length(failures) + 1L]] <- data.frame(
      replication = rep_id,
      attempt = as.integer(attempt),
      seed = bootstrap_attempt_seed(common$seed, rep_id, attempt),
      error = conditionMessage(checkpoint),
      stringsAsFactors = FALSE
    )
  }
  list(
    status = "error", path = NA_character_, replication = rep_id,
    attempt = common$control$bootstrap_max_retries,
    cached = FALSE,
    failures = do.call(rbind, failures)
  )
}

initialize_bootstrap_worker <- function(common) {
  .cf_bootstrap_worker_state$common <- common
  invisible(NULL)
}

bootstrap_worker_replication <- function(rep_id) {
  if (is.null(.cf_bootstrap_worker_state$common)) {
    stop("bootstrap worker state was not initialized", call. = FALSE)
  }
  bootstrap_replication(rep_id, .cf_bootstrap_worker_state$common)
}

report_bootstrap_progress <- function(done, total, started) {
  elapsed <- proc.time()[["elapsed"]] - started
  rate <- if (done > 0L) elapsed / done else NA_real_
  eta <- if (is.finite(rate)) rate * (total - done) else NA_real_
  message(sprintf(
    "Bootstrap: %d/%d complete (%.1f%%), elapsed %.1f min, ETA %.1f min",
    done, total, 100 * done / total, elapsed / 60, eta / 60
  ))
}

bind_bootstrap_failures <- function(outcomes) {
  tables <- lapply(outcomes, `[[`, "failures")
  tables <- tables[vapply(tables, nrow, integer(1L)) > 0L]
  if (!length(tables)) {
    return(data.frame(
      replication = integer(), attempt = integer(), seed = integer(),
      error = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, tables)
}

run_cf_bootstrap <- function(
    prepared, point, model, solver, control, reps, workers,
    draw_point_workers = 1L, checkpoint_dir, seed) {
  workers <- min(workers, reps)
  draw_point_workers <- min(assert_scalar_integer(
    draw_point_workers, "draw_point_workers", 1L
  ), 2L)
  data_fingerprint <- bootstrap_data_fingerprint(prepared)
  bootstrap_engine <- resolve_qr_bootstrap_engine(
    control$qr_bootstrap_engine, model, point, control$bootstrap_scheme
  )
  if (workers > 1L || bootstrap_engine != "standard") {
    draw_point_workers <- 1L
  }
  signature <- bootstrap_signature(
    prepared, model, solver, control, seed, bootstrap_engine,
    draw_point_workers, data_fingerprint, point
  )
  run_dir <- file.path(checkpoint_dir, paste0("run_", signature))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  common <- list(
    prepared = prepared,
    point = point,
    model = model,
    solver = solver,
    control = control,
    seed = seed,
    signature = signature,
    data_fingerprint = data_fingerprint,
    bootstrap_engine = bootstrap_engine,
    draw_point_workers = draw_point_workers,
    run_dir = run_dir
  )
  rep_ids <- seq_len(reps)
  progress_started <- proc.time()[["elapsed"]]
  outcomes <- if (workers == 1L) {
    completed <- vector("list", reps)
    for (i in seq_along(rep_ids)) {
      completed[[i]] <- bootstrap_replication(rep_ids[[i]], common)
      if (isTRUE(control$bootstrap_progress)) {
        report_bootstrap_progress(i, reps, progress_started)
      }
    }
    completed
  } else {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    package_worker_init(cluster)
    parallel::clusterCall(cluster, function(worker_common) {
      initialize <- get(
        "initialize_bootstrap_worker",
        envir = asNamespace("scalableCounterfactual")
      )
      initialize(worker_common)
    }, common)
    batches <- split(rep_ids, ceiling(seq_along(rep_ids) / workers))
    completed <- vector("list", reps)
    offset <- 0L
    for (batch in batches) {
      batch_results <- parallel::parLapplyLB(cluster, batch, function(rep_id) {
        run_replication <- get(
          "bootstrap_worker_replication",
          envir = asNamespace("scalableCounterfactual")
        )
        run_replication(rep_id)
      })
      completed[seq.int(offset + 1L, offset + length(batch_results))] <-
        batch_results
      offset <- offset + length(batch_results)
      if (isTRUE(control$bootstrap_progress)) {
        report_bootstrap_progress(offset, reps, progress_started)
      }
    }
    completed
  }
  failure_table <- bind_bootstrap_failures(outcomes)
  if (nrow(failure_table)) {
    data.table::fwrite(
      failure_table, file.path(run_dir, "bootstrap_failures.csv")
    )
  }
  final_failures <- vapply(outcomes, function(x) x$status != "ok", logical(1L))
  if (any(final_failures)) {
    stop(
      sum(final_failures), " bootstrap replication(s) failed after ",
      control$bootstrap_max_retries + 1L, " attempts; diagnostics: ",
      file.path(run_dir, "bootstrap_failures.csv"),
      call. = FALSE
    )
  }
  paths <- vapply(outcomes, `[[`, character(1L), "path")
  checkpoints <- lapply(seq_along(paths), function(index) {
    valid_bootstrap_checkpoint(paths[[index]], common, rep_ids[[index]])
  })
  if (any(vapply(checkpoints, is.null, logical(1L)))) {
    stop("bootstrap checkpoint failed final integrity validation", call. = FALSE)
  }
  effect_names <- c("structure", "composition", "total")
  effect_draws <- stats::setNames(lapply(effect_names, function(effect) {
    do.call(rbind, lapply(checkpoints, function(x) x$effects[effect, ]))
  }), effect_names)
  list(
    signature = signature,
    data_fingerprint = data_fingerprint,
    engine = bootstrap_engine,
    scheme = control$bootstrap_scheme,
    checkpoint_dir = run_dir,
    replications = reps,
    draw_point_workers = draw_point_workers,
    retries = sum(vapply(checkpoints, `[[`, integer(1L), "attempt")),
    failures = failure_table,
    effects = effect_draws,
    warnings = unique(unlist(lapply(checkpoints, `[[`, "warnings"))),
    dropped_design_columns = unique(unlist(lapply(
      checkpoints,
      `[[`,
      "dropped_design_columns"
    ))),
    resources = data.frame(
      replication = vapply(checkpoints, `[[`, integer(1L), "replication"),
      attempt = vapply(checkpoints, `[[`, integer(1L), "attempt"),
      seed = vapply(checkpoints, `[[`, integer(1L), "seed"),
      cached = vapply(outcomes, `[[`, logical(1L), "cached"),
      elapsed_seconds = vapply(checkpoints, `[[`, numeric(1L), "elapsed_seconds"),
      peak_r_heap_mb = vapply(checkpoints, `[[`, numeric(1L), "peak_r_heap_mb")
    )
  )
}
