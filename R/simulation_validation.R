simulation_scenarios <- function() {
  list(
    location_shift = list(
      mu0 = 0, mu1 = 0, a0 = 0, a1 = 0.20,
      b0 = 0.50, b1 = 0.50, sigma = 1, x_sd = 1
    ),
    composition_shift = list(
      mu0 = 0, mu1 = 0.50, a0 = 0, a1 = 0,
      b0 = 0.60, b1 = 0.60, sigma = 1, x_sd = 1
    ),
    combined_scale_shift = list(
      mu0 = 0, mu1 = 0.50, a0 = 0, a1 = 0.15,
      b0 = 0.50, b1 = 0.80, sigma = 0.80, x_sd = 1
    )
  )
}

gaussian_true_decomposition <- function(parameters, quantiles) {
  z <- stats::qnorm(quantiles)
  scale0 <- sqrt((parameters$b0 * parameters$x_sd)^2 + parameters$sigma^2)
  scale1 <- sqrt((parameters$b1 * parameters$x_sd)^2 + parameters$sigma^2)
  q00 <- parameters$a0 + parameters$b0 * parameters$mu0 + z * scale0
  q01 <- parameters$a0 + parameters$b0 * parameters$mu1 + z * scale0
  q11 <- parameters$a1 + parameters$b1 * parameters$mu1 + z * scale1
  rbind(
    structure = q11 - q01,
    composition = q01 - q00,
    total = q11 - q00
  )
}

simulation_checkpoint_signature <- function(task) {
  object_md5(list(
    schema_version = 1L,
    package_version = as.character(
      utils::packageVersion("scalableCounterfactual")
    ),
    scenario = task$scenario,
    parameters = task$parameters,
    replication = task$replication,
    n_per_group = task$n_per_group,
    solver = task$solver,
    nreg = task$nreg,
    trimming = task$trimming,
    reported_quantiles = task$reported_quantiles,
    weighted = task$weighted,
    bootstrap_reps = task$bootstrap_reps,
    bootstrap_scheme = task$bootstrap_scheme,
    master_seed = task$master_seed,
    max_retries = task$max_retries,
    task_timeout_seconds = task$task_timeout_seconds,
    extension_registry = extension_registry_fingerprint(list(
      qr = intersect(task$solver, names(registry_snapshot()$qr)),
      linear = character(), distribution = character()
    ))
  ))
}

run_simulation_attempt <- function(task) {
  if (is.finite(task$task_timeout_seconds)) {
    setTimeLimit(elapsed = task$task_timeout_seconds, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  }
  run_simulation_replication(task)
}

run_simulation_task <- function(task) {
  signature <- simulation_checkpoint_signature(task)
  path <- task$checkpoint_path
  if (!is.null(path) && file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(existing) && identical(existing$signature, signature)) {
      existing$cached <- TRUE
      return(existing)
    }
  }
  failures <- list()
  for (attempt in 0:task$max_retries) {
    attempt_task <- task
    attempt_task$seed <- bootstrap_attempt_seed(
      task$master_seed, task$task_id, attempt, stream = 3L
    )
    result <- tryCatch(
      run_simulation_attempt(attempt_task),
      error = function(error) list(
        status = "error", error = conditionMessage(error), elapsed = NA_real_,
        rows = NULL, identity_residual = NA_real_
      )
    )
    result$attempt <- as.integer(attempt)
    result$seed <- attempt_task$seed
    result$cached <- FALSE
    if (result$status == "ok") {
      result$signature <- signature
      result$attempt_failures <- if (length(failures)) {
        do.call(rbind, failures)
      } else data.frame()
      if (!is.null(path)) atomic_save_rds(result, path, compress = FALSE)
      return(result)
    }
    failures[[length(failures) + 1L]] <- data.frame(
      scenario = task$scenario,
      replication = task$replication,
      attempt = as.integer(attempt),
      seed = attempt_task$seed,
      error = result$error,
      stringsAsFactors = FALSE
    )
  }
  result$signature <- signature
  result$attempt_failures <- do.call(rbind, failures)
  if (!is.null(path)) atomic_save_rds(result, path, compress = FALSE)
  result
}

run_simulation_replication <- function(task) {
  set.seed(task$seed)
  parameters <- task$parameters
  n <- task$n_per_group
  group <- rep(0:1, each = n)
  x <- c(
    stats::rnorm(n, parameters$mu0, parameters$x_sd),
    stats::rnorm(n, parameters$mu1, parameters$x_sd)
  )
  intercept <- ifelse(group == 0L, parameters$a0, parameters$a1)
  slope <- ifelse(group == 0L, parameters$b0, parameters$b1)
  y <- intercept + slope * x + stats::rnorm(2L * n, 0, parameters$sigma)
  weights <- if (isTRUE(task$weighted)) {
    stats::runif(2L * n, 0.5, 1.5)
  } else {
    rep(1, 2L * n)
  }
  data <- data.frame(y = y, x = x, group = group, weight = weights)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(counterfactual_decompose(
    y ~ x,
    data = data,
    group = "group",
    weights = "weight",
    model = "qr",
    solver = task$solver,
    control = cf_control(
      nreg = task$nreg,
      trimming = task$trimming,
      reported_quantiles = task$reported_quantiles,
      bootstrap_scheme = task$bootstrap_scheme,
      bootstrap_progress = FALSE,
      crossing_diagnostics = FALSE,
      marginal_method = "matrix"
    ),
    bootstrap_reps = task$bootstrap_reps,
    point_workers = 1L,
    bootstrap_workers = 1L,
    seed = task$seed
  ), error = identity)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(fit, "error")) {
    return(list(
      status = "error", error = conditionMessage(fit), elapsed = elapsed,
      rows = NULL, identity_residual = NA_real_
    ))
  }
  truth <- gaussian_true_decomposition(
    parameters, task$reported_quantiles
  )
  rows <- do.call(rbind, lapply(
    c("structure", "composition", "total"), function(effect) {
      estimated <- fit$results[fit$results$effect == effect, , drop = FALSE]
      true_value <- as.numeric(truth[effect, ])
      data.frame(
        scenario = task$scenario,
        replication = task$replication,
        quantile = task$reported_quantiles,
        effect = effect,
        truth = true_value,
        estimate = estimated$estimate,
        error = estimated$estimate - true_value,
        squared_error = (estimated$estimate - true_value)^2,
        pointwise_covered = if (task$bootstrap_reps > 1L) {
          estimated$pointwise_lower <= true_value &
            estimated$pointwise_upper >= true_value
        } else NA,
        uniform_curve_covered = if (task$bootstrap_reps > 1L) {
          rep(all(
            estimated$uniform_lower <= true_value &
              estimated$uniform_upper >= true_value
          ), length(true_value))
        } else NA,
        stringsAsFactors = FALSE
      )
    }
  ))
  list(
    status = "ok", error = NA_character_, elapsed = unname(elapsed),
    rows = rows,
    identity_residual = max(abs(fit$point$identity_residual))
  )
}

#' Monte Carlo validation of counterfactual QR decomposition
#'
#' Simulates Gaussian linear designs with analytically known structure,
#' composition, and total quantile effects. The built-in scenarios separately
#' exercise pure structure, pure composition, and combined location/scale
#' changes. Optional nested bootstrap runs also report pointwise and uniform
#' empirical coverage.
#'
#' @param replications Monte Carlo replications per scenario.
#' @param n_per_group Observations generated in each group.
#' @param scenarios Any subset of `location_shift`, `composition_shift`, and
#'   `combined_scale_shift`.
#' @param solver Registered QR solver.
#' @param nreg Conditional-quantile grid size.
#' @param trimming Conditional QR trimming.
#' @param reported_quantiles Marginal quantiles evaluated.
#' @param weighted Generate independent positive sampling weights.
#' @param bootstrap_reps Bootstrap replications within each Monte Carlo draw.
#' @param bootstrap_scheme Bootstrap sampling scheme.
#' @param workers Parallel Monte Carlo workers. Nested bootstrap workers remain
#'   fixed at one.
#' @param seed Master simulation seed.
#' @param progress Print completed simulation tasks.
#' @param checkpoint_dir Optional directory for replication-level checkpoints.
#'   Completed tasks are reused when the full simulation signature matches.
#' @param max_retries Deterministic replacement attempts after a failed task.
#' @param task_timeout_seconds Best-effort elapsed-time limit for one task, or
#'   `Inf`. Compiled solvers that do not check interrupts may exceed this limit.
#' @return An object of class `cf_simulation_validation` containing raw errors,
#'   aggregate performance, failures, and run settings.
#' @export
simulate_counterfactual_validation <- function(
    replications = 100L,
    n_per_group = 1000L,
    scenarios = names(simulation_scenarios()),
    solver = "fn",
    nreg = 49L,
    trimming = 0.02,
    reported_quantiles = c(0.1, 0.5, 0.9),
    weighted = FALSE,
    bootstrap_reps = 0L,
    bootstrap_scheme = "counterfactual",
    workers = 1L,
    seed = 20260809L,
    progress = interactive(),
    checkpoint_dir = NULL,
    max_retries = 1L,
    task_timeout_seconds = Inf) {
  replications <- assert_scalar_integer(replications, "replications", 1L)
  n_per_group <- assert_scalar_integer(n_per_group, "n_per_group", 20L)
  nreg <- assert_scalar_integer(nreg, "nreg", 3L)
  bootstrap_reps <- assert_scalar_integer(
    bootstrap_reps, "bootstrap_reps", 0L
  )
  workers <- assert_scalar_integer(workers, "workers", 1L)
  seed <- assert_scalar_integer(seed, "seed", 1L)
  max_retries <- assert_scalar_integer(max_retries, "max_retries", 0L)
  if (length(task_timeout_seconds) != 1L || is.na(task_timeout_seconds) ||
      task_timeout_seconds <= 0) {
    stop("task_timeout_seconds must be positive or Inf", call. = FALSE)
  }
  weighted <- assert_scalar_logical(weighted, "weighted")
  progress <- assert_scalar_logical(progress, "progress")
  solver <- match.arg(solver, supported_qr_solvers())
  bootstrap_scheme <- match.arg(
    bootstrap_scheme, c("counterfactual", "empirical", "multiplier")
  )
  trimming <- assert_probability(trimming, "trimming", open = FALSE)
  if (trimming >= 0.5) stop("trimming must be below 0.5", call. = FALSE)
  reported_quantiles <- sort(unique(as.numeric(reported_quantiles)))
  if (!length(reported_quantiles) || any(!is.finite(reported_quantiles)) ||
      any(reported_quantiles <= 0 | reported_quantiles >= 1)) {
    stop("reported_quantiles must lie strictly between 0 and 1", call. = FALSE)
  }
  definitions <- simulation_scenarios()
  scenarios <- unique(as.character(scenarios))
  invalid <- setdiff(scenarios, names(definitions))
  if (!length(scenarios) || length(invalid)) {
    stop("unsupported simulation scenario(s): ", paste(invalid, collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(checkpoint_dir)) {
    if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L ||
        is.na(checkpoint_dir) || !nzchar(checkpoint_dir)) {
      stop("checkpoint_dir must be NULL or one nonempty path", call. = FALSE)
    }
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  tasks <- list()
  task_id <- 0L
  for (scenario in scenarios) {
    for (replication in seq_len(replications)) {
      task_id <- task_id + 1L
      tasks[[task_id]] <- list(
        scenario = scenario,
        parameters = definitions[[scenario]],
        replication = replication,
        n_per_group = n_per_group,
        solver = solver,
        nreg = nreg,
        trimming = trimming,
        reported_quantiles = reported_quantiles,
        weighted = weighted,
        bootstrap_reps = bootstrap_reps,
        bootstrap_scheme = bootstrap_scheme,
        master_seed = seed,
        task_id = task_id,
        max_retries = max_retries,
        task_timeout_seconds = as.numeric(task_timeout_seconds),
        checkpoint_path = if (is.null(checkpoint_dir)) NULL else file.path(
          checkpoint_dir,
          sprintf("%s_rep_%05d.rds", scenario, replication)
        )
      )
    }
  }
  started <- proc.time()[["elapsed"]]
  if (workers == 1L) {
    results <- lapply(seq_along(tasks), function(i) {
      result <- run_simulation_task(tasks[[i]])
      if (progress) message("Simulation ", i, "/", length(tasks))
      result
    })
  } else {
    cluster <- parallel::makeCluster(min(workers, length(tasks)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    package_worker_init(cluster)
    batches <- split(seq_along(tasks), ceiling(seq_along(tasks) /
      min(workers, length(tasks))))
    results <- vector("list", length(tasks))
    completed <- 0L
    for (indices in batches) {
      batch_results <- parallel::parLapplyLB(cluster, tasks[indices], function(task) {
        runner <- get(
          "run_simulation_task",
          envir = asNamespace("scalableCounterfactual")
        )
        runner(task)
      })
      results[indices] <- batch_results
      completed <- completed + length(indices)
      if (progress) message("Simulation ", completed, "/", length(tasks))
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  successful <- vapply(results, function(x) x$status == "ok", logical(1L))
  raw <- data.table::rbindlist(
    lapply(results[successful], `[[`, "rows"), fill = TRUE
  )
  failure_tables <- lapply(results, `[[`, "attempt_failures")
  failure_tables <- failure_tables[vapply(failure_tables, nrow, integer(1L)) > 0L]
  failures <- if (length(failure_tables)) do.call(rbind, failure_tables) else {
    data.frame(
      scenario = character(), replication = integer(), attempt = integer(),
      seed = integer(), error = character(), stringsAsFactors = FALSE
    )
  }
  summary <- if (nrow(raw)) raw[, .(
    successful_replications = .N,
    truth = mean(truth),
    mean_estimate = mean(estimate),
    bias = mean(error),
    absolute_bias = abs(mean(error)),
    monte_carlo_sd = stats::sd(estimate),
    rmse = sqrt(mean(squared_error)),
    pointwise_coverage = if (all(is.na(pointwise_covered))) {
      NA_real_
    } else mean(pointwise_covered, na.rm = TRUE)
  ), by = .(scenario, effect, quantile)] else data.table::data.table()
  curve_coverage <- if (nrow(raw)) {
    curve_rows <- unique(raw[, .(
      scenario, replication, effect, uniform_curve_covered
    )])
    curve_rows[, .(
      successful_replications = .N,
      uniform_coverage = if (all(is.na(uniform_curve_covered))) {
        NA_real_
      } else mean(uniform_curve_covered, na.rm = TRUE)
    ), by = .(scenario, effect)]
  } else data.table::data.table()
  resources <- data.frame(
    tasks = length(tasks),
    successful_tasks = sum(successful),
    failed_tasks = sum(!successful),
    elapsed_seconds = unname(elapsed),
    mean_task_seconds = mean(vapply(results, `[[`, numeric(1L), "elapsed")),
    cached_tasks = sum(vapply(results, `[[`, logical(1L), "cached")),
    retry_attempts = sum(vapply(results, `[[`, integer(1L), "attempt")),
    maximum_identity_residual = if (any(successful)) {
      max(vapply(
        results[successful], `[[`, numeric(1L), "identity_residual"
      ))
    } else NA_real_,
    stringsAsFactors = FALSE
  )
  structure(list(
    summary = as.data.frame(summary),
    curve_coverage = as.data.frame(curve_coverage),
    raw = as.data.frame(raw),
    failures = failures,
    resources = resources,
    settings = list(
      replications = replications,
      n_per_group = n_per_group,
      scenarios = scenarios,
      solver = solver,
      nreg = nreg,
      trimming = trimming,
      reported_quantiles = reported_quantiles,
      weighted = weighted,
      bootstrap_reps = bootstrap_reps,
      bootstrap_scheme = bootstrap_scheme,
      workers = workers,
      seed = seed,
      checkpoint_dir = if (is.null(checkpoint_dir)) NA_character_ else
        normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE),
      max_retries = max_retries,
      task_timeout_seconds = as.numeric(task_timeout_seconds),
      package_version = as.character(
        utils::packageVersion("scalableCounterfactual")
      )
    )
  ), class = "cf_simulation_validation")
}

#' @export
print.cf_simulation_validation <- function(x, ...) {
  cat("Counterfactual decomposition simulation validation\n")
  cat("  scenarios:", paste(x$settings$scenarios, collapse = ", "), "\n")
  cat("  replications per scenario:", x$settings$replications, "\n")
  cat("  observations per group:", x$settings$n_per_group, "\n")
  cat("  solver:", x$settings$solver, "\n")
  cat("  failed tasks:", x$resources$failed_tasks, "\n")
  cat("  cached tasks:", x$resources$cached_tasks, "\n")
  cat("  retry attempts:", x$resources$retry_attempts, "\n")
  cat("  elapsed seconds:", round(x$resources$elapsed_seconds, 2), "\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Write simulation-validation outputs
#'
#' @param object A `cf_simulation_validation` object.
#' @param output_dir Destination directory.
#' @return Invisibly returns `output_dir`.
#' @export
write_simulation_validation <- function(object, output_dir) {
  if (!inherits(object, "cf_simulation_validation")) {
    stop("object must be a cf_simulation_validation result", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(object$summary, file.path(output_dir, "summary.csv"))
  data.table::fwrite(object$raw, file.path(output_dir, "raw.csv"))
  data.table::fwrite(object$resources, file.path(output_dir, "resources.csv"))
  data.table::fwrite(
    object$curve_coverage, file.path(output_dir, "curve_coverage.csv")
  )
  failures_path <- file.path(output_dir, "failures.csv")
  if (file.exists(failures_path)) unlink(failures_path)
  if (nrow(object$failures)) {
    data.table::fwrite(object$failures, failures_path)
  }
  saveRDS(object, file.path(output_dir, "simulation_validation.rds"),
          compress = "xz")
  invisible(output_dir)
}
