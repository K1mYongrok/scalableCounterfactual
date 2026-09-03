isolated_qr_benchmark_task <- function(task_path, result_path) {
  task <- readRDS(task_path)
  restore_registry_snapshot(task$registry_snapshot)
  result <- tryCatch({
    validate_execution_parallelism(
      "qr", task$solver, task$control, task$point_workers
    )
    measured <- measure_resources(function() {
      estimate_point_prepared(
        task$prepared,
        model = "qr",
        solver = task$solver,
        control = task$control,
        point_workers = task$point_workers,
        point_seed = task$seed,
        keep_fits = TRUE
      )
    })
    point <- measured$value
    fits <- point$fits
    list(
      status = "ok",
      error = NA_character_,
      elapsed_seconds = point$elapsed_seconds,
      measured_wrapper_seconds = measured$elapsed_seconds,
      peak_r_heap_mb = max(
        measured$peak_r_heap_mb,
        max(point$resources$peak_r_heap_mb, na.rm = TRUE)
      ),
      effects = point$effects,
      group0_coefficients = fits$group0$coefficients,
      group1_coefficients = fits$group1$coefficients,
      resolved_group0_solver = fits$group0$solver,
      resolved_group1_solver = fits$group1$solver,
      group0_exact = fits$group0$solver_exact,
      group1_exact = fits$group1$solver_exact,
      warning_count = length(unique(c(fits$group0$warnings, fits$group1$warnings))),
      warnings = paste(unique(c(fits$group0$warnings, fits$group1$warnings)),
                       collapse = " | ")
    )
  }, error = function(e) {
    list(status = "error", error = conditionMessage(e))
  })
  atomic_save_rds(result, result_path, compress = FALSE)
  invisible(result_path)
}

run_isolated_qr_process <- function(
    task, poll_interval_ms = 25L, timeout_seconds = Inf) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop(
      "benchmark_qr_scaling() requires the suggested processx package",
      call. = FALSE
    )
  }
  script <- system.file(
    "scripts", "cfbenchmark_worker.R", package = "scalableCounterfactual"
  )
  if (!nzchar(script)) {
    stop("installed benchmark worker script was not found", call. = FALSE)
  }
  directory <- tempfile("scalablecf_isolated_")
  dir.create(directory, recursive = TRUE)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  task_path <- file.path(directory, "task.rds")
  result_path <- file.path(directory, "result.rds")
  stdout_path <- file.path(directory, "stdout.txt")
  stderr_path <- file.path(directory, "stderr.txt")
  task$registry_snapshot <- task$registry_snapshot %||% registry_snapshot()
  saveRDS(task, task_path, compress = FALSE)
  command <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") {
    "Rscript.exe"
  } else {
    "Rscript"
  })
  process <- processx::process$new(
    command,
    c("--vanilla", script, "--task", task_path, "--result", result_path),
    stdout = stdout_path,
    stderr = stderr_path,
    env = c(
      Sys.getenv(),
      R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)
    ),
    cleanup_tree = TRUE,
    windows_hide_window = TRUE
  )
  started <- proc.time()[["elapsed"]]
  peak_rss <- NA_real_
  repeat {
    alive <- process$is_alive()
    if (alive) {
      memory <- tryCatch(process$get_memory_info(), error = function(e) NULL)
      if (!is.null(memory) && "rss" %in% names(memory)) {
        peak_rss <- max(c(peak_rss, as.numeric(memory[["rss"]])), na.rm = TRUE)
      }
    }
    elapsed <- proc.time()[["elapsed"]] - started
    if (is.finite(timeout_seconds) && elapsed > timeout_seconds && alive) {
      process$kill_tree()
      process$wait(timeout = 5000L)
      return(list(
        status = "timeout", error = "isolated benchmark timed out",
        process_elapsed_seconds = unname(elapsed),
        peak_process_rss_mb = peak_rss / 1024^2,
        process_exit_status = process$get_exit_status(),
        stdout = if (file.exists(stdout_path)) paste(readLines(stdout_path), collapse = "\n") else "",
        stderr = if (file.exists(stderr_path)) paste(readLines(stderr_path), collapse = "\n") else ""
      ))
    }
    if (!alive) break
    Sys.sleep(poll_interval_ms / 1000)
  }
  process$wait()
  elapsed <- proc.time()[["elapsed"]] - started
  exit_status <- process$get_exit_status()
  stdout <- if (file.exists(stdout_path)) paste(readLines(stdout_path), collapse = "\n") else ""
  stderr <- if (file.exists(stderr_path)) paste(readLines(stderr_path), collapse = "\n") else ""
  if (!identical(exit_status, 0L) || !file.exists(result_path)) {
    return(list(
      status = "error",
      error = if (nzchar(stderr)) stderr else paste("worker exit status", exit_status),
      process_elapsed_seconds = unname(elapsed),
      peak_process_rss_mb = peak_rss / 1024^2,
      process_exit_status = exit_status,
      stdout = stdout,
      stderr = stderr
    ))
  }
  result <- readRDS(result_path)
  result$process_elapsed_seconds <- unname(elapsed)
  result$peak_process_rss_mb <- peak_rss / 1024^2
  result$process_exit_status <- exit_status
  result$stdout <- stdout
  result$stderr <- stderr
  result
}

summarize_scaling_benchmark <- function(raw, reference_solver) {
  finite_median <- function(values) {
    values <- values[is.finite(values)]
    if (length(values)) stats::median(values) else NA_real_
  }
  finite_max <- function(values) {
    values <- values[is.finite(values)]
    if (length(values)) max(values) else NA_real_
  }
  keys <- unique(raw[c("sample_n", "solver")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    selected <- raw$sample_n == keys$sample_n[[i]] &
      raw$solver == keys$solver[[i]]
    all_rows <- raw[selected, , drop = FALSE]
    ok <- all_rows[all_rows$status == "ok", , drop = FALSE]
    data.frame(
      sample_n = keys$sample_n[[i]],
      solver = keys$solver[[i]],
      repetitions = nrow(ok),
      failed_repetitions = nrow(all_rows) - nrow(ok),
      median_seconds = if (nrow(ok)) stats::median(ok$elapsed_seconds) else NA_real_,
      median_process_seconds = if (nrow(ok)) {
        stats::median(ok$process_elapsed_seconds)
      } else NA_real_,
      median_peak_process_rss_mb = finite_median(ok$peak_process_rss_mb),
      max_peak_process_rss_mb = finite_max(ok$peak_process_rss_mb),
      max_observed_process_rss_mb = finite_max(all_rows$peak_process_rss_mb),
      max_abs_effect_difference = if (nrow(ok) &&
          any(is.finite(ok$max_abs_effect_difference))) {
        max(ok$max_abs_effect_difference, na.rm = TRUE)
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  reference_lookup <- stats::setNames(
    summary$median_seconds[summary$solver == reference_solver],
    summary$sample_n[summary$solver == reference_solver]
  )
  summary$reference_solver <- reference_solver
  summary$solver_speed_relative_to_reference <- unname(
    reference_lookup[as.character(summary$sample_n)] / summary$median_seconds
  )
  summary$reference_speedup_over_solver <- unname(
    summary$median_seconds / reference_lookup[as.character(summary$sample_n)]
  )
  smallest <- min(summary$sample_n)
  baseline_time <- stats::setNames(
    summary$median_seconds[summary$sample_n == smallest],
    summary$solver[summary$sample_n == smallest]
  )
  baseline_rss <- stats::setNames(
    summary$median_peak_process_rss_mb[summary$sample_n == smallest],
    summary$solver[summary$sample_n == smallest]
  )
  summary$time_growth_vs_smallest <- summary$median_seconds /
    unname(baseline_time[summary$solver])
  summary$rss_growth_vs_smallest <- summary$median_peak_process_rss_mb /
    unname(baseline_rss[summary$solver])
  summary$time_scaling_exponent <- ifelse(
    summary$sample_n == smallest | !is.finite(summary$time_growth_vs_smallest),
    NA_real_,
    log(summary$time_growth_vs_smallest) /
      log(summary$sample_n / smallest)
  )
  summary <- summary[order(summary$sample_n, summary$median_seconds), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

valid_scaling_checkpoint <- function(
    checkpoint, signature, sample_sizes, solvers, repetitions,
    reported_quantiles, conditional_quantiles, design_columns,
    point_workers_requested, point_workers) {
  tryCatch({
  if (!is.list(checkpoint) || !identical(checkpoint$signature, signature) ||
      !is.list(checkpoint$rows) || !is.list(checkpoint$fits)) {
    return(FALSE)
  }
  fit_keys <- names(checkpoint$fits)
  if (is.null(fit_keys)) fit_keys <- character()
  if (any(!nzchar(fit_keys)) || anyDuplicated(fit_keys)) return(FALSE)
  if (!length(checkpoint$rows) && !length(checkpoint$fits)) return(TRUE)
  if (length(checkpoint$rows) != length(checkpoint$fits)) return(FALSE)
  required_row_fields <- c(
    "sample_n", "solver", "repetition", "run_order", "status",
    "elapsed_seconds", "process_elapsed_seconds", "peak_r_heap_mb",
    "peak_process_rss_mb", "point_workers_requested", "point_workers"
  )
  scalar_number <- function(value, finite = TRUE, minimum = -Inf) {
    valid <- is.numeric(value) && length(value) == 1L && !is.na(value)
    if (finite) valid <- valid && is.finite(value)
    valid && value >= minimum
  }
  missing_or_nonnegative <- function(value) {
    is.numeric(value) && length(value) == 1L &&
      (is.na(value) || (is.finite(value) && value >= 0))
  }
  integer_scalar <- function(value, minimum = 1L) {
    scalar_number(value, minimum = minimum) && value == floor(value)
  }
  row_keys <- vapply(checkpoint$rows, function(row) {
    if (!is.data.frame(row) || nrow(row) != 1L ||
        !all(required_row_fields %in% names(row)) ||
        !is.numeric(row$sample_n) || !is.numeric(row$repetition) ||
        !is.character(row$solver) || !is.character(row$status) ||
        length(row$sample_n) != 1L || length(row$repetition) != 1L ||
        length(row$solver) != 1L || length(row$status) != 1L ||
        !row$sample_n %in% sample_sizes || !row$solver %in% solvers ||
        !row$repetition %in% seq_len(repetitions) ||
        !row$status %in% c("ok", "error", "timeout") ||
        !integer_scalar(row$sample_n) ||
        !integer_scalar(row$repetition) ||
        !integer_scalar(row$run_order) || row$run_order > length(solvers) ||
        !integer_scalar(row$point_workers_requested) ||
        row$point_workers_requested != point_workers_requested ||
        !integer_scalar(row$point_workers) ||
        row$point_workers != point_workers ||
        !missing_or_nonnegative(row$peak_process_rss_mb) ||
        !scalar_number(row$process_elapsed_seconds, minimum = 0) ||
        (row$status == "ok" && (
          !scalar_number(row$elapsed_seconds, minimum = 0) ||
          !scalar_number(row$peak_r_heap_mb, minimum = 0)
        ))) {
      return(NA_character_)
    }
    paste(row$sample_n, row$repetition, row$solver, sep = "::")
  }, character(1L))
  if (anyNA(row_keys) || anyDuplicated(row_keys) ||
      !setequal(row_keys, fit_keys)) return(FALSE)
  all(vapply(fit_keys, function(key) {
    fit <- checkpoint$fits[[key]]
    row <- checkpoint$rows[[match(key, row_keys)]]
    if (!is.list(fit) || !is.character(fit$status) ||
        length(fit$status) != 1L || !identical(fit$status, row$status)) {
      return(FALSE)
    }
    if (fit$status != "ok") {
      return(is.character(fit$error) && length(fit$error) == 1L &&
        !is.na(fit$error) && nzchar(fit$error))
    }
    is.matrix(fit$effects) && is.numeric(fit$effects) &&
      identical(dim(fit$effects), c(3L, length(reported_quantiles))) &&
      all(is.finite(fit$effects)) &&
      is.matrix(fit$group0_coefficients) &&
      is.numeric(fit$group0_coefficients) &&
      identical(dim(fit$group0_coefficients), c(
        as.integer(design_columns), length(conditional_quantiles)
      )) && all(is.finite(fit$group0_coefficients)) &&
      is.matrix(fit$group1_coefficients) &&
      is.numeric(fit$group1_coefficients) &&
      identical(dim(fit$group1_coefficients), c(
        as.integer(design_columns), length(conditional_quantiles)
      )) && all(is.finite(fit$group1_coefficients))
  }, logical(1L)))
  }, error = function(error) FALSE)
}

#' Isolated-process QR scaling benchmark
#'
#' Runs each solver and repetition in a fresh R process. In addition to elapsed
#' time and numerical differences, it samples the operating system resident set
#' size (RSS) of the worker process. This avoids contamination from memory left
#' behind by solvers run earlier in the same R session.
#'
#' @inheritParams benchmark_qr_solvers
#' @param sample_sizes One or more stratified sample sizes.
#' @param repetitions Recorded repetitions per solver and sample size.
#' @param warmup Discarded warmup repetitions per solver and sample size.
#' @param randomize_order Randomize solver order within each size/repetition.
#' @param rss_poll_interval_ms RSS sampling interval in milliseconds.
#' @param timeout_seconds Per-process timeout; `Inf` disables it.
#' @param checkpoint_path Optional RDS checkpoint updated after every recorded
#'   solver run. An interrupted benchmark resumes completed runs by default.
#' @param resume Reuse a compatible `checkpoint_path` when it exists.
#' @return A list containing repeat-level `raw` results and a `summary`.
#' @export
benchmark_qr_scaling <- function(
    formula,
    data,
    group,
    weights = NULL,
    solvers = c("fn", "qfnb", "pfnb", "proqreg", "profn"),
    reference_solver = "pfnb",
    control = cf_control(crossing_diagnostics = FALSE),
    sample_sizes = c(5000L, 20000L, 100000L),
    seed = 20260719L,
    point_workers = 1L,
    repetitions = 3L,
    warmup = 1L,
    randomize_order = TRUE,
    rss_poll_interval_ms = 25L,
    timeout_seconds = Inf,
    checkpoint_path = NULL,
    resume = TRUE) {
  validate_cf_control(control)
  solvers <- unique(as.character(solvers))
  invalid <- setdiff(solvers, supported_qr_solvers())
  if (length(invalid)) {
    stop("unsupported solvers: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (!reference_solver %in% solvers) {
    stop("reference_solver must be included in solvers", call. = FALSE)
  }
  sample_sizes <- sort(unique(vapply(
    sample_sizes, assert_scalar_integer, integer(1L), name = "sample_sizes",
    minimum = 1L
  )))
  repetitions <- assert_scalar_integer(repetitions, "repetitions", 1L)
  warmup <- assert_scalar_integer(warmup, "warmup", 0L)
  point_workers_requested <- assert_scalar_integer(
    point_workers, "point_workers", 1L
  )
  point_workers <- min(point_workers_requested, 2L)
  for (solver in solvers) {
    validate_execution_parallelism(
      "qr", solver, control, point_workers_requested
    )
  }
  rss_poll_interval_ms <- assert_scalar_integer(
    rss_poll_interval_ms, "rss_poll_interval_ms", 1L
  )
  if (length(timeout_seconds) != 1L || is.na(timeout_seconds) ||
      timeout_seconds <= 0) {
    stop("timeout_seconds must be positive or Inf", call. = FALSE)
  }
  resume <- assert_scalar_logical(resume, "resume")
  prepared_full <- prepare_cf_data(formula, data, group, weights)
  if (any(sample_sizes > prepared_full$n)) {
    stop("sample_sizes cannot exceed the prepared analysis sample", call. = FALSE)
  }
  checkpoint_signature <- object_md5(list(
    schema_version = 3L,
    package_version = as.character(
      utils::packageVersion("scalableCounterfactual")
    ),
    data_fingerprint = bootstrap_data_fingerprint(prepared_full),
    extension_registry = extension_registry_fingerprint(list(
      qr = intersect(solvers, names(registry_snapshot()$qr)),
      linear = character(),
      distribution = character()
    )),
    solvers = sort(solvers),
    reference_solver = reference_solver,
    control = control,
    sample_sizes = sample_sizes,
    seed = seed,
    point_workers_requested = point_workers_requested,
    point_workers = point_workers,
    repetitions = repetitions,
    warmup = warmup,
    randomize_order = randomize_order,
    rss_poll_interval_ms = rss_poll_interval_ms,
    timeout_seconds = timeout_seconds,
    runtime = stats::setNames(lapply(sort(solvers), function(solver) {
      bootstrap_runtime_identity("qr", solver, control)
    }), sort(solvers))
  ))
  checkpoint <- NULL
  if (!is.null(checkpoint_path) && isTRUE(resume) && file.exists(checkpoint_path)) {
    checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(error) NULL)
    if (is.null(checkpoint)) {
      warning("scaling benchmark checkpoint is unreadable; recomputing it",
              call. = FALSE)
    }
    if (!is.null(checkpoint) &&
        !identical(checkpoint$signature, checkpoint_signature)) {
      stop("scaling benchmark checkpoint is incompatible with this run",
           call. = FALSE)
    }
    if (!is.null(checkpoint) && !valid_scaling_checkpoint(
      checkpoint, checkpoint_signature, sample_sizes, solvers, repetitions,
      control$reported_quantiles, control$conditional_quantiles,
      ncol(prepared_full$X0), point_workers_requested, point_workers
    )) {
      warning("scaling benchmark checkpoint is malformed; recomputing it",
              call. = FALSE)
      checkpoint <- NULL
    }
  }
  rows <- checkpoint$rows %||% list()
  fits <- checkpoint$fits %||% list()
  completed_keys <- names(fits) %||% character()
  save_checkpoint <- function() {
    if (is.null(checkpoint_path)) return(invisible(NULL))
    atomic_save_rds(list(
      signature = checkpoint_signature,
      rows = rows,
      fits = fits,
      updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ), checkpoint_path, compress = FALSE)
  }
  run_id <- 0L
  for (sample_n in sample_sizes) {
    prepared <- if (sample_n == prepared_full$n) {
      prepared_full
    } else {
      subsample_prepared_data(prepared_full, sample_n, seed + sample_n)
    }
    has_recorded_run <- any(startsWith(
      completed_keys, paste0(sample_n, "::")
    ))
    if (!has_recorded_run && warmup > 0L) {
      for (iteration in seq_len(warmup)) {
        set.seed(seed + sample_n + iteration * 10000L)
        order <- if (isTRUE(randomize_order)) sample(solvers) else solvers
        for (solver in order) {
          run_id <- run_id + 1L
          solver_seed <- seed + sample_n +
            match(solver, supported_qr_solvers()) * 1000L
          run_isolated_qr_process(
            list(
              prepared = prepared, solver = solver, control = control,
              point_workers = point_workers, seed = solver_seed
            ),
            poll_interval_ms = rss_poll_interval_ms,
            timeout_seconds = timeout_seconds
          )
        }
      }
    }
    for (repetition in seq_len(repetitions)) {
      iteration <- warmup + repetition
      set.seed(seed + sample_n + iteration * 10000L)
      order <- if (isTRUE(randomize_order)) sample(solvers) else solvers
      for (solver in order) {
        key <- paste(sample_n, repetition, solver, sep = "::")
        if (key %in% completed_keys) next
        run_id <- run_id + 1L
        solver_seed <- seed + sample_n +
          match(solver, supported_qr_solvers()) * 1000L
        result <- run_isolated_qr_process(
          list(
            prepared = prepared, solver = solver, control = control,
            point_workers = point_workers, seed = solver_seed
          ),
          poll_interval_ms = rss_poll_interval_ms,
          timeout_seconds = timeout_seconds
        )
        fits[[key]] <- result
        rows[[length(rows) + 1L]] <- data.frame(
          sample_n = sample_n,
          solver = solver,
          repetition = repetition,
          run_order = match(solver, order),
          status = result$status,
          elapsed_seconds = result$elapsed_seconds %||% NA_real_,
          process_elapsed_seconds = result$process_elapsed_seconds %||% NA_real_,
          peak_r_heap_mb = result$peak_r_heap_mb %||% NA_real_,
          peak_process_rss_mb = result$peak_process_rss_mb %||% NA_real_,
          rss_poll_interval_ms = rss_poll_interval_ms,
          point_workers_requested = point_workers_requested,
          point_workers = point_workers,
          error = result$error %||% NA_character_,
          warning_count = result$warning_count %||% NA_integer_,
          warnings = result$warnings %||% NA_character_,
          package_version = as.character(
            utils::packageVersion("scalableCounterfactual")
          ),
          R_version = R.version.string,
          quantreg_version = as.character(utils::packageVersion("quantreg")),
          platform = R.version$platform,
          stringsAsFactors = FALSE
        )
        completed_keys <- c(completed_keys, key)
        save_checkpoint()
      }
    }
  }
  raw <- do.call(rbind, rows)
  raw$reference_solver <- reference_solver
  raw$max_abs_effect_difference <- NA_real_
  raw$max_abs_group0_coefficient_difference <- NA_real_
  raw$max_abs_group1_coefficient_difference <- NA_real_
  for (i in seq_len(nrow(raw))) {
    key <- paste(raw$sample_n[[i]], raw$repetition[[i]], raw$solver[[i]],
                 sep = "::")
    reference_key <- paste(
      raw$sample_n[[i]], raw$repetition[[i]], reference_solver, sep = "::"
    )
    fit <- fits[[key]]
    reference <- fits[[reference_key]]
    if (!is.null(fit) && !is.null(reference) && fit$status == "ok" &&
        reference$status == "ok") {
      raw$max_abs_effect_difference[[i]] <- max(
        abs(fit$effects - reference$effects)
      )
      raw$max_abs_group0_coefficient_difference[[i]] <- safe_max_abs_difference(
        fit$group0_coefficients, reference$group0_coefficients
      )
      raw$max_abs_group1_coefficient_difference[[i]] <- safe_max_abs_difference(
        fit$group1_coefficients, reference$group1_coefficients
      )
    }
  }
  list(
    raw = raw,
    summary = summarize_scaling_benchmark(raw, reference_solver),
    checkpoint_path = checkpoint_path,
    checkpoint_signature = checkpoint_signature
  )
}
