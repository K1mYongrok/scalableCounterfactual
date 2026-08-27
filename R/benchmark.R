extract_qr_coefficients <- function(point, group) {
  fit <- point$fits[[group]]
  if (!inherits(fit, "cf_qr_fit")) return(NULL)
  fit$coefficients
}

#' Benchmark QR solvers on one frozen design
#'
#' Every solver receives the same prepared design matrices, outcomes, weights,
#' quantile grid, and point-worker setting. Timing excludes CSV reading and model
#' matrix preparation.
#'
#' @param formula Outcome and covariate formula.
#' @param data Analysis data.
#' @param group Binary group column or vector.
#' @param weights Optional positive sampling weights.
#' @param solvers Solvers to benchmark.
#' @param reference_solver Solver used for numerical differences.
#' @param control A [cf_control()] object.
#' @param sample_n Optional stratified benchmark sample size.
#' @param seed Sampling seed.
#' @param point_workers Identical worker count for every solver.
#' @return A data frame of timing, memory, convergence, and numerical parity.
#' @export
benchmark_qr_solvers <- function(
    formula,
    data,
    group,
    weights = NULL,
    solvers = supported_qr_solvers(),
    reference_solver = "br",
    control = cf_control(),
    sample_n = NULL,
    seed = 20260719L,
    point_workers = 1L) {
  validate_cf_control(control)
  solvers <- unique(as.character(solvers))
  invalid <- setdiff(solvers, supported_qr_solvers())
  if (length(invalid)) {
    stop("unsupported solvers: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (!reference_solver %in% solvers) {
    stop("reference_solver must be included in solvers", call. = FALSE)
  }
  prepared <- prepare_cf_data(formula, data, group, weights)
  if (!is.null(sample_n)) {
    prepared <- subsample_prepared_data(prepared, sample_n, seed)
  }
  if (length(prepared$dropped_design_columns)) {
    warning(
      "Dropped columns that were unidentified in at least one group: ",
      paste(prepared$dropped_design_columns, collapse = ", "),
      call. = FALSE
    )
  }
  fits <- stats::setNames(vector("list", length(solvers)), solvers)
  rows <- vector("list", length(solvers))
  for (i in seq_along(solvers)) {
    solver <- solvers[[i]]
    solver_seed <- seed + match(solver, supported_qr_solvers()) * 1000L
    result <- tryCatch({
      measured <- measure_resources(function() {
        estimate_point_prepared(
          prepared,
          model = "qr",
          solver = solver,
          control = control,
          point_workers = point_workers,
          point_seed = solver_seed,
          keep_fits = TRUE
        )
      })
      fits[[solver]] <- measured$value
      fit_objects <- measured$value$fits
      resolved_solvers <- vapply(
        fit_objects, function(x) x$solver, character(1L)
      )
      solver_exact <- vapply(
        fit_objects, function(x) isTRUE(x$solver_exact), logical(1L)
      )
      process_fallbacks <- vapply(
        fit_objects, function(x) length(x$process_fallback_taus), integer(1L)
      )
      convergence_flags <- unlist(lapply(
        fit_objects,
        function(x) x$convergence_flag
      ))
      convergence_available <- any(!is.na(convergence_flags))
      solver_warnings <- unique(unlist(lapply(
        fit_objects,
        function(x) x$warnings
      )))
      condition_estimates <- vapply(
        fit_objects,
        function(x) x$preconditioning_condition_estimate,
        numeric(1L)
      )
      data.frame(
        solver = solver,
        resolved_group0_solver = resolved_solvers[[1L]],
        resolved_group1_solver = resolved_solvers[[2L]],
        group0_exact = solver_exact[[1L]],
        group1_exact = solver_exact[[2L]],
        group0_process_fallbacks = process_fallbacks[[1L]],
        group1_process_fallbacks = process_fallbacks[[2L]],
        status = "ok",
        observations = prepared$n,
        design_columns = ncol(prepared$X0),
        conditional_quantiles = length(control$conditional_quantiles),
        point_workers = point_workers,
        elapsed_seconds = measured$value$elapsed_seconds,
        measured_wrapper_seconds = measured$elapsed_seconds,
        peak_r_heap_mb = max(
          measured$peak_r_heap_mb,
          measured$value$resources$peak_r_heap_mb
        ),
        convergence_diagnostics_available = convergence_available,
        nonzero_convergence_flags = if (convergence_available) {
          sum(convergence_flags != 0, na.rm = TRUE)
        } else {
          NA_integer_
        },
        warning_count = length(solver_warnings),
        warnings = paste(solver_warnings, collapse = " | "),
        preconditioning_method = paste(unique(vapply(
          fit_objects,
          function(x) x$preconditioning_method,
          character(1L)
        )), collapse = ","),
        max_condition_estimate = if (any(is.finite(condition_estimates))) {
          max(condition_estimates[is.finite(condition_estimates)])
        } else {
          NA_real_
        },
        error = NA_character_,
        package_version = as.character(
          utils::packageVersion("scalableCounterfactual")
        ),
        R_version = R.version.string,
        quantreg_version = as.character(utils::packageVersion("quantreg")),
        platform = R.version$platform,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        solver = solver,
        resolved_group0_solver = NA_character_,
        resolved_group1_solver = NA_character_,
        group0_exact = NA,
        group1_exact = NA,
        group0_process_fallbacks = NA_integer_,
        group1_process_fallbacks = NA_integer_,
        status = "error",
        observations = prepared$n,
        design_columns = ncol(prepared$X0),
        conditional_quantiles = length(control$conditional_quantiles),
        point_workers = point_workers,
        elapsed_seconds = NA_real_,
        measured_wrapper_seconds = NA_real_,
        peak_r_heap_mb = NA_real_,
        convergence_diagnostics_available = NA,
        nonzero_convergence_flags = NA_integer_,
        warning_count = NA_integer_,
        warnings = NA_character_,
        preconditioning_method = NA_character_,
        max_condition_estimate = NA_real_,
        error = conditionMessage(e),
        package_version = as.character(
          utils::packageVersion("scalableCounterfactual")
        ),
        R_version = R.version.string,
        quantreg_version = as.character(utils::packageVersion("quantreg")),
        platform = R.version$platform,
        stringsAsFactors = FALSE
      )
    })
    rows[[i]] <- result
  }
  reference <- fits[[reference_solver]]
  output <- do.call(rbind, rows)
  output$reference_solver <- reference_solver
  output$max_abs_effect_difference <- vapply(solvers, function(solver) {
    if (is.null(reference) || is.null(fits[[solver]])) return(NA_real_)
    max(abs(fits[[solver]]$effects - reference$effects))
  }, numeric(1L))
  output$max_abs_group0_coefficient_difference <- vapply(solvers, function(solver) {
    if (is.null(reference) || is.null(fits[[solver]])) return(NA_real_)
    safe_max_abs_difference(
      extract_qr_coefficients(fits[[solver]], "group0"),
      extract_qr_coefficients(reference, "group0")
    )
  }, numeric(1L))
  output$max_abs_group1_coefficient_difference <- vapply(solvers, function(solver) {
    if (is.null(reference) || is.null(fits[[solver]])) return(NA_real_)
    safe_max_abs_difference(
      extract_qr_coefficients(fits[[solver]], "group1"),
      extract_qr_coefficients(reference, "group1")
    )
  }, numeric(1L))
  attr(output, "fits") <- fits
  output
}

#' Repeated QR solver benchmark
#'
#' Repeats [benchmark_qr_solvers()] on the same stratified sample. Solver order
#' can be randomized while solver-specific seeds remain fixed.
#'
#' @inheritParams benchmark_qr_solvers
#' @param repetitions Number of recorded repetitions.
#' @param warmup Number of discarded warmup repetitions.
#' @param randomize_order Randomize solver execution order in each repetition.
#' @return A list with repeat-level `raw` results and a solver-level `summary`.
#' @export
benchmark_qr_solvers_repeated <- function(
    formula,
    data,
    group,
    weights = NULL,
    solvers = supported_qr_solvers(),
    reference_solver = "br",
    control = cf_control(),
    sample_n = NULL,
    seed = 20260719L,
    point_workers = 1L,
    repetitions = 5L,
    warmup = 0L,
    randomize_order = TRUE) {
  validate_cf_control(control)
  repetitions <- assert_scalar_integer(repetitions, "repetitions", 1L)
  warmup <- assert_scalar_integer(warmup, "warmup", 0L)
  solvers <- unique(as.character(solvers))

  run_once <- function(rep_id, recorded) {
    order_seed <- seed + rep_id * 10000L + if (recorded) 0L else 500000L
    set.seed(order_seed)
    solver_order <- if (isTRUE(randomize_order)) sample(solvers) else solvers
    result <- benchmark_qr_solvers(
      formula = formula,
      data = data,
      group = group,
      weights = weights,
      solvers = solver_order,
      reference_solver = reference_solver,
      control = control,
      sample_n = sample_n,
      seed = seed,
      point_workers = point_workers
    )
    result$repetition <- if (recorded) rep_id else NA_integer_
    result$run_order <- match(result$solver, solver_order)
    result
  }

  if (warmup > 0L) {
    invisible(lapply(seq_len(warmup), run_once, recorded = FALSE))
  }
  raw <- do.call(rbind, lapply(
    seq_len(repetitions), run_once, recorded = TRUE
  ))
  summary <- do.call(rbind, lapply(solvers, function(solver_name) {
    all_rows <- raw[raw$solver == solver_name, , drop = FALSE]
    rows <- all_rows[all_rows$status == "ok", , drop = FALSE]
    if (!nrow(rows)) {
      return(data.frame(
        solver = solver_name,
        repetitions = 0L,
        failed_repetitions = nrow(all_rows),
        mean_seconds = NA_real_,
        median_seconds = NA_real_,
        min_seconds = NA_real_,
        max_seconds = NA_real_,
        sd_seconds = NA_real_,
        mean_peak_r_heap_mb = NA_real_,
        max_abs_effect_difference = NA_real_,
        warning_runs = 0L,
        process_fallback_runs = 0L,
        stringsAsFactors = FALSE
      ))
    }
    numerical_differences <- rows$max_abs_effect_difference
    data.frame(
      solver = rows$solver[[1L]],
      repetitions = nrow(rows),
      failed_repetitions = nrow(all_rows) - nrow(rows),
      mean_seconds = mean(rows$elapsed_seconds),
      median_seconds = stats::median(rows$elapsed_seconds),
      min_seconds = min(rows$elapsed_seconds),
      max_seconds = max(rows$elapsed_seconds),
      sd_seconds = if (nrow(rows) > 1L) stats::sd(rows$elapsed_seconds) else NA_real_,
      mean_peak_r_heap_mb = mean(rows$peak_r_heap_mb),
      max_abs_effect_difference = if (any(is.finite(numerical_differences))) {
        max(numerical_differences[is.finite(numerical_differences)])
      } else {
        NA_real_
      },
      warning_runs = sum(rows$warning_count > 0L, na.rm = TRUE),
      process_fallback_runs = sum(
        rows$group0_process_fallbacks > 0L |
          rows$group1_process_fallbacks > 0L,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  reference_median <- summary$median_seconds[
    summary$solver == reference_solver
  ]
  if (length(reference_median) != 1L) {
    stop("reference solver failed in every recorded repetition", call. = FALSE)
  }
  summary$relative_to_reference_median <-
    summary$median_seconds / reference_median
  summary <- summary[order(is.na(summary$median_seconds), summary$median_seconds),
                     , drop = FALSE]
  list(raw = raw, summary = summary)
}
