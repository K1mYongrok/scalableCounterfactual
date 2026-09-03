#' Benchmark non-QR conditional-model backends
#'
#' Every backend receives the same prepared outcomes, design matrices, group
#' split, weights, thresholds, and execution settings. Use
#' [benchmark_qr_solvers()] for quantile-regression solvers.
#'
#' @inheritParams benchmark_qr_solvers
#' @param model One of `loc`, `locsca`, `logit`, `probit`, `cloglog`, or `lpm`.
#' @param backends Linear backends for `loc`, `locsca`, and `lpm`, or binary
#'   distribution-regression backends for `logit`, `probit`, and `cloglog`.
#' @param reference_backend Backend used for numerical differences.
#' @return A data frame of timing, memory, diagnostics, and numerical parity.
#' @export
benchmark_conditional_backends <- function(
    formula,
    data,
    group,
    weights = NULL,
    model = c("loc", "locsca", "logit", "probit", "cloglog", "lpm"),
    backends = NULL,
    reference_backend = NULL,
    control = cf_control(),
    sample_n = NULL,
    seed = 20260719L,
    point_workers = 1L) {
  validate_cf_control(control)
  model <- match.arg(model)
  linear_model <- model %in% c("loc", "locsca", "lpm")
  supported <- if (linear_model) {
    supported_linear_backends()
  } else {
    supported_dr_backends()
  }
  if (is.null(backends)) backends <- supported
  backends <- unique(as.character(backends))
  invalid <- setdiff(backends, supported)
  if (length(invalid)) {
    stop("unsupported backends: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (is.null(reference_backend)) {
    reference_backend <- if (linear_model) "qr" else "glm"
  }
  if (!reference_backend %in% backends) {
    stop("reference_backend must be included in backends", call. = FALSE)
  }
  point_workers_requested <- assert_scalar_integer(
    point_workers, "point_workers", 1L
  )
  point_workers <- min(point_workers_requested, 2L)
  for (backend in backends) {
    backend_control <- control
    if (linear_model) {
      backend_control$linear_backend <- backend
    } else {
      backend_control$dr_backend <- backend
    }
    validate_execution_parallelism(
      model, NA_character_, backend_control, point_workers_requested
    )
  }
  prepared <- prepare_cf_data(formula, data, group, weights)
  if (!is.null(sample_n)) {
    prepared <- subsample_prepared_data(prepared, sample_n, seed)
  }

  fits <- stats::setNames(vector("list", length(backends)), backends)
  rows <- lapply(seq_along(backends), function(index) {
    backend <- backends[[index]]
    backend_control <- control
    if (linear_model) {
      backend_control$linear_backend <- backend
    } else {
      backend_control$dr_backend <- backend
    }
    result <- tryCatch({
      measured <- measure_resources(function() {
        estimate_point_prepared(
          prepared,
          model = model,
          solver = NA_character_,
          control = backend_control,
          point_workers = point_workers,
          point_seed = seed,
          keep_fits = TRUE
        )
      })
      fits[[backend]] <<- measured$value
      resolved <- vapply(
        measured$value$fits, function(fit) fit$backend, character(1L)
      )
      grid_sizes <- vapply(measured$value$fits, function(fit) {
        if (!is.null(fit$thresholds)) length(fit$thresholds) else
          length(fit$taus)
      }, integer(1L))
      effective_threshold_workers <- if (linear_model) {
        c(group0 = 1L, group1 = 1L)
      } else {
        vapply(
          measured$value$fits, `[[`, integer(1L), "threshold_workers"
        )
      }
      data.frame(
        model = model,
        backend = backend,
        resolved_group0_backend = resolved[[1L]],
        resolved_group1_backend = resolved[[2L]],
        status = "ok",
        observations = prepared$n,
        design_columns = ncol(prepared$X0),
        conditional_grid_size = if (
          identical(grid_sizes[[1L]], grid_sizes[[2L]])
        ) grid_sizes[[1L]] else NA_integer_,
        conditional_grid_size_group0 = grid_sizes[[1L]],
        conditional_grid_size_group1 = grid_sizes[[2L]],
        point_workers_requested = point_workers_requested,
        point_workers = point_workers,
        threshold_workers_requested = if (linear_model) 1L else
          control$dr_workers,
        threshold_workers = if (identical(
          effective_threshold_workers[[1L]],
          effective_threshold_workers[[2L]]
        )) effective_threshold_workers[[1L]] else NA_integer_,
        threshold_workers_group0 = effective_threshold_workers[[1L]],
        threshold_workers_group1 = effective_threshold_workers[[2L]],
        elapsed_seconds = measured$value$elapsed_seconds,
        measured_wrapper_seconds = measured$elapsed_seconds,
        peak_r_heap_mb = max(
          measured$peak_r_heap_mb,
          measured$value$resources$peak_r_heap_mb
        ),
        warning_count = length(measured$value$warnings),
        warnings = paste(measured$value$warnings, collapse = " | "),
        error = NA_character_,
        package_version = as.character(
          utils::packageVersion("scalableCounterfactual")
        ),
        stringsAsFactors = FALSE
      )
    }, error = function(error) {
      data.frame(
        model = model,
        backend = backend,
        resolved_group0_backend = NA_character_,
        resolved_group1_backend = NA_character_,
        status = "error",
        observations = prepared$n,
        design_columns = ncol(prepared$X0),
        conditional_grid_size = NA_integer_,
        conditional_grid_size_group0 = NA_integer_,
        conditional_grid_size_group1 = NA_integer_,
        point_workers_requested = point_workers_requested,
        point_workers = point_workers,
        threshold_workers_requested = if (linear_model) 1L else
          control$dr_workers,
        threshold_workers = NA_integer_,
        threshold_workers_group0 = NA_integer_,
        threshold_workers_group1 = NA_integer_,
        elapsed_seconds = NA_real_,
        measured_wrapper_seconds = NA_real_,
        peak_r_heap_mb = NA_real_,
        warning_count = NA_integer_,
        warnings = NA_character_,
        error = conditionMessage(error),
        package_version = as.character(
          utils::packageVersion("scalableCounterfactual")
        ),
        stringsAsFactors = FALSE
      )
    })
    result
  })
  output <- do.call(rbind, rows)
  reference <- fits[[reference_backend]]
  output$reference_backend <- reference_backend
  output$max_abs_effect_difference <- vapply(backends, function(backend) {
    if (is.null(reference) || is.null(fits[[backend]])) return(NA_real_)
    max(abs(fits[[backend]]$effects - reference$effects))
  }, numeric(1L))
  attr(output, "fits") <- fits
  output
}
