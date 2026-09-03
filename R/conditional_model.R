fit_conditional_model <- function(
    X, y, weights, model, solver, control, quantile_frequency = NULL,
    censoring = NULL, event = NULL) {
  model <- match.arg(model, supported_cf_models())
  if (is.null(quantile_frequency)) quantile_frequency <- rep(1, length(y))
  if (!is.numeric(quantile_frequency) ||
      length(quantile_frequency) != length(y) ||
      any(!is.finite(quantile_frequency)) || any(quantile_frequency <= 0)) {
    stop("quantile_frequency must contain one finite positive value per row",
         call. = FALSE)
  }
  if (model == "qr") {
    base_weights <- weights / quantile_frequency
    return(fit_weighted_qr(
      X,
      y,
      base_weights,
      control$conditional_quantiles,
      solver,
      precondition = control$qr_precondition,
      onestep_first_solver = control$onestep_first_solver,
      onestep_bandwidth = control$onestep_bandwidth,
      gpu_control = control,
      frequency = quantile_frequency
    ))
  }
  if (model == "cqr") {
    return(fit_weighted_cqr(
      X, y, weights, censoring, control$full_conditional_quantiles, solver,
      right = control$cqr_right,
      nsteps = control$cqr_nsteps,
      first_cut = control$cqr_first_cut,
      later_cut = control$cqr_later_cut,
      precondition = control$qr_precondition,
      dr_backend = control$dr_backend,
      dr_maxit = control$dr_maxit,
      dr_tolerance = control$dr_tolerance,
      quantile_frequency = quantile_frequency
    ))
  }
  if (model == "loc") {
    return(fit_location_model(
      X,
      y,
      weights,
      control$full_conditional_quantiles,
      linear_backend = control$linear_backend,
      quantile_frequency = quantile_frequency
    ))
  }
  if (model == "locsca") {
    return(fit_location_scale_model(
      X,
      y,
      weights,
      control$full_conditional_quantiles,
      linear_backend = control$linear_backend,
      quantile_frequency = quantile_frequency
    ))
  }
  if (model == "cox") {
    base_weights <- weights / quantile_frequency
    return(fit_weighted_cox(
      X, y, base_weights, event, frequency = quantile_frequency
    ))
  }
  fit_distribution_regression(
    X,
    y,
    weights,
    model,
    control$nreg,
    dr_backend = control$dr_backend,
    linear_backend = control$linear_backend,
    dr_workers = control$dr_workers,
    warm_start = control$dr_warm_start,
    maxit = control$dr_maxit,
    tolerance = control$dr_tolerance,
    precondition = control$dr_precondition,
    control = control
  )
}

fit_group_task <- function(task) {
  measured <- measure_resources(function() {
    fit_conditional_model(
      task$X,
      task$y,
      task$weights,
      task$model,
      task$solver,
      task$control,
      quantile_frequency = task$quantile_frequency,
      censoring = task$censoring,
      event = task$event
    )
  })
  list(
    fit = measured$value,
    elapsed_seconds = measured$elapsed_seconds,
    peak_r_heap_mb = measured$peak_r_heap_mb
  )
}
