fit_conditional_model <- function(
    X, y, weights, model, solver, control, censoring = NULL, event = NULL) {
  model <- match.arg(model, supported_cf_models())
  if (model == "qr") {
    return(fit_weighted_qr(
      X,
      y,
      weights,
      control$conditional_quantiles,
      solver,
      precondition = control$qr_precondition,
      onestep_first_solver = control$onestep_first_solver,
      onestep_bandwidth = control$onestep_bandwidth,
      gpu_control = control
    ))
  }
  if (model == "cqr") {
    return(fit_weighted_cqr(
      X, y, weights, censoring, control$conditional_quantiles, solver,
      right = control$cqr_right,
      nsteps = control$cqr_nsteps,
      first_cut = control$cqr_first_cut,
      later_cut = control$cqr_later_cut,
      precondition = control$qr_precondition,
      dr_backend = control$dr_backend,
      dr_maxit = control$dr_maxit,
      dr_tolerance = control$dr_tolerance
    ))
  }
  if (model == "loc") {
    return(fit_location_model(
      X,
      y,
      weights,
      control$conditional_quantiles,
      linear_backend = control$linear_backend
    ))
  }
  if (model == "locsca") {
    return(fit_location_scale_model(
      X,
      y,
      weights,
      control$conditional_quantiles,
      linear_backend = control$linear_backend
    ))
  }
  if (model == "cox") {
    return(fit_weighted_cox(X, y, weights, event))
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
