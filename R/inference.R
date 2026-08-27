point_effect_table <- function(point, reps, alpha) {
  effect_names <- rownames(point$effects)
  do.call(rbind, lapply(effect_names, function(effect) {
    data.frame(
      model = point$model,
      solver = point$solver,
      quantile = point$quantiles,
      effect = effect,
      estimate = as.numeric(point$effects[effect, ]),
      identified = is.finite(as.numeric(point$effects[effect, ])),
      std_error = NA_real_,
      pointwise_lower = NA_real_,
      pointwise_upper = NA_real_,
      uniform_lower = NA_real_,
      uniform_upper = NA_real_,
      bootstrap_reps = reps,
      bootstrap_reps_effective = 0L,
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  }))
}

effect_inference <- function(effect, estimate, draws, control, model, solver) {
  identified <- is.finite(estimate)
  effective_reps <- colSums(is.finite(draws))
  standard_error <- if (control$robust_se) {
    apply(draws, 2L, function(x) {
      diff(stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)) / 1.34
    })
  } else {
    apply(draws, 2L, stats::sd, na.rm = TRUE)
  }
  standard_error[!identified | effective_reps < 2L] <- NA_real_
  point_critical <- stats::qnorm(1 - control$alpha / 2)
  safe_se <- ifelse(standard_error > 0, standard_error, NA_real_)
  studentized <- abs(sweep(draws, 2L, estimate, "-") /
    matrix(safe_se, nrow(draws), ncol(draws), byrow = TRUE))
  usable <- identified & apply(studentized, 2L, function(x) all(is.finite(x)))
  uniform_critical <- if (any(usable)) {
    maximum_t <- apply(studentized[, usable, drop = FALSE], 1L, max)
    as.numeric(stats::quantile(
      maximum_t, 1 - control$alpha, names = FALSE, na.rm = TRUE
    ))
  } else {
    NA_real_
  }
  data.frame(
    model = model,
    solver = solver,
    quantile = control$reported_quantiles,
    effect = effect,
    estimate = estimate,
    identified = identified,
    std_error = standard_error,
    pointwise_lower = estimate - point_critical * standard_error,
    pointwise_upper = estimate + point_critical * standard_error,
    uniform_lower = estimate - uniform_critical * standard_error,
    uniform_upper = estimate + uniform_critical * standard_error,
    bootstrap_reps = nrow(draws),
    bootstrap_reps_effective = effective_reps,
    alpha = control$alpha,
    stringsAsFactors = FALSE
  )
}

bootstrap_inference <- function(point, bootstrap, control) {
  do.call(rbind, lapply(c("structure", "composition", "total"), function(effect) {
    effect_inference(
      effect = effect,
      estimate = as.numeric(point$effects[effect, ]),
      draws = bootstrap$effects[[effect]],
      control = control,
      model = point$model,
      solver = point$solver
    )
  }))
}
