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
  effective_reps <- as.integer(colSums(is.finite(draws)))
  standard_error <- vapply(seq_len(ncol(draws)), function(index) {
    values <- draws[is.finite(draws[, index]), index]
    if (length(values) < 2L) return(NA_real_)
    if (isTRUE(control$robust_se)) {
      diff(stats::quantile(values, c(0.25, 0.75))) / 1.34
    } else {
      stats::sd(values)
    }
  }, numeric(1L))
  standard_error[!identified | effective_reps < 2L] <- NA_real_
  point_critical <- stats::qnorm(1 - control$alpha / 2)
  uniform_support <- rep(FALSE, length(estimate))
  if (nrow(draws) >= 2L) {
    uniform_support <- identified & is.finite(standard_error) &
      effective_reps == nrow(draws)
  }
  positive_scale <- uniform_support & standard_error > 0
  uniform_critical <- if (any(positive_scale)) {
    studentized <- abs(sweep(
      draws[, positive_scale, drop = FALSE],
      2L,
      estimate[positive_scale],
      "-"
    ) / matrix(
      standard_error[positive_scale],
      nrow(draws),
      sum(positive_scale),
      byrow = TRUE
    ))
    maximum_t <- apply(studentized, 1L, max)
    as.numeric(stats::quantile(
      maximum_t, 1 - control$alpha, names = FALSE
    ))
  } else if (any(uniform_support)) {
    0
  } else {
    NA_real_
  }
  uniform_lower <- uniform_upper <- rep(NA_real_, length(estimate))
  if (is.finite(uniform_critical)) {
    uniform_lower[uniform_support] <- estimate[uniform_support] -
      uniform_critical * standard_error[uniform_support]
    uniform_upper[uniform_support] <- estimate[uniform_support] +
      uniform_critical * standard_error[uniform_support]
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
    uniform_lower = uniform_lower,
    uniform_upper = uniform_upper,
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
