# Cox fitting is delegated to survival; weighted Breslow and Kaplan-Meier
# distribution calculations are implemented here. See inst/provenance/METHODS.md.
fit_weighted_cox <- function(X, y, weights, event = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  weights <- normalize_weights(weights)
  if (is.null(event)) event <- rep(1L, length(y))
  event <- as.integer(event)
  if (nrow(X) != length(y) || length(weights) != length(y) ||
      length(event) != length(y)) {
    stop("Cox inputs have incompatible sizes", call. = FALSE)
  }
  if (any(y < 0)) {
    stop("Cox duration regression requires a nonnegative outcome", call. = FALSE)
  }
  if (!all(event %in% c(0L, 1L)) || !any(event == 1L)) {
    stop("Cox event must contain 0/1 and at least one event", call. = FALSE)
  }
  intercept <- match("(Intercept)", colnames(X))
  covariate_columns <- setdiff(seq_len(ncol(X)), intercept)
  if (!length(covariate_columns)) {
    stop("Cox regression requires at least one non-intercept covariate",
         call. = FALSE)
  }
  Z <- X[, covariate_columns, drop = FALSE]
  fit_warnings <- character()
  cox_control <- survival::coxph.control()
  fitted <- withCallingHandlers(
    survival::coxph.fit(
      x = Z,
      y = survival::Surv(y, event),
      strata = NULL,
      offset = rep(0, length(y)),
      init = NULL,
      control = cox_control,
      weights = weights,
      method = "breslow",
      rownames = as.character(seq_along(y)),
      resid = FALSE
    ),
    warning = function(condition) {
      fit_warnings <<- c(fit_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  slope <- as.numeric(fitted$coefficients)
  if (length(slope) != ncol(Z) || any(!is.finite(slope))) {
    stop("Cox regression returned invalid coefficients", call. = FALSE)
  }
  nonconvergence_warning <- any(grepl(
    "did not converge|ran out of iterations|failed to converge",
    fit_warnings, ignore.case = TRUE
  ))
  cox_iterations <- if (length(fitted$iter)) {
    max(as.integer(fitted$iter), na.rm = TRUE)
  } else {
    NA_integer_
  }
  iteration_limit <- is.finite(cox_iterations) &&
    cox_iterations >= cox_control$iter.max
  if (nonconvergence_warning || iteration_limit) {
    stop(
      "Cox regression failed to converge after ", cox_iterations,
      " iteration(s)",
      if (length(fit_warnings)) paste0(": ", paste(fit_warnings, collapse = " | "))
      else "",
      call. = FALSE
    )
  }
  infinite_warning <- any(grepl(
    "coefficient.*infinite|loglik converged before|infinite coefficient",
    fit_warnings, ignore.case = TRUE
  ))
  variance_problem <- is.null(fitted$var) || any(!is.finite(fitted$var))
  convergence_flag <- if (infinite_warning || variance_problem) 2L else 0L
  if (variance_problem) {
    fit_warnings <- c(fit_warnings, "non-finite Cox coefficient variance")
  }
  coefficients <- numeric(ncol(X))
  coefficients[covariate_columns] <- slope
  raw_lp <- drop(X %*% coefficients)
  lp_shift <- max(raw_lp)
  relative_risk <- exp(raw_lp - lp_shift)

  ordering <- order(y)
  sorted_time <- y[ordering]
  risk_from_time <- rev(cumsum(rev(weights[ordering] * relative_risk[ordering])))
  event_times <- sort(unique(y[event == 1L]))
  event_weight <- vapply(event_times, function(time) {
    sum(weights[event == 1L & y == time])
  }, numeric(1L))
  risk_set <- risk_from_time[match(event_times, sorted_time)]
  increments <- event_weight / risk_set
  if (any(!is.finite(increments)) || any(increments <= 0)) {
    stop("Cox baseline-hazard calculation failed", call. = FALSE)
  }
  baseline_hazard <- cumsum(increments)
  structure(list(
    model = "cox",
    coefficients = coefficients,
    event_times = event_times,
    baseline_hazard = baseline_hazard,
    linear_predictor_shift = lp_shift,
    backend = "survival::coxph.fit(breslow)",
    iterations = cox_iterations,
    converged = convergence_flag == 0L,
    convergence_flag = convergence_flag,
    event_count = sum(event == 1L),
    warnings = unique(fit_warnings)
  ), class = c("cf_cox_fit", "cf_conditional_fit"))
}

cox_relative_risk <- function(fit, X) {
  predictor <- drop(as.matrix(X) %*% fit$coefficients) -
    fit$linear_predictor_shift
  exp(pmin(700, predictor))
}

cox_boundary_result <- function(
    probs, upper_cdf, identified, final_time, policy, context) {
  boundary <- probs > upper_cdf
  if (any(boundary) && identical(policy, "error")) {
    stop(
      context, " quantile(s) above the identified CDF limit ",
      signif(upper_cdf, 6), ": ",
      paste(signif(probs[boundary], 6), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(boundary)) {
    identified[boundary] <- if (identical(policy, "cap")) final_time else NA_real_
  }
  list(values = identified, boundary = boundary)
}

weighted_km_quantiles <- function(
    y, event, weights, probs, boundary = c("na", "error", "cap")) {
  boundary <- match.arg(boundary)
  y <- as.numeric(y)
  event <- as.integer(event)
  weights <- normalize_weights(weights)
  ordering <- order(y)
  sorted_time <- y[ordering]
  risk_weight <- rev(cumsum(rev(weights[ordering])))
  event_times <- sort(unique(y[event == 1L]))
  event_weight <- vapply(event_times, function(time) {
    sum(weights[event == 1L & y == time])
  }, numeric(1L))
  risk_set <- risk_weight[match(event_times, sorted_time)]
  survival <- cumprod(pmax(0, 1 - event_weight / risk_set))
  cdf <- 1 - survival
  upper_cdf <- cdf[[length(cdf)]]
  identified <- vapply(probs, function(probability) {
    hit <- which(cdf >= probability)[1L]
    if (is.na(hit)) NA_real_ else event_times[[hit]]
  }, numeric(1L))
  resolved <- cox_boundary_result(
    probs, upper_cdf, identified, event_times[[length(event_times)]],
    boundary, "Kaplan-Meier"
  )
  attr(resolved$values, "marginal_diagnostics") <- list(
    method = "weighted_kaplan_meier",
    passes = 1L, histogram_bins = NA_integer_, candidate_draws = NA_real_,
    estimated_matrix_mb = NA_real_,
    boundary_quantiles = sum(resolved$boundary), identified_cdf_max = upper_cdf
  )
  resolved$values
}

cox_marginal_quantiles <- function(
    fit, X, weights, probs, boundary = c("na", "error", "cap")) {
  boundary <- match.arg(boundary)
  weights <- normalize_weights(weights)
  risk <- cox_relative_risk(fit, X)
  hazards <- fit$baseline_hazard
  times <- fit$event_times
  cache <- new.env(parent = emptyenv())
  evaluate_index <- function(index) {
    key <- as.character(index)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      probability <- -expm1(-hazards[[index]] * risk)
      assign(key, sum(weights * probability) / sum(weights), envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
  upper_cdf <- evaluate_index(length(hazards))
  identified <- vapply(probs, function(probability) {
    if (probability > upper_cdf) return(NA_real_)
    lower <- 1L
    upper <- length(hazards)
    while (lower < upper) {
      middle <- as.integer(floor((lower + upper) / 2))
      if (evaluate_index(middle) >= probability) upper <- middle else {
        lower <- middle + 1L
      }
    }
    times[[lower]]
  }, numeric(1L))
  resolved <- cox_boundary_result(
    probs, upper_cdf, identified, times[[length(times)]], boundary,
    "Cox marginal"
  )
  result <- resolved$values
  attr(result, "marginal_diagnostics") <- list(
    method = "cox_cdf_binary_search",
    passes = length(ls(cache, all.names = TRUE)),
    histogram_bins = NA_integer_,
    candidate_draws = NA_real_,
    estimated_matrix_mb = as.numeric(nrow(X)) * length(hazards) * 8 / 1024^2,
    boundary_quantiles = sum(resolved$boundary), identified_cdf_max = upper_cdf
  )
  result
}
