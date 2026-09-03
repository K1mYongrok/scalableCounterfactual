library(scalableCounterfactual)

atomic_path <- tempfile("cf_atomic_", fileext = ".rds")
scalableCounterfactual:::atomic_save_rds(list(value = 1L), atomic_path)
scalableCounterfactual:::atomic_save_rds(list(value = 2L), atomic_path)
stopifnot(readRDS(atomic_path)$value == 2L)
unlink(atomic_path)

set.seed(42)
n <- 700L
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.4)
group <- rbinom(n, 1, plogis(0.2 + 0.3 * x1))
weights <- runif(n, 0.5, 2)
y <- 1 + 0.2 * group + 0.7 * x1 - 0.3 * x2 +
  (0.5 + 0.2 * group) * rnorm(n)
test_data <- data.frame(y, x1, x2, group, weights)
control <- cf_control(
  nreg = 9L,
  trimming = 0.05,
  reported_quantiles = c(0.25, 0.5, 0.75)
)
stopifnot(identical(
  scalableCounterfactual:::validate_model_solver("qr", NULL)$solver,
  "auto"
))
stopifnot(identical(formals(fit_weighted_qr)$solver, "auto"))
invalid_integer <- try(cf_control(nreg = Inf), silent = TRUE)
stopifnot(inherits(invalid_integer, "try-error"))

invalid_control <- try(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "loc", control = list()
), silent = TRUE)
stopifnot(inherits(invalid_control, "try-error"))

non_qr_solver <- try(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "logit", solver = "pfnb", control = control
), silent = TRUE)
stopifnot(inherits(non_qr_solver, "try-error"))

fractional_group <- test_data
fractional_group$group[[1L]] <- 0.2
fractional_error <- try(counterfactual_decompose(
  y ~ x1 + x2, fractional_group, "group", "weights",
  model = "loc", control = control
), silent = TRUE)
stopifnot(inherits(fractional_error, "try-error"))

deprecated_control <- suppressWarnings(cf_control(weighted_bootstrap = FALSE))
stopifnot(identical(deprecated_control$bootstrap_scheme, "empirical"))

br <- counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  model = "qr",
  solver = "br",
  control = control
)
chunked_br <- counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  model = "qr",
  solver = "br",
  control = cf_control(
    nreg = 9L,
    trimming = 0.05,
    reported_quantiles = c(0.25, 0.5, 0.75),
    marginal_method = "chunked",
    marginal_chunk_rows = 73L,
    marginal_histogram_bins = 256L,
    marginal_candidate_max = 100000L
  )
)
stopifnot(max(abs(br$point$effects - chunked_br$point$effects)) < 1e-10)
stopifnot(all(chunked_br$point$marginal_diagnostics$method == "chunked"))
stopifnot(nrow(chunked_br$point$crossing_diagnostics) == 3L)
stopifnot(all(
  chunked_br$point$crossing_diagnostics$crossing_pair_share >= 0 &
    chunked_br$point$crossing_diagnostics$crossing_pair_share <= 1
))
discrete_fit <- structure(
  list(coefficients = 0, residual_quantiles = 1:4),
  class = c("cf_loc_fit", "cf_conditional_fit")
)
discrete_X <- matrix(1, 3, 1)
discrete_probs <- c(0.25, 0.5, 0.75)
matrix_quantiles <- scalableCounterfactual:::marginal_quantiles(
  discrete_fit, discrete_X, rep(1, 3), discrete_probs,
  cf_control(
    reported_quantiles = discrete_probs, marginal_method = "matrix",
    crossing_diagnostics = FALSE
  )
)
chunked_quantiles <- scalableCounterfactual:::marginal_quantiles(
  discrete_fit, discrete_X, rep(1, 3), discrete_probs,
  cf_control(
    reported_quantiles = discrete_probs, marginal_method = "chunked",
    marginal_chunk_rows = 2L, marginal_histogram_bins = 256L,
    crossing_diagnostics = FALSE
  )
)
stopifnot(identical(
  as.numeric(matrix_quantiles), as.numeric(chunked_quantiles)
))
diagnostic_output <- tempfile("cf_diagnostic_output_")
write_cf_outputs(chunked_br, diagnostic_output)
stopifnot(all(file.exists(file.path(diagnostic_output, c(
  "marginalization_diagnostics.csv",
  "quantile_crossing_diagnostics.csv"
)))))
unlink(diagnostic_output, recursive = TRUE)
stopifnot(identical(br$metadata$group, "group"))
stopifnot(identical(br$metadata$weights, "weights"))
stopifnot(!any(c(
  "prepare_wage_gap_data", "wage_gap_formula"
) %in% getNamespaceExports("scalableCounterfactual")))
pfnb <- counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  model = "qr",
  solver = "pfnb",
  control = control
)
stopifnot(max(abs(br$point$effects - pfnb$point$effects)) < 1e-3)
stopifnot(max(abs(pfnb$point$identity_residual)) < 1e-10)

process_control <- cf_control(
  nreg = 21L,
  trimming = 0.02,
  reported_quantiles = c(0.25, 0.5, 0.75)
)
br_process <- counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  solver = "br", control = process_control
)
proqreg <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  solver = "proqreg", control = process_control
))
profn <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  solver = "profn", control = process_control
))
onestep <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  solver = "onestep", control = process_control
))
automatic <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  solver = "auto", control = process_control
))
stopifnot(max(abs(br_process$point$effects - proqreg$point$effects)) < 1e-6)
stopifnot(max(abs(br_process$point$effects - profn$point$effects)) < 1e-4)
stopifnot(max(abs(onestep$point$identity_residual)) < 1e-10)
stopifnot(!onestep$metadata$solver_group0_exact)
stopifnot(!onestep$point$fits$group0$preconditioned)
stopifnot(identical(
  onestep$metadata$stata_source_version,
  "qrprocess 1.1.3"
))
stopifnot(identical(
  onestep$metadata$stata_source_commit,
  "ec56830ef9c84ce54411ad59c5ce94535847d9df"
))
stopifnot(all(is.finite(unlist(
  onestep$point$fits$group0$onestep_inverse_jacobians[-11L]
))))
stopifnot(identical(automatic$metadata$solver_group0, "qfnb"))
stopifnot(identical(automatic$metadata$solver_group1, "qfnb"))

coarse_onestep <- try(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  solver = "onestep", control = control
), silent = TRUE)
stopifnot(inherits(coarse_onestep, "try-error"))
stopifnot(identical(
  scalableCounterfactual:::stata_weighted_quantile_type2(
    1:4, rep(1, 4), 0.5
  ),
  2.5
))
stopifnot(identical(
  scalableCounterfactual:::stata_weighted_quantile_type2(
    1:3, c(1, 1, 2), 0.5
  ),
  2.5
))

br_rescaled_weights <- counterfactual_decompose(
  y ~ x1 + x2,
  transform(test_data, weights = weights * 1000),
  group = "group",
  weights = "weights",
  model = "qr",
  solver = "br",
  control = control
)
stopifnot(max(abs(
  br$point$effects - br_rescaled_weights$point$effects
)) < 1e-12)
stopifnot(max(abs(
  br$point$diagnostics - br_rescaled_weights$point$diagnostics
)) < 1e-12)

br_unscaled <- counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  model = "qr",
  solver = "br",
  control = cf_control(
    nreg = 9L,
    trimming = 0.05,
    reported_quantiles = c(0.25, 0.5, 0.75),
    qr_precondition = FALSE
  )
)
stopifnot(max(abs(br$point$effects - br_unscaled$point$effects)) < 1e-6)

for (model in c("loc", "locsca", "logit", "probit", "cloglog", "lpm")) {
  fit <- counterfactual_decompose(
    y ~ x1 + x2,
    test_data,
    group = "group",
    weights = "weights",
    model = model,
    control = control
  )
  stopifnot(nrow(fit$results) == 9L)
  stopifnot(max(abs(fit$point$identity_residual)) < 1e-10)
}

crossing_fit <- structure(list(
  coefficients = matrix(c(3, 1, 2), nrow = 1L),
  taus = c(0.25, 0.5, 0.75)
), class = c("cf_qr_fit", "cf_conditional_fit"))
crossing_X <- matrix(1, nrow = 4L, ncol = 1L)
raw_draws <- scalableCounterfactual:::predict_conditional_draws(
  crossing_fit, crossing_X,
  cf_control(legacy_qr_shift = FALSE, quantile_noncrossing = "none")
)
rearranged_draws <- scalableCounterfactual:::predict_conditional_draws(
  crossing_fit, crossing_X,
  cf_control(legacy_qr_shift = FALSE, quantile_noncrossing = "rearrange")
)
stopifnot(any(apply(raw_draws, 1L, is.unsorted)))
stopifnot(!any(apply(rearranged_draws, 1L, is.unsorted)))
stopifnot(identical(sort(as.vector(raw_draws)), sort(as.vector(rearranged_draws))))

rearranged_br <- counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "qr", solver = "br",
  control = cf_control(
    nreg = 9L, trimming = 0.05,
    reported_quantiles = c(0.25, 0.5, 0.75),
    quantile_noncrossing = "rearrange"
  )
)
stopifnot(max(abs(br$point$effects - rearranged_br$point$effects)) < 1e-12)
stopifnot(nrow(rearranged_br$point$crossing_diagnostics) == 6L)
stopifnot(all(
  rearranged_br$point$crossing_diagnostics$crossing_pairs[
    rearranged_br$point$crossing_diagnostics$stage == "corrected"
  ] == 0
))

cqr_x2 <- rnorm(n)
latent_cqr <- 0.8 + 0.2 * group + 0.7 * x1 - 0.2 * cqr_x2 +
  rnorm(n, sd = 0.7)
censoring_point <- 0
cqr_data <- transform(
  test_data,
  cqr_x2 = cqr_x2,
  censored_y = pmax(censoring_point, latent_cqr)
)
missing_censoring <- try(counterfactual_decompose(
  censored_y ~ x1 + cqr_x2, cqr_data, "group", "weights",
  model = "cqr", solver = "fn", control = control
), silent = TRUE)
stopifnot(inherits(missing_censoring, "try-error"))
cqr_control <- cf_control(
  nreg = 5L, trimming = 0.2,
  reported_quantiles = c(0.35, 0.5, 0.65),
  cqr_nsteps = 3L,
  quantile_noncrossing = "rearrange"
)
cqr_fit <- suppressWarnings(counterfactual_decompose(
  censored_y ~ x1 + cqr_x2,
  cqr_data,
  group = "group",
  weights = "weights",
  model = "cqr",
  solver = "fn",
  control = cqr_control,
  censoring = censoring_point
))
stopifnot(all(is.finite(cqr_fit$point$effects)))
stopifnot(max(abs(cqr_fit$point$identity_residual)) < 1e-10)
stopifnot(identical(cqr_fit$metadata$solver_group0, "fn"))
stopifnot(all(
  cqr_fit$point$crossing_diagnostics$crossing_pairs[
    cqr_fit$point$crossing_diagnostics$stage == "corrected"
  ] == 0
))
cqr_data$censoring_column <- censoring_point
prepared_cqr_column <- scalableCounterfactual:::prepare_cf_data(
  censored_y ~ x1 + cqr_x2, cqr_data, "group", "weights",
  model = "cqr", censoring = "censoring_column"
)
stopifnot(all(prepared_cqr_column$censoring0 == censoring_point))

right_limit <- 2
cqr_data$right_censored_y <- pmin(right_limit, latent_cqr)
right_cqr_fit <- suppressWarnings(counterfactual_decompose(
  right_censored_y ~ x1 + cqr_x2,
  cqr_data,
  group = "group",
  weights = "weights",
  model = "cqr",
  solver = "fn",
  control = cf_control(
    nreg = 5L, trimming = 0.2,
    reported_quantiles = c(0.4, 0.5, 0.6),
    cqr_right = TRUE
  ),
  censoring = right_limit
))
stopifnot(all(is.finite(right_cqr_fit$point$effects)))

set.seed(4242)
cox_rate <- exp(-0.25 * group + 0.35 * x1 - 0.15 * x2)
event_time <- rexp(n, rate = cox_rate)
censor_time <- rexp(n, rate = 0.18)
cox_data <- transform(
  test_data,
  duration = pmin(event_time, censor_time),
  event = as.integer(event_time <= censor_time)
)
cox_fit <- counterfactual_decompose(
  duration ~ x1 + x2,
  cox_data,
  group = "group",
  weights = "weights",
  model = "cox",
  event = "event",
  control = cf_control(reported_quantiles = c(0.25, 0.5, 0.75))
)
stopifnot(all(is.finite(cox_fit$point$effects)))
stopifnot(max(abs(cox_fit$point$identity_residual)) < 1e-10)
stopifnot(identical(
  cox_fit$metadata$conditional_backend_group0,
  "survival::coxph.fit(breslow)"
))
cox_group0 <- cox_data$group == 0L
cox_reference_data <- cox_data[cox_group0, ]
cox_reference_data$reference_weight <-
  cox_reference_data$weights / sum(cox_reference_data$weights)
cox_reference <- survival::coxph(
  survival::Surv(duration, event) ~ x1 + x2,
  data = cox_reference_data,
  weights = reference_weight,
  ties = "breslow"
)
stopifnot(max(abs(
  unname(stats::coef(cox_reference)) -
    cox_fit$point$fits$group0$coefficients[-1L]
)) < 1e-7)
stopifnot(cox_fit$point$fits$group0$convergence_flag == 0L)
stopifnot(isTRUE(cox_fit$point$fits$group0$converged))
km_reference <- scalableCounterfactual:::weighted_km_quantiles(
  cox_reference_data$duration,
  cox_reference_data$event,
  cox_reference_data$reference_weight,
  c(0.25, 0.5, 0.75),
  "error"
)
stopifnot(max(abs(
  km_reference - cox_fit$point$diagnostics["reference_observed", ]
)) < 1e-12)

strong_censor_time <- 1:10
strong_censor_event <- c(1L, 1L, rep(0L, 8L))
strong_censor_weight <- rep(1, 10L)
km_na <- scalableCounterfactual:::weighted_km_quantiles(
  strong_censor_time, strong_censor_event, strong_censor_weight,
  c(0.1, 0.5, 0.9), "na"
)
stopifnot(is.finite(km_na[[1L]]), all(is.na(km_na[2:3])))
km_cap <- scalableCounterfactual:::weighted_km_quantiles(
  strong_censor_time, strong_censor_event, strong_censor_weight,
  c(0.5, 0.9), "cap"
)
stopifnot(identical(as.numeric(km_cap), c(2, 2)))
km_error <- try(scalableCounterfactual:::weighted_km_quantiles(
  strong_censor_time, strong_censor_event, strong_censor_weight,
  0.9, "error"
), silent = TRUE)
stopifnot(inherits(km_error, "try-error"))

boundary_inference <- scalableCounterfactual:::effect_inference(
  "total", c(0.2, NA_real_),
  matrix(c(0.18, 0.22, 0.19, NA, NA, NA), nrow = 3L),
  cf_control(reported_quantiles = c(0.5, 0.9)), "cox", NA_character_
)
stopifnot(
  identical(boundary_inference$identified, c(TRUE, FALSE)),
  is.na(boundary_inference$std_error[[2L]]),
  boundary_inference$bootstrap_reps_effective[[2L]] == 0L
)

cox_cuda_error <- try(counterfactual_decompose(
  duration ~ x1 + x2, cox_data, "group", "weights",
  model = "cox", event = "event",
  control = cf_control(gpu_backend = "cuda")
), silent = TRUE)
stopifnot(inherits(cox_cuda_error, "try-error"))
loc_cuda_dr_error <- try(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights", model = "loc",
  control = cf_control(dr_backend = "cuda")
), silent = TRUE)
stopifnot(inherits(loc_cuda_dr_error, "try-error"))

cqr_bootstrap_dir <- tempfile("cf_cqr_bootstrap_")
cqr_bootstrap <- suppressWarnings(counterfactual_decompose(
  censored_y ~ x1 + cqr_x2,
  cqr_data,
  group = "group",
  weights = "weights",
  model = "cqr",
  solver = "fn",
  control = cf_control(
    nreg = 5L, trimming = 0.2,
    reported_quantiles = c(0.4, 0.5, 0.6),
    bootstrap_progress = FALSE
  ),
  censoring = censoring_point,
  bootstrap_reps = 2L,
  checkpoint_dir = cqr_bootstrap_dir,
  seed = 404L
))
stopifnot(all(is.finite(cqr_bootstrap$results$std_error)))
unlink(cqr_bootstrap_dir, recursive = TRUE)

cox_bootstrap_dir <- tempfile("cf_cox_bootstrap_")
cox_bootstrap <- counterfactual_decompose(
  duration ~ x1 + x2,
  cox_data,
  group = "group",
  weights = "weights",
  model = "cox",
  event = "event",
  control = cf_control(
    reported_quantiles = c(0.3, 0.5, 0.7),
    bootstrap_progress = FALSE
  ),
  bootstrap_reps = 2L,
  checkpoint_dir = cox_bootstrap_dir,
  seed = 405L
)
stopifnot(all(is.finite(cox_bootstrap$results$std_error)))
unlink(cox_bootstrap_dir, recursive = TRUE)

linear_reference <- counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "loc", control = cf_control(linear_backend = "qr")
)
linear_chol <- counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "loc", control = cf_control(linear_backend = "chol")
)
stopifnot(max(abs(
  linear_reference$point$effects - linear_chol$point$effects
)) < 1e-10)
if (requireNamespace("fastglm", quietly = TRUE)) {
  linear_fastglm <- counterfactual_decompose(
    y ~ x1 + x2, test_data, "group", "weights",
    model = "loc", control = cf_control(linear_backend = "fastglm")
  )
  stopifnot(max(abs(
    linear_reference$point$effects - linear_fastglm$point$effects
  )) < 1e-7)
}

dr_reference <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "logit",
  control = cf_control(nreg = 9L, dr_backend = "glm", dr_warm_start = FALSE)
))
fallback_X <- stats::model.matrix(~ x1 + x2, test_data)
fallback_threshold <- stats::median(test_data$y)
fallback_process <- list(
  coefficients = matrix(0, ncol(fallback_X), 1L),
  converged = FALSE,
  boundary = FALSE,
  iterations = 1L,
  backend = "cuda",
  fallback_reason = NA_character_
)
fallback_result <- scalableCounterfactual:::refit_cuda_dr_failures(
  fallback_process, fallback_X, test_data$y, test_data$weights,
  fallback_threshold, "logit", 100L, 1e-8, c(1, 0, 0)
)
stopifnot(
  identical(fallback_result$fallback, 1L),
  identical(fallback_result$process$backend, "glm"),
  identical(
    fallback_result$process$fallback_reason, "cuda_nonconvergence"
  ),
  isTRUE(fallback_result$process$converged)
)
dr_warm <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "logit",
  control = cf_control(nreg = 9L, dr_backend = "glm", dr_warm_start = TRUE)
))
stopifnot(max(abs(dr_reference$point$effects - dr_warm$point$effects)) < 1e-6)
stopifnot(isTRUE(dr_warm$metadata$dr_warm_start_effective_group0))
stopifnot(isTRUE(dr_warm$metadata$dr_warm_start_effective_group1))
for (group_name in c("group0", "group1")) {
  group_design <- if (group_name == "group0") {
    scalableCounterfactual:::prepare_cf_data(
      y ~ x1 + x2, test_data, "group", "weights"
    )$X0
  } else {
    scalableCounterfactual:::prepare_cf_data(
      y ~ x1 + x2, test_data, "group", "weights"
    )$X1
  }
  endpoint_fit <- dr_warm$point$fits[[group_name]]
  endpoint_probability <- stats::plogis(drop(
    group_design %*% endpoint_fit$coefficients[, ncol(endpoint_fit$coefficients)]
  ))
  stopifnot(max(abs(endpoint_probability - (1 - 1e-8))) < 1e-10)
}

dr_parallel <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "logit",
  control = cf_control(
    nreg = 9L, dr_backend = "glm", dr_workers = 2L,
    dr_warm_start = FALSE
  )
))
stopifnot(max(abs(
  dr_reference$point$effects - dr_parallel$point$effects
)) < 1e-10)

if (requireNamespace("fastglm", quietly = TRUE)) {
  dr_fastglm <- suppressWarnings(counterfactual_decompose(
    y ~ x1 + x2, test_data, "group", "weights",
    model = "logit",
    control = cf_control(
      nreg = 9L, dr_backend = "fastglm", dr_warm_start = FALSE
    )
  ))
  stopifnot(max(abs(
    dr_reference$point$effects - dr_fastglm$point$effects
  )) < 1e-5)
  stopifnot(identical(
    dr_fastglm$metadata$conditional_backend_group0, "fastglm"
  ))

  dr_bootstrap_dir <- tempfile("cf_dr_fastglm_")
  dr_fastglm_bootstrap <- suppressWarnings(counterfactual_decompose(
    y ~ x1 + x2, test_data, "group", "weights",
    model = "logit",
    control = cf_control(
      nreg = 7L, dr_backend = "fastglm", dr_warm_start = TRUE
    ),
    bootstrap_reps = 2L,
    bootstrap_workers = 1L,
    checkpoint_dir = dr_bootstrap_dir,
    seed = 808L
  ))
  stopifnot(all(is.finite(dr_fastglm_bootstrap$results$std_error)))
  unlink(dr_bootstrap_dir, recursive = TRUE)
}
if (requireNamespace("speedglm", quietly = TRUE)) {
  dr_speedglm <- suppressWarnings(counterfactual_decompose(
    y ~ x1 + x2, test_data, "group", "weights",
    model = "logit",
    control = cf_control(
      nreg = 9L, dr_backend = "speedglm", dr_warm_start = FALSE
    )
  ))
  stopifnot(max(abs(
    dr_reference$point$effects - dr_speedglm$point$effects
  )) < 1e-5)
}

parallel_layer_error <- try(counterfactual_decompose(
  y ~ x1 + x2, test_data, "group", "weights",
  model = "logit",
  control = cf_control(nreg = 9L, dr_workers = 2L),
  point_workers = 2L
), silent = TRUE)
stopifnot(inherits(parallel_layer_error, "try-error"))
stopifnot(all(c(
  "model", "backend", "implementation", "objective_preserving"
) %in% names(conditional_backend_registry())))

constant_data <- transform(test_data, constant_covariate = 1)
constant_fit <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2 + constant_covariate,
  constant_data,
  group = "group",
  weights = "weights",
  model = "loc",
  control = control
))
stopifnot(identical(
  constant_fit$metadata$dropped_design_columns,
  "constant_covariate"
))

factor_error <- try(counterfactual_decompose(
  factor(x2) ~ x1,
  test_data,
  group = "group",
  weights = "weights",
  model = "loc",
  control = control
), silent = TRUE)
stopifnot(inherits(factor_error, "try-error"))

checkpoint_dir <- tempfile("cf_bootstrap_")
boot <- counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  model = "loc",
  control = control,
  bootstrap_reps = 2L,
  bootstrap_workers = if (parallel::detectCores() >= 2L) 2L else 1L,
  checkpoint_dir = checkpoint_dir,
  seed = 99L
)
stopifnot(all(is.finite(boot$results$std_error)))
stopifnot(all(boot$bootstrap$resources$attempt == 0L))
stopifnot(nrow(boot$bootstrap$failures) == 0L)
checkpoint_files <- sort(list.files(
  boot$bootstrap$checkpoint_dir,
  pattern = "\\.rds$",
  full.names = TRUE
))
checkpoint_hashes <- unname(tools::md5sum(checkpoint_files))
boot_extended <- counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  model = "loc",
  control = control,
  bootstrap_reps = 3L,
  bootstrap_workers = if (parallel::detectCores() >= 2L) 2L else 1L,
  checkpoint_dir = checkpoint_dir,
  seed = 99L
)
stopifnot(identical(boot$bootstrap$signature, boot_extended$bootstrap$signature))
stopifnot(identical(
  normalizePath(boot_extended$metadata$checkpoint_dir, winslash = "/"),
  normalizePath(boot_extended$bootstrap$checkpoint_dir, winslash = "/")
))
stopifnot(length(list.files(
  boot_extended$bootstrap$checkpoint_dir,
  pattern = "\\.rds$"
)) == 3L)
stopifnot(identical(
  checkpoint_hashes,
  unname(tools::md5sum(checkpoint_files))
))

changed_data <- transform(test_data, y = y + seq_len(nrow(test_data)) * 1e-8)
boot_changed <- counterfactual_decompose(
  y ~ x1 + x2,
  changed_data,
  group = "group",
  weights = "weights",
  model = "loc",
  control = control,
  bootstrap_reps = 1L,
  bootstrap_workers = 1L,
  checkpoint_dir = checkpoint_dir,
  seed = 99L
)
stopifnot(!identical(
  boot_extended$bootstrap$signature,
  boot_changed$bootstrap$signature
))
unlink(checkpoint_dir, recursive = TRUE)

retry_dir <- tempfile("cf_retry_")
dir.create(retry_dir)
retry_namespace <- asNamespace("scalableCounterfactual")
retry_name <- "bootstrap_replication_attempt"
retry_original <- get(retry_name, envir = retry_namespace)
unlockBinding(retry_name, retry_namespace)
assign(retry_name, function(rep_id, attempt, common) {
  if (attempt == 0L) stop("intentional first-attempt failure")
  list(
    signature = common$signature,
    data_fingerprint = "test",
    replication = as.integer(rep_id),
    attempt = as.integer(attempt),
    seed = scalableCounterfactual:::bootstrap_attempt_seed(
      common$seed, rep_id, attempt
    ),
    attempt_failures = data.frame()
  )
}, envir = retry_namespace)
lockBinding(retry_name, retry_namespace)
retry_result <- tryCatch(
  scalableCounterfactual:::bootstrap_replication(1L, list(
    run_dir = retry_dir,
    signature = "retry_test",
    seed = 11L,
    control = cf_control(bootstrap_max_retries = 1L)
  )),
  finally = {
    unlockBinding(retry_name, retry_namespace)
    assign(retry_name, retry_original, envir = retry_namespace)
    lockBinding(retry_name, retry_namespace)
  }
)
stopifnot(identical(retry_result$status, "ok"))
stopifnot(retry_result$attempt == 1L)
stopifnot(nrow(retry_result$failures) == 1L)
unlink(retry_dir, recursive = TRUE)

xy_standard_dir <- tempfile("cf_xy_standard_")
xy_preprocess_dir <- tempfile("cf_xy_preprocess_")
xy_standard <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  solver = "profn",
  control = cf_control(
    nreg = 9L,
    trimming = 0.05,
    reported_quantiles = c(0.25, 0.5, 0.75),
    qr_bootstrap_engine = "standard"
  ),
  bootstrap_reps = 1L,
  checkpoint_dir = xy_standard_dir,
  seed = 123L
))
xy_preprocess <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  solver = "profn",
  control = cf_control(
    nreg = 9L,
    trimming = 0.05,
    reported_quantiles = c(0.25, 0.5, 0.75),
    qr_bootstrap_engine = "xy_preprocess"
  ),
  bootstrap_reps = 1L,
  checkpoint_dir = xy_preprocess_dir,
  seed = 123L
))
stopifnot(identical(xy_preprocess$bootstrap$engine, "xy_preprocess"))
stopifnot(max(abs(
  xy_standard$bootstrap$effects$total -
    xy_preprocess$bootstrap$effects$total
)) < 1e-5)
unlink(c(xy_standard_dir, xy_preprocess_dir), recursive = TRUE)

multiplier_xy_error <- try(counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  solver = "profn",
  control = cf_control(
    nreg = 9L,
    trimming = 0.05,
    bootstrap_scheme = "multiplier",
    qr_bootstrap_engine = "xy_preprocess"
  ),
  bootstrap_reps = 1L,
  checkpoint_dir = tempfile("cf_xy_multiplier_")
), silent = TRUE)
stopifnot(inherits(multiplier_xy_error, "try-error"))

imbalanced <- data.frame(
  y = rnorm(1000L),
  x = rnorm(1000L),
  group = c(rep(0L, 5L), rep(1L, 995L)),
  weight = 1
)
imbalanced_prepared <- scalableCounterfactual:::prepare_cf_data(
  y ~ x, imbalanced, "group", "weight"
)
imbalanced_sample <- scalableCounterfactual:::subsample_prepared_data(
  imbalanced_prepared, 50L, 101L
)
stopifnot(imbalanced_sample$n == 50L)
stopifnot(imbalanced_sample$n0 > ncol(imbalanced_sample$X0))
stopifnot(imbalanced_sample$n1 > ncol(imbalanced_sample$X1))

separation_block <- data.frame(
  y = c(seq(-2, -0.1, length.out = 30), seq(0.1, 2, length.out = 30)),
  x = c(rep(-1, 30), rep(1, 30))
)
separated <- rbind(
  transform(separation_block, group = 0L, weight = 1),
  transform(separation_block, group = 1L, weight = 1)
)
separation_fit <- suppressWarnings(counterfactual_decompose(
  y ~ x, separated, "group", "weight",
  model = "logit",
  control = cf_control(
    nreg = 5L, reported_quantiles = 0.5, dr_backend = "glm"
  )
))
stopifnot(nchar(separation_fit$metadata$conditional_model_warnings) > 0L)
stopifnot(any(
  separation_fit$point$fits$group0$convergence_flag == 2L |
    separation_fit$point$fits$group1$convergence_flag == 2L
))

onestep_bootstrap_dir <- tempfile("cf_onestep_")
onestep_bootstrap <- suppressWarnings(counterfactual_decompose(
  y ~ x1 + x2,
  test_data,
  group = "group",
  weights = "weights",
  solver = "onestep",
  control = cf_control(
    nreg = 21L,
    trimming = 0.02,
    reported_quantiles = c(0.25, 0.5, 0.75),
    qr_bootstrap_engine = "auto"
  ),
  bootstrap_reps = 2L,
  bootstrap_workers = 1L,
  checkpoint_dir = onestep_bootstrap_dir,
  seed = 321L
))
stopifnot(identical(onestep_bootstrap$bootstrap$engine, "onestep"))
stopifnot(all(is.finite(onestep_bootstrap$results$std_error)))
stopifnot(all(onestep_bootstrap$bootstrap$resources$elapsed_seconds > 0))
functional_tests <- functional_effect_tests(
  onestep_bootstrap, constants = 0.1, quantile_range = c(0.25, 0.75)
)
stopifnot(nrow(functional_tests) == 15L)
stopifnot(all(functional_tests$effect %in% c(
  "structure", "composition", "total"
)))
stopifnot(all(functional_tests$ks_p_value >= 0 &
  functional_tests$ks_p_value <= 1))
stopifnot(all(functional_tests$cms_p_value >= 0 &
  functional_tests$cms_p_value <= 1))
functional_output <- tempfile("cf_functional_output_")
write_cf_outputs(onestep_bootstrap, functional_output)
stopifnot(file.exists(file.path(
  functional_output, "functional_effect_tests.csv"
)))
unlink(functional_output, recursive = TRUE)
unlink(onestep_bootstrap_dir, recursive = TRUE)
