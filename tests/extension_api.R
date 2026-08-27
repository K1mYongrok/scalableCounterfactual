library(scalableCounterfactual)

set.seed(91)
n <- 300L
data <- data.frame(
  y = rnorm(n),
  x = rnorm(n),
  group = rep(0:1, each = n / 2L),
  weight = runif(n, 0.5, 1.5)
)

custom_fn <- function(X, y, weights, taus, control) {
  weighted_X <- X * weights
  weighted_y <- y * weights
  coefficients <- vapply(taus, function(tau) {
    as.numeric(quantreg::rq.fit.fnb(weighted_X, weighted_y, tau)$coefficients)
  }, numeric(ncol(X)))
  list(
    coefficients = coefficients,
    flag = rep(0L, length(taus)),
    warnings = character()
  )
}
register_qr_solver(
  "testfn", custom_fn, exact = TRUE,
  description = "test wrapper around rq.fit.fnb"
)
stopifnot("testfn" %in% qr_solver_registry()$solver)
reference <- counterfactual_decompose(
  y ~ x, data, "group", "weight", solver = "fn",
  control = cf_control(
    nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
    crossing_diagnostics = FALSE
  )
)
checkpoint <- tempfile("cf_custom_solver_")
custom <- counterfactual_decompose(
  y ~ x, data, "group", "weight", solver = "testfn",
  control = cf_control(
    nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
    crossing_diagnostics = FALSE, bootstrap_progress = FALSE
  ),
  bootstrap_reps = 2L,
  point_workers = 2L,
  bootstrap_workers = 2L,
  checkpoint_dir = checkpoint,
  seed = 12L
)
stopifnot(max(abs(custom$point$effects - reference$point$effects)) < 1e-7)
stopifnot(all(is.finite(custom$bootstrap$effects$total)))
stopifnot(identical(custom$metadata$bootstrap_point_workers, 1L))
unlink(checkpoint, recursive = TRUE)
stopifnot(unregister_qr_solver("testfn"))

solver_factory <- function(multiplier) {
  force(multiplier)
  function(X, y, weights, taus, control) {
    coefficients <- vapply(taus, function(tau) {
      as.numeric(quantreg::rq.fit.fnb(
        X * weights, (y * multiplier) * weights, tau
      )$coefficients)
    }, numeric(ncol(X)))
    list(coefficients = coefficients)
  }
}
register_qr_solver("closuretest", solver_factory(1), exact = TRUE)
fingerprint_one <- scalableCounterfactual:::custom_registry_entry(
  "qr", "closuretest"
)$function_fingerprint
register_qr_solver(
  "closuretest", solver_factory(2), exact = TRUE, overwrite = TRUE
)
fingerprint_two <- scalableCounterfactual:::custom_registry_entry(
  "qr", "closuretest"
)$function_fingerprint
stopifnot(!identical(fingerprint_one, fingerprint_two))
stopifnot(unregister_qr_solver("closuretest"))

global_multiplier <- 1
global_solver <- function(X, y, weights, taus, control) {
  coefficients <- vapply(taus, function(tau) {
    as.numeric(quantreg::rq.fit.fnb(
      X * weights, (y * global_multiplier) * weights, tau
    )$coefficients)
  }, numeric(ncol(X)))
  list(coefficients = coefficients)
}
register_qr_solver("globalsafe", global_solver, exact = TRUE)
global_parallel <- counterfactual_decompose(
  y ~ x, data, "group", "weight", solver = "globalsafe",
  point_workers = 2L,
  control = cf_control(
    nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
    crossing_diagnostics = FALSE
  )
)
stopifnot(max(abs(
  global_parallel$point$effects - reference$point$effects
)) < 1e-7)
stopifnot(unregister_qr_solver("globalsafe"))

signature_control <- cf_control(
  nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
  marginal_method = "matrix", crossing_diagnostics = FALSE
)
legacy_signature_control <- cf_control(
  nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
  marginal_method = "matrix", legacy_weighted_quantile = TRUE,
  crossing_diagnostics = FALSE
)
prepared <- scalableCounterfactual:::prepare_cf_data(
  y ~ x, data, "group", "weight"
)
signature_default <- scalableCounterfactual:::bootstrap_signature(
  prepared, "qr", "fn", signature_control, 12L, "standard",
  point = reference$point
)
signature_legacy <- scalableCounterfactual:::bootstrap_signature(
  prepared, "qr", "fn", legacy_signature_control, 12L, "standard",
  point = reference$point
)
stopifnot(!identical(signature_default, signature_legacy))

active_before <- scalableCounterfactual:::active_extension_fingerprint(
  "qr", "fn", signature_control, reference$point
)
register_qr_solver("unusedsolver", custom_fn, exact = TRUE)
active_after <- scalableCounterfactual:::active_extension_fingerprint(
  "qr", "fn", signature_control, reference$point
)
stopifnot(identical(active_before, active_after))
stopifnot(unregister_qr_solver("unusedsolver"))

custom_lm <- function(X, y, weights) {
  fit <- stats::lm.wfit(X, y, weights)
  list(coefficients = fit$coefficients, residuals = fit$residuals)
}
register_conditional_backend(
  "testlm", "linear", custom_lm,
  description = "test lm.wfit backend"
)
linear_reference <- counterfactual_decompose(
  y ~ x, data, "group", "weight", model = "loc",
  control = cf_control(
    linear_backend = "qr", reported_quantiles = c(0.25, 0.5, 0.75)
  )
)
linear_custom <- counterfactual_decompose(
  y ~ x, data, "group", "weight", model = "loc",
  control = cf_control(
    linear_backend = "testlm", reported_quantiles = c(0.25, 0.5, 0.75)
  )
)
stopifnot(max(abs(
  linear_reference$point$effects - linear_custom$point$effects
)) < 1e-8)
stopifnot(unregister_conditional_backend("testlm", "linear"))

custom_glm <- function(
    X, response, weights, model, start, maxit, tolerance) {
  fit <- stats::glm.fit(
    X, response, weights = weights,
    family = stats::binomial(link = model),
    intercept = FALSE, start = start,
    control = stats::glm.control(epsilon = tolerance, maxit = maxit)
  )
  list(
    coefficients = fit$coefficients,
    fitted.values = fit$fitted.values,
    converged = fit$converged,
    boundary = fit$boundary,
    iterations = fit$iter
  )
}
register_conditional_backend(
  "testglm", "distribution", custom_glm,
  description = "test glm.fit backend"
)
dr_reference <- suppressWarnings(counterfactual_decompose(
  y ~ x, data, "group", "weight", model = "logit",
  control = cf_control(
    nreg = 7L, dr_backend = "glm", dr_warm_start = FALSE,
    reported_quantiles = c(0.25, 0.5, 0.75)
  )
))
dr_custom <- suppressWarnings(counterfactual_decompose(
  y ~ x, data, "group", "weight", model = "logit",
  control = cf_control(
    nreg = 7L, dr_backend = "testglm", dr_warm_start = FALSE,
    reported_quantiles = c(0.25, 0.5, 0.75)
  )
))
stopifnot(max(abs(dr_reference$point$effects - dr_custom$point$effects)) < 1e-7)
stopifnot(unregister_conditional_backend("testglm", "distribution"))

custom_glm_without_diagnostics <- function(
    X, response, weights, model, start, maxit, tolerance) {
  fit <- stats::glm.fit(
    X, response, weights = weights,
    family = stats::binomial(link = model), intercept = FALSE,
    control = stats::glm.control(epsilon = tolerance, maxit = maxit)
  )
  list(coefficients = fit$coefficients)
}
register_conditional_backend(
  "testglmnodiag", "distribution", custom_glm_without_diagnostics
)
dr_without_diagnostics <- suppressWarnings(counterfactual_decompose(
  y ~ x, data, "group", "weight", model = "logit",
  control = cf_control(
    nreg = 7L, dr_backend = "testglmnodiag", dr_warm_start = FALSE,
    reported_quantiles = c(0.25, 0.5, 0.75)
  )
))
stopifnot(anyNA(
  dr_without_diagnostics$point$fits$group0$convergence_flag
))
stopifnot(unregister_conditional_backend("testglmnodiag", "distribution"))

invalid_probability_backend <- function(
    X, response, weights, model, start, maxit, tolerance) {
  list(
    coefficients = rep(0, ncol(X)), fitted.values = 2,
    converged = TRUE, boundary = FALSE, iterations = 1L
  )
}
register_conditional_backend(
  "invalidprob", "distribution", invalid_probability_backend
)
invalid_probability_result <- tryCatch(
  counterfactual_decompose(
    y ~ x, data, "group", "weight", model = "logit",
    control = cf_control(
      nreg = 7L, dr_backend = "invalidprob", dr_warm_start = FALSE,
      reported_quantiles = c(0.25, 0.5, 0.75)
    )
  ),
  error = conditionMessage
)
stopifnot(is.character(invalid_probability_result))
stopifnot(grepl("invalid fitted.values", invalid_probability_result, fixed = TRUE))
stopifnot(unregister_conditional_backend("invalidprob", "distribution"))

summarized <- summary(custom, quantiles = c(0.25, 0.75))
stopifnot(inherits(summarized, "summary.cfdecomp"))
stopifnot(nrow(summarized$effects) == 6L)
plot_path <- tempfile("cf_plot_", fileext = ".pdf")
grDevices::pdf(plot_path)
plotted <- plot(custom, effects = c("structure", "total"), interval = "uniform")
grDevices::dev.off()
stopifnot(length(plotted) == 2L, file.exists(plot_path))
unlink(plot_path)

simulation <- simulate_counterfactual_validation(
  replications = 2L,
  n_per_group = 120L,
  solver = "fn",
  nreg = 9L,
  reported_quantiles = c(0.25, 0.5, 0.75),
  workers = 2L,
  progress = FALSE,
  seed = 33L
)
stopifnot(inherits(simulation, "cf_simulation_validation"))
stopifnot(simulation$resources$failed_tasks == 0L)
stopifnot(nrow(simulation$summary) == 27L)
stopifnot(simulation$resources$maximum_identity_residual < 1e-10)

simulation_inference <- simulate_counterfactual_validation(
  replications = 2L,
  n_per_group = 80L,
  scenarios = "location_shift",
  solver = "fn",
  nreg = 7L,
  reported_quantiles = c(0.25, 0.5, 0.75),
  bootstrap_reps = 2L,
  workers = 1L,
  progress = FALSE,
  seed = 34L
)
stopifnot(all(!is.na(simulation_inference$raw$pointwise_covered)))
stopifnot(all(!is.na(simulation_inference$raw$uniform_curve_covered)))
stopifnot(all(!is.na(simulation_inference$curve_coverage$uniform_coverage)))

simulation_checkpoint <- tempfile("cf_simulation_checkpoint_")
simulation_cached_first <- simulate_counterfactual_validation(
  replications = 2L, n_per_group = 80L, scenarios = "location_shift",
  solver = "fn", nreg = 7L, reported_quantiles = c(0.25, 0.5, 0.75),
  workers = 2L, progress = FALSE, checkpoint_dir = simulation_checkpoint,
  seed = 35L
)
simulation_cached_second <- simulate_counterfactual_validation(
  replications = 2L, n_per_group = 80L, scenarios = "location_shift",
  solver = "fn", nreg = 7L, reported_quantiles = c(0.25, 0.5, 0.75),
  workers = 2L, progress = FALSE, checkpoint_dir = simulation_checkpoint,
  seed = 35L
)
stopifnot(simulation_cached_first$resources$cached_tasks == 0L)
stopifnot(simulation_cached_second$resources$cached_tasks == 2L)
stopifnot(identical(simulation_cached_first$raw, simulation_cached_second$raw))
unlink(simulation_checkpoint, recursive = TRUE)

large_seed_simulation <- simulate_counterfactual_validation(
  replications = 1L, n_per_group = 30L, scenarios = "location_shift",
  solver = "fn", nreg = 5L, reported_quantiles = 0.5,
  workers = 1L, progress = FALSE, seed = .Machine$integer.max
)
stopifnot(large_seed_simulation$resources$successful_tasks == 1L)

simulation_output <- tempfile("cf_simulation_output_")
dir.create(simulation_output)
writeLines("stale", file.path(simulation_output, "failures.csv"))
write_simulation_validation(simulation, simulation_output)
stopifnot(!file.exists(file.path(simulation_output, "failures.csv")))
stopifnot(file.exists(file.path(simulation_output, "curve_coverage.csv")))
unlink(simulation_output, recursive = TRUE)

decomposition_output <- tempfile("cf_decomposition_output_")
write_cf_outputs(custom, decomposition_output)
stopifnot(file.exists(file.path(decomposition_output, "bootstrap_resources.csv")))
write_cf_outputs(reference, decomposition_output)
stopifnot(!file.exists(file.path(decomposition_output, "bootstrap_resources.csv")))
stopifnot(!file.exists(file.path(decomposition_output, "bootstrap_failures.csv")))
unlink(decomposition_output, recursive = TRUE)
