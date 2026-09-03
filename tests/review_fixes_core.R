library(scalableCounterfactual)

default_control <- cf_control()
stopifnot(
  identical(default_control$legacy_qr_shift, TRUE),
  identical(default_control$dr_noncrossing, "cummax")
)
stopifnot(identical(
  cf_control(dr_noncrossing = "isotonic")$dr_noncrossing,
  "isotonic"
))

offset_data <- data.frame(
  y = rnorm(20), x = rnorm(20), off = rnorm(20), group = rep(0:1, each = 10)
)
offset_error <- try(
  scalableCounterfactual:::prepare_cf_data(
    y ~ x + offset(off), offset_data, "group"
  ),
  silent = TRUE
)
stopifnot(
  inherits(offset_error, "try-error"),
  grepl("offset\\(\\) terms are not supported", as.character(offset_error))
)

finite_data <- data.frame(
  y = rnorm(30), x = rnorm(30), q = rnorm(30),
  group = rep(0:1, each = 15)
)
finite_data$x[[3L]] <- Inf
finite_prepared <- scalableCounterfactual:::prepare_cf_data(
  y ~ x + q, finite_data, "group"
)
stopifnot(finite_prepared$n == 29L, finite_prepared$omitted_rows == 1L)

# An intercept-only specification is identified and has a zero composition
# effect. It is useful for unconditional two-group distribution comparisons.
set.seed(4401)
intercept_data <- data.frame(
  y = c(rnorm(50), 0.3 + rnorm(50)),
  group = rep(0:1, each = 50),
  weight = runif(100, 0.5, 2)
)
for (intercept_model in c(
    "qr", "loc", "locsca", "logit", "probit", "cloglog", "lpm")) {
  intercept_fit <- suppressWarnings(counterfactual_decompose(
    y ~ 1, intercept_data, "group", "weight", model = intercept_model,
    solver = if (intercept_model == "qr") "fn" else NULL,
    control = cf_control(
      nreg = 5L, reported_quantiles = 0.5,
      crossing_diagnostics = FALSE
    )
  ))
  stopifnot(
    is.finite(intercept_fit$point$effects["composition", 1L]),
    abs(intercept_fit$point$effects["composition", 1L]) < 1e-12
  )
}
intercept_cqr_data <- intercept_data
set.seed(4402)
intercept_cqr_data$y <- c(
  rep(0, 2L), exp(rnorm(48L)),
  rep(0, 2L), exp(0.2 + rnorm(48L))
)
intercept_cqr <- suppressWarnings(counterfactual_decompose(
  y ~ 1, intercept_cqr_data, "group", "weight", model = "cqr",
  solver = "fn", censoring = 0,
  control = cf_control(
    nreg = 5L, reported_quantiles = 0.5,
    crossing_diagnostics = FALSE
  )
))
stopifnot(abs(intercept_cqr$point$effects["composition", 1L]) < 1e-12)

set.seed(44)
X0 <- cbind(
  "(Intercept)" = 1, tiny = 1e-12 * rnorm(30), q = rnorm(30)
)
X1 <- cbind(
  "(Intercept)" = 1, tiny = 1e-12 * rnorm(30), q = rnorm(30)
)
scale_reduced <- scalableCounterfactual:::reduce_common_design(X0, X1)
stopifnot(identical(
  scale_reduced$retained, c("(Intercept)", "tiny", "q")
))

a0 <- 1:8
c0 <- c(-2, 1, 0, 3, -1, 4, 2, 5)
b1 <- c(2, -1, 3, 0, 5, 1, 4, -2)
a1 <- 2:9
X0 <- cbind("(Intercept)" = 1, a = a0, b = a0, c = c0)
X1 <- cbind("(Intercept)" = 1, a = a1, b = b1, c = a1)
joint_reduced <- scalableCounterfactual:::reduce_common_design(X0, X1)
stopifnot(
  identical(joint_reduced$retained, c("(Intercept)", "b", "c")),
  qr(joint_reduced$X0)$rank == 3L,
  qr(joint_reduced$X1)$rank == 3L
)

no_intercept_reduced <- scalableCounterfactual:::reduce_common_design(
  X0[, -1L, drop = FALSE], X1[, -1L, drop = FALSE]
)
stopifnot(
  ncol(no_intercept_reduced$X0) == 2L,
  qr(no_intercept_reduced$X0)$rank == 2L,
  qr(no_intercept_reduced$X1)$rank == 2L
)

# Frequency-compressed resamples can have exactly p stored, full-rank rows
# while still representing more than p sampled observations.
frequency_design0 <- cbind("(Intercept)" = 1, x = c(-1, 1))
frequency_design1 <- cbind("(Intercept)" = 1, x = c(-2, 2))
frequency_reduced <- scalableCounterfactual:::reduce_common_design(
  frequency_design0, frequency_design1,
  effective_n0 = 3, effective_n1 = 4
)
stopifnot(
  identical(ncol(frequency_reduced$X0), 2L),
  inherits(try(
    scalableCounterfactual:::reduce_common_design(
      frequency_design0, frequency_design1
    ), silent = TRUE
  ), "try-error")
)

single_column <- matrix(1:5, ncol = 1L)
single_rearranged <- scalableCounterfactual:::rearrange_quantile_rows(
  single_column
)
stopifnot(
  identical(dim(single_rearranged), c(5L, 1L)),
  identical(as.numeric(single_rearranged), as.numeric(single_column))
)

set.seed(1)
n <- 22L
X <- cbind("(Intercept)" = 1, x = rnorm(n))
location_fit <- structure(
  list(coefficients = c(0, 1), residual_quantiles = c(-1, 0, 1)),
  class = c("cf_loc_fit", "cf_conditional_fit")
)
extreme_weights <- 10^runif(n, -16, 16)
probability <- 0.33
matrix_control <- cf_control(
  reported_quantiles = probability, marginal_method = "matrix"
)
chunked_control <- cf_control(
  reported_quantiles = probability, marginal_method = "chunked",
  marginal_chunk_rows = 5L, marginal_histogram_bins = 256L,
  marginal_candidate_max = 100000L
)
matrix_quantile <- scalableCounterfactual:::marginal_quantiles(
  location_fit, X, extreme_weights, probability, matrix_control
)
chunked_quantile <- scalableCounterfactual:::marginal_quantiles(
  location_fit, X, extreme_weights, probability, chunked_control
)
stopifnot(max(abs(matrix_quantile - chunked_quantile)) < 1e-14)

# Integer bootstrap frequencies retain the expanded sample's type-7 rank
# positions in both matrix and bounded-memory marginalization.
frequency_X <- cbind("(Intercept)" = 1, x = c(0, 10))
frequency_fit <- structure(list(
  model = "loc", coefficients = c(0, 1), residual_quantiles = 0,
  taus = 0.5
), class = c("cf_loc_fit", "cf_conditional_fit"))
frequency_matrix <- scalableCounterfactual:::marginal_quantiles(
  frequency_fit, frequency_X, c(2, 1), 0.5, matrix_control,
  normalization_rows = 3
)
frequency_chunked <- scalableCounterfactual:::marginal_quantiles(
  frequency_fit, frequency_X, c(2, 1), 0.5, chunked_control,
  normalization_rows = 3
)
stopifnot(
  identical(as.numeric(frequency_matrix), 0),
  identical(as.numeric(frequency_chunked), 0)
)

# The package-wide default uses the same stable weighted type-7 rule. Tiny
# positive weights that do not advance a floating-point cumulative sum must
# not make matrix and vector paths disagree.
vector_quantile <- scalableCounterfactual:::weighted_quantile(
  as.vector(scalableCounterfactual:::predict_conditional_draws(
    location_fit, X, matrix_control
  )),
  rep(extreme_weights, times = length(location_fit$residual_quantiles)),
  probability
)
stopifnot(max(abs(matrix_quantile - vector_quantile)) < 1e-14)
equal_weight_values <- c(-2, -1, 0, 4, 9)
stopifnot(isTRUE(all.equal(
  scalableCounterfactual:::weighted_quantile(
    equal_weight_values, rep(1, 5), c(0.25, 0.5, 0.75)
  ),
  as.numeric(stats::quantile(
    equal_weight_values, c(0.25, 0.5, 0.75), type = 7
  ))
)))
moderate_values <- c(-3.2, -1.1, -0.4, 0.2, 0.8, 2.7, 4.1)
moderate_weights <- c(0.7, 1.3, 0.9, 1.8, 0.6, 1.1, 1.5)
moderate_probs <- c(0.1, 0.5, 0.9)
stopifnot(isTRUE(all.equal(
  scalableCounterfactual:::weighted_quantile(
    moderate_values, moderate_weights, moderate_probs
  ),
  as.numeric(Hmisc::wtd.quantile(
    moderate_values,
    weights = moderate_weights,
    probs = moderate_probs,
    normwt = TRUE
  )),
  tolerance = 0
)))

set.seed(8)
X <- cbind("(Intercept)" = 1, matrix(rnorm(300), 100, 3))
dr_fit <- structure(
  list(
    coefficients = matrix(rnorm(4 * 9), 4, 9),
    thresholds = seq_len(9), model = "logit"
  ),
  class = c("cf_dr_fit", "cf_conditional_fit")
)
dr_weights <- runif(100, 0.1, 3)
dr_matrix <- scalableCounterfactual:::cpu_dr_marginal_cdf(
  X, dr_fit, dr_weights,
  cf_control(marginal_method = "matrix")
)
dr_chunked <- scalableCounterfactual:::cpu_dr_marginal_cdf(
  X, dr_fit, dr_weights,
  cf_control(marginal_method = "chunked", marginal_chunk_rows = 13L)
)
dr_auto <- scalableCounterfactual:::cpu_dr_marginal_cdf(
  X, dr_fit, dr_weights,
  cf_control(marginal_method = "auto", marginal_matrix_max_mb = 0.001,
             marginal_chunk_rows = 13L)
)
stopifnot(
  max(abs(dr_matrix$cdf - dr_chunked$cdf)) < 1e-12,
  identical(dr_matrix$method, "cdf_matrix"),
  identical(dr_chunked$method, "cdf_chunked"),
  identical(dr_auto$method, "cdf_chunked"),
  dr_chunked$passes == 8L
)

set.seed(202)
dr_data <- data.frame(
  group = rep(0:1, each = 100L),
  x = rnorm(200L),
  weights = runif(200L, 0.2, 2)
)
dr_data$y <- 0.4 * dr_data$group + 0.6 * dr_data$x + rnorm(200L)
dr_matrix_fit <- suppressWarnings(counterfactual_decompose(
  y ~ x, dr_data, "group", "weights", model = "logit",
  control = cf_control(
    nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
    marginal_method = "matrix", dr_backend = "glm"
  )
))
dr_chunked_fit <- suppressWarnings(counterfactual_decompose(
  y ~ x, dr_data, "group", "weights", model = "logit",
  control = cf_control(
    nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
    marginal_method = "chunked", marginal_chunk_rows = 17L,
    dr_backend = "glm"
  )
))
stopifnot(
  max(abs(
    dr_matrix_fit$point$effects - dr_chunked_fit$point$effects
  )) < 1e-12,
  all(dr_chunked_fit$point$marginal_diagnostics$method ==
        "cdf_chunked_cummax"),
  identical(dr_chunked_fit$metadata$dr_noncrossing, "cummax"),
  is.numeric(dr_chunked_fit$metadata$dr_raw_crossing_pairs_max),
  is.numeric(dr_chunked_fit$metadata$dr_max_noncrossing_adjustment)
)

rearranged_cdf <- scalableCounterfactual:::inverse_step_cdf(
  1:3, c(0.2, 0.8, 0.4), 0.5, "rearrange"
)
default_cdf <- scalableCounterfactual:::inverse_step_cdf(
  1:3, c(0.2, 0.8, 0.4), 0.5
)
uncorrected_cdf <- scalableCounterfactual:::inverse_step_cdf(
  1:3, c(0.2, 0.8, 0.4), 0.5, "none"
)
rearrangement_diagnostics <- attr(
  rearranged_cdf, "dr_noncrossing_diagnostics"
)
stopifnot(
  identical(as.numeric(rearranged_cdf), 3),
  identical(as.numeric(default_cdf), 2),
  identical(as.numeric(uncorrected_cdf), 2),
  rearrangement_diagnostics$dr_raw_crossing_pairs == 1L,
  rearrangement_diagnostics$dr_corrected_crossing_pairs == 0L,
  rearrangement_diagnostics$dr_max_adjustment == 0.4
)

# LPM CDF values are probability-bounded even when monotonicity correction is
# disabled. Diagnostics distinguish that bounding from crossing correction and
# measure the total change from the original fitted CDF.
lpm_cdf <- scalableCounterfactual:::inverse_step_cdf(
  1:4, c(-0.2, 0.7, 1.3, 0.6), 0.5, "none"
)
lpm_diagnostics <- attr(lpm_cdf, "dr_noncrossing_diagnostics")
stopifnot(
  identical(as.numeric(lpm_cdf), 2),
  lpm_diagnostics$dr_raw_crossing_pairs == 1L,
  lpm_diagnostics$dr_out_of_bounds_values == 2L,
  isTRUE(all.equal(lpm_diagnostics$dr_max_bound_adjustment, 0.3)),
  isTRUE(all.equal(lpm_diagnostics$dr_max_adjustment, 0.3))
)

set.seed(123)
n <- 400L
x <- rnorm(n)
latent <- 0.3 + 1.2 * x + rnorm(n, sd = 0.7)
X <- cbind("(Intercept)" = 1, x = x)
left_y <- pmax(0, latent)
nonconverged_cqr <- try(
  scalableCounterfactual:::fit_weighted_cqr(
    X, left_y, rep(1, n), rep(0, n), 0.5, "fn",
    dr_backend = "glm", dr_maxit = 1L
  ),
  silent = TRUE
)
stopifnot(
  inherits(nonconverged_cqr, "try-error"),
  grepl("selection model failed to converge", as.character(nonconverged_cqr))
)

invalid_left <- left_y
invalid_left[[1L]] <- -1
stopifnot(inherits(try(
  scalableCounterfactual:::fit_weighted_cqr(
    X, invalid_left, rep(1, n), rep(0, n), 0.5, "fn"
  ), silent = TRUE
), "try-error"))

right_limit <- 2
right_y <- pmin(right_limit, latent)
right_fit <- suppressWarnings(scalableCounterfactual:::fit_weighted_cqr(
  X, right_y, rep(1, n), rep(right_limit, n), c(0.3, 0.5, 0.7),
  "fn", right = TRUE
))
stopifnot(
  identical(rownames(right_fit$selection_sizes),
            paste0("tau_", signif(right_fit$taus, 6))),
  is.list(right_fit$selection_fit),
  identical(right_fit$selection_diagnostics$converged, TRUE),
  all(c(
    "converged", "boundary", "iterations", "backend", "fallback_used",
    "initial_backend", "initial_converged", "initial_boundary", "warnings"
  ) %in%
        names(right_fit$selection_diagnostics))
)

invalid_right <- right_y
invalid_right[[1L]] <- right_limit + 1
stopifnot(inherits(try(
  scalableCounterfactual:::fit_weighted_cqr(
    X, invalid_right, rep(1, n), rep(right_limit, n), 0.5, "fn",
    right = TRUE
  ), silent = TRUE
), "try-error"))

separated_x <- c(rep(-10, 100), rep(10, 100))
separated_X <- cbind("(Intercept)" = 1, x = separated_x)
separated_y <- c(rep(0, 100), rep(10, 100))
boundary_cqr <- try(
  scalableCounterfactual:::fit_weighted_cqr(
    separated_X, separated_y, rep(1, 200), rep(0, 200), 0.5, "fn",
    dr_backend = "glm", dr_maxit = 100L
  ),
  silent = TRUE
)
stopifnot(
  inherits(boundary_cqr, "try-error"),
  grepl(
    "boundary probabilities|failed to converge", as.character(boundary_cqr)
  )
)

# A large coefficient caused only by covariate units is not separation.
set.seed(991)
n <- 2000L
small_x <- rnorm(n, sd = 0.001)
uncensored <- rbinom(n, 1L, stats::plogis(0.1 + 800 * small_x))
scaled_X <- cbind("(Intercept)" = 1, x = small_x)
scaled_selection <- scalableCounterfactual:::fit_binary_threshold(
  scaled_X, uncensored, rep(1, n), "logit", backend = "glm"
)
stopifnot(
  abs(scaled_selection$coefficients[[2L]]) > 100,
  isTRUE(scaled_selection$converged),
  identical(scaled_selection$boundary, FALSE)
)
scaled_y <- ifelse(uncensored == 1L, 1 + 0.1 * rnorm(n), 0)
scaled_cqr <- suppressWarnings(scalableCounterfactual:::fit_weighted_cqr(
  scaled_X, scaled_y, rep(1, n), rep(0, n), 0.5, "fn",
  dr_backend = "glm"
))
stopifnot(
  all(is.finite(scaled_cqr$coefficients)),
  identical(scaled_cqr$selection_diagnostics$boundary, FALSE)
)

# A converged selection model is not separated merely because a large design
# range produces a few legitimate fitted probabilities below 1e-8.
set.seed(1771)
strong_x <- seq(-40, 40, length.out = 5000L)
strong_uncensored <- stats::rbinom(
  length(strong_x), 1L, stats::plogis(strong_x / 2)
)
strong_X <- cbind("(Intercept)" = 1, x = strong_x)
strong_selection <- scalableCounterfactual:::fit_binary_threshold(
  strong_X, strong_uncensored, rep(1, length(strong_x)), "logit",
  backend = "glm"
)
stopifnot(
  isTRUE(strong_selection$converged),
  identical(strong_selection$boundary, FALSE),
  min(stats::plogis(drop(strong_X %*% strong_selection$coefficients))) < 1e-8
)
strong_y <- ifelse(strong_uncensored == 1L, exp(rnorm(length(strong_x))), 0)
strong_cqr <- suppressWarnings(scalableCounterfactual:::fit_weighted_cqr(
  strong_X, strong_y, rep(1, length(strong_x)), rep(0, length(strong_x)),
  0.5, "fn", dr_backend = "glm"
))
stopifnot(all(is.finite(strong_cqr$coefficients)))

# Collapsing bootstrap duplicates must not make a full-rank CQR selection look
# underidentified merely because it contains exactly p distinct rows.  The
# effective selection has more than p observations once frequencies are used.
local({
  backend_name <- "reviewcqrfrequency"
  register_conditional_backend(
    backend_name, "distribution",
    function(X, response, weights, model, start, maxit, tolerance) {
      list(
        coefficients = c(0, 1),
        fitted.values = stats::plogis(drop(X %*% c(0, 1))),
        converged = TRUE,
        boundary = FALSE,
        iterations = 1L
      )
    },
    description = "deterministic CQR frequency regression-test backend"
  )
  on.exit(unregister_conditional_backend(
    backend_name, "distribution", quiet = TRUE
  ), add = TRUE)

  frequency_X <- cbind("(Intercept)" = 1, x = c(-2, -1, 1, 2))
  frequency_y <- c(0, 1, 1, 2)
  frequency <- c(1, 1, 2, 2)
  collapsed_fit <- suppressWarnings(
    scalableCounterfactual:::fit_weighted_cqr(
      frequency_X, frequency_y, frequency, rep(0, 4), 0.5, "fn",
      nsteps = 2L, dr_backend = backend_name,
      quantile_frequency = frequency
    )
  )
  expanded_index <- rep.int(seq_len(nrow(frequency_X)), frequency)
  expanded_fit <- suppressWarnings(
    scalableCounterfactual:::fit_weighted_cqr(
      frequency_X[expanded_index, , drop = FALSE],
      frequency_y[expanded_index], rep(1, length(expanded_index)),
      rep(0, length(expanded_index)), 0.5, "fn", nsteps = 2L,
      dr_backend = backend_name,
      quantile_frequency = rep(1, length(expanded_index))
    )
  )
  stopifnot(
    identical(collapsed_fit$selection_sizes[[1L]], 4L),
    identical(expanded_fit$selection_sizes[[1L]], 4L),
    isTRUE(all.equal(
      collapsed_fit$coefficients, expanded_fit$coefficients,
      tolerance = 1e-8, check.attributes = FALSE
    ))
  )
})

# Cox case weights use Breslow replication semantics. Frequency-compressed
# samples therefore retain the expanded sample size and event count, and
# floating-point noise at an exact CDF jump cannot change the selected time.
local({
  cox_X <- cbind(
    "(Intercept)" = 1,
    x = c(-1.4, -0.8, -0.2, 0.2, 0.6, 1, 1.3, 1.8)
  )
  cox_y <- c(1, 1, 2, 2, 3, 3, 4, 4)
  cox_event <- c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L)
  cox_weight <- c(0.7, 1.3, 2, 0.9, 1.5, 0.8, 1.2, 1.8)
  cox_frequency <- c(2L, 0L, 2L, 1L, 0L, 1L, 1L, 1L)
  active <- cox_frequency > 0L
  expanded_index <- rep.int(seq_along(cox_y), cox_frequency)

  collapsed_fit <- scalableCounterfactual:::fit_weighted_cox(
    cox_X[active, , drop = FALSE], cox_y[active], cox_weight[active],
    cox_event[active], frequency = cox_frequency[active]
  )
  expanded_fit <- scalableCounterfactual:::fit_weighted_cox(
    cox_X[expanded_index, , drop = FALSE], cox_y[expanded_index],
    cox_weight[expanded_index], cox_event[expanded_index]
  )
  stopifnot(
    identical(collapsed_fit$event_count, expanded_fit$event_count),
    identical(collapsed_fit$frequency_effective_n, 8),
    isTRUE(all.equal(
      collapsed_fit$coefficients, expanded_fit$coefficients,
      tolerance = 1e-12, check.attributes = FALSE
    )),
    isTRUE(all.equal(
      collapsed_fit$baseline_hazard, expanded_fit$baseline_hazard,
      tolerance = 1e-12, check.attributes = FALSE
    ))
  )

  jump_probability <- 0.69853019456604781
  collapsed_quantile <- scalableCounterfactual:::cox_marginal_quantiles(
    collapsed_fit, cox_X[active, , drop = FALSE],
    cox_weight[active] * cox_frequency[active], jump_probability
  )
  expanded_quantile <- scalableCounterfactual:::cox_marginal_quantiles(
    expanded_fit, cox_X[expanded_index, , drop = FALSE],
    cox_weight[expanded_index], jump_probability
  )
  stopifnot(
    identical(as.numeric(collapsed_quantile), 3),
    identical(as.numeric(collapsed_quantile), as.numeric(expanded_quantile))
  )

  km_y <- 1:2
  km_event <- c(1L, 1L)
  km_weight <- c(3.9960193827538748, 0.98867032871115956)
  km_frequency <- c(1L, 4L)
  km_probability <- 0.50259964935679968
  km_index <- rep.int(seq_along(km_y), km_frequency)
  collapsed_km <- scalableCounterfactual:::weighted_km_quantiles(
    km_y, km_event, km_weight * km_frequency, km_probability
  )
  expanded_km <- scalableCounterfactual:::weighted_km_quantiles(
    km_y[km_index], km_event[km_index], km_weight[km_index], km_probability
  )
  stopifnot(
    identical(as.numeric(collapsed_km), 1),
    identical(as.numeric(collapsed_km), as.numeric(expanded_km))
  )
})

# Quasi-complete separation can converge without a backend boundary flag:
# observations tied at x=0 contain both outcomes, while all remaining outcomes
# lie on the correct side of the fitted hyperplane.  The shared signed-margin
# diagnostic must still reject it.
quasi_x <- c(rep(-1, 100L), rep(0, 40L), rep(1, 100L))
quasi_response <- c(rep(0, 100L), rep(c(0, 1), each = 20L), rep(1, 100L))
quasi_X <- cbind("(Intercept)" = 1, x = quasi_x)
quasi_selection <- suppressWarnings(
  scalableCounterfactual:::fit_binary_threshold(
    quasi_X, quasi_response, rep(1, length(quasi_x)), "logit",
    backend = "glm"
  )
)
stopifnot(isTRUE(quasi_selection$boundary))

# Tail trimming belongs to QR only. CQR, location, and location-scale models
# use the complete midpoint grid, matching Counterfactual 1.2.
trimmed_control <- cf_control(
  nreg = 5L, trimming = 0.2, reported_quantiles = c(0.25, 0.5, 0.75),
  crossing_diagnostics = FALSE
)
untrimmed_control <- cf_control(
  nreg = 5L, trimming = 0, reported_quantiles = c(0.25, 0.5, 0.75),
  crossing_diagnostics = FALSE
)
stopifnot(
  identical(trimmed_control$conditional_quantiles, c(0.3, 0.5, 0.7)),
  identical(
    trimmed_control$full_conditional_quantiles,
    c(0.1, 0.3, 0.5, 0.7, 0.9)
  )
)
grid_data <- data.frame(
  y = c(seq_len(50L) / 10, 0.2 + seq_len(50L) / 10),
  x = rep(seq(-1, 1, length.out = 50L), 2L),
  group = rep(0:1, each = 50L),
  weight = 1
)
for (grid_model in c("loc", "locsca")) {
  trimmed_fit <- counterfactual_decompose(
    y ~ x, grid_data, "group", "weight", model = grid_model,
    control = trimmed_control
  )
  untrimmed_fit <- counterfactual_decompose(
    y ~ x, grid_data, "group", "weight", model = grid_model,
    control = untrimmed_control
  )
  stopifnot(
    identical(trimmed_fit$point$fits$group0$taus, c(0.1, 0.3, 0.5, 0.7, 0.9)),
    isTRUE(all.equal(
      trimmed_fit$point$effects, untrimmed_fit$point$effects,
      tolerance = 1e-12, check.attributes = FALSE
    ))
  )
}

# The default DR correction stays backward compatible with 1.0. Increasing
# rearrangement is the policy corresponding to Counterfactual 1.2's sorted
# marginal CDF implementation.
stopifnot(
  identical(cf_control()$dr_noncrossing, "cummax"),
  identical(
    as.numeric(scalableCounterfactual:::inverse_step_cdf(
      1:3, c(0.2, 0.8, 0.4), 0.5, "rearrange"
    )),
    3
  )
)
