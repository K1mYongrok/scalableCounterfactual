library(scalableCounterfactual)

set.seed(7)
n <- 500L
d <- data.frame(
  y = rnorm(n),
  x = rnorm(n),
  group = rbinom(n, 1, 0.5),
  weight = runif(n, 0.5, 1.5)
)
result <- benchmark_qr_solvers(
  y ~ x,
  d,
  group = "group",
  weights = "weight",
  solvers = c(
    "br", "fn", "pfn", "qfnb", "pfnb",
    "proqreg", "profn", "onestep", "auto"
  ),
  reference_solver = "br",
  control = cf_control(nreg = 21L, trimming = 0.02),
  point_workers = 1L
)
stopifnot(all(result$status == "ok"))
stopifnot(result$max_abs_effect_difference[[1L]] == 0)
stopifnot(is.na(result$nonzero_convergence_flags[result$solver == "br"]))
stopifnot(isTRUE(result$convergence_diagnostics_available[result$solver == "pfnb"]))
stopifnot(result$group0_exact[result$solver == "profn"])
stopifnot(!result$group0_exact[result$solver == "onestep"])
stopifnot(result$resolved_group0_solver[result$solver == "auto"] == "qfnb")
stopifnot(all(result$package_version ==
  as.character(utils::packageVersion("scalableCounterfactual"))))

repeated <- benchmark_qr_solvers_repeated(
  y ~ x,
  d,
  group = "group",
  weights = "weight",
  solvers = c("fn", "pfnb"),
  reference_solver = "pfnb",
  control = cf_control(nreg = 9L, trimming = 0.05),
  point_workers = 1L,
  repetitions = 2L,
  warmup = 0L,
  randomize_order = TRUE
)
stopifnot(nrow(repeated$raw) == 4L)
stopifnot(nrow(repeated$summary) == 2L)

linear_backends <- benchmark_conditional_backends(
  y ~ x,
  d,
  group = "group",
  weights = "weight",
  model = "loc",
  backends = c("qr", "chol"),
  reference_backend = "qr",
  control = cf_control(),
  sample_n = 300L,
  point_workers = 1L
)
stopifnot(all(linear_backends$status == "ok"))
stopifnot(max(linear_backends$max_abs_effect_difference) < 1e-10)

if (requireNamespace("fastglm", quietly = TRUE)) {
  dr_backends <- suppressWarnings(benchmark_conditional_backends(
    y ~ x,
    d,
    group = "group",
    weights = "weight",
    model = "logit",
    backends = c("glm", "fastglm"),
    reference_backend = "glm",
    control = cf_control(nreg = 9L, dr_warm_start = FALSE),
    sample_n = 300L,
    point_workers = 1L
  ))
  stopifnot(all(dr_backends$status == "ok"))
  stopifnot(max(dr_backends$max_abs_effect_difference) < 1e-5)
}

if (requireNamespace("processx", quietly = TRUE)) {
  scaling_checkpoint <- tempfile("cf_scaling_checkpoint_", fileext = ".rds")
  isolated <- benchmark_qr_scaling(
    y ~ x,
    d,
    group = "group",
    weights = "weight",
    solvers = c("fn", "qfnb"),
    reference_solver = "fn",
    control = cf_control(
      nreg = 9L,
      trimming = 0.05,
      reported_quantiles = c(0.25, 0.5, 0.75),
      crossing_diagnostics = FALSE
    ),
    sample_sizes = 300L,
    point_workers = 1L,
    repetitions = 1L,
    warmup = 0L,
    rss_poll_interval_ms = 10L,
    checkpoint_path = scaling_checkpoint
  )
  stopifnot(all(isolated$raw$status == "ok"))
  stopifnot(all(is.finite(isolated$raw$peak_process_rss_mb)))
  stopifnot(all(isolated$raw$peak_process_rss_mb > 0))
  stopifnot(nrow(isolated$summary) == 2L)
  stopifnot(isolated$raw$max_abs_effect_difference[
    isolated$raw$solver == "fn"
  ] == 0)
  isolated_resumed <- benchmark_qr_scaling(
    y ~ x,
    d,
    group = "group",
    weights = "weight",
    solvers = c("fn", "qfnb"),
    reference_solver = "fn",
    control = cf_control(
      nreg = 9L,
      trimming = 0.05,
      reported_quantiles = c(0.25, 0.5, 0.75),
      crossing_diagnostics = FALSE
    ),
    sample_sizes = 300L,
    point_workers = 1L,
    repetitions = 1L,
    warmup = 0L,
    rss_poll_interval_ms = 10L,
    checkpoint_path = scaling_checkpoint
  )
  stopifnot(identical(isolated$raw, isolated_resumed$raw))
  unlink(scaling_checkpoint)
}
