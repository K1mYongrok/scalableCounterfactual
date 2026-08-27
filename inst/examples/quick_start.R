# Self-contained CPU example. It runs without project-specific data.

library(scalableCounterfactual)

set.seed(42)
n <- 600L
example_data <- data.frame(
  outcome = rnorm(n),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.4),
  group = rep(0:1, each = n / 2L),
  sampling_weight = runif(n, 0.5, 2)
)

fit <- counterfactual_decompose(
  outcome ~ x1 + x2,
  data = example_data,
  group = "group",
  weights = "sampling_weight",
  model = "qr",
  solver = "pfnb",
  control = cf_control(
    nreg = 9L,
    reported_quantiles = c(0.25, 0.5, 0.75),
    bootstrap_progress = FALSE
  ),
  bootstrap_reps = 3L,
  point_workers = 1L,
  bootstrap_workers = 1L,
  seed = 42L
)

print(summary(fit))

output_dir <- file.path(getwd(), "scalableCounterfactual_quick_start_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
grDevices::png(
  file.path(output_dir, "decomposition.png"),
  width = 1800, height = 1200, res = 180
)
plot(fit)
grDevices::dev.off()
write_cf_outputs(fit, output_dir)
message("Example outputs: ", normalizePath(output_dir, winslash = "/"))
