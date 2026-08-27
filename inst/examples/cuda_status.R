# Optional CUDA check. CPU installation remains fully supported.

library(scalableCounterfactual)

status <- gpu_backend_status()
print(status)

if (!isTRUE(status$available)) {
  message("CUDA is unavailable; no GPU example was run.")
} else {
  set.seed(42)
  n <- 1000L
  example_data <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(0:1, each = n / 2L),
    weight = runif(n, 0.5, 2)
  )
  fit <- counterfactual_decompose(
    y ~ x,
    data = example_data,
    group = "group",
    weights = "weight",
    model = "qr",
    solver = "pfnb",
    control = cf_control(
      nreg = 19L,
      reported_quantiles = c(0.1, 0.5, 0.9),
      gpu_backend = "cuda",
      gpu_precision = "float64",
      bootstrap_progress = FALSE
    ),
    bootstrap_reps = 0L,
    seed = 42L
  )
  print(summary(fit))
}
