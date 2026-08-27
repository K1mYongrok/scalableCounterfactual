#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
value_after <- function(flag, default) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) return(default)
  args[[position + 1L]]
}
output_dir <- value_after(
  "--output", file.path("output", "release_validation")
)
workers <- as.integer(value_after(
  "--workers", max(1L, min(4L, parallel::detectCores(logical = FALSE)))
))
if (!is.finite(workers) || workers < 1L) stop("invalid --workers", call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(scalableCounterfactual)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
simulation_dir <- file.path(output_dir, "simulation")
checkpoint_dir <- file.path(simulation_dir, "checkpoints")
simulation <- simulate_counterfactual_validation(
  replications = if (quick) 10L else 100L,
  n_per_group = if (quick) 300L else 1000L,
  scenarios = c(
    "location_shift", "composition_shift", "combined_scale_shift"
  ),
  solver = "fn",
  nreg = if (quick) 19L else 49L,
  trimming = 0.02,
  reported_quantiles = c(0.1, 0.5, 0.9),
  weighted = TRUE,
  bootstrap_reps = 0L,
  workers = workers,
  seed = 20260810L,
  progress = TRUE,
  checkpoint_dir = checkpoint_dir,
  max_retries = 1L
)
write_simulation_validation(simulation, simulation_dir)

set.seed(20260810L)
n <- if (quick) 1000L else 5000L
i <- seq_len(n)
benchmark_data <- data.frame(
  y = 0.2 * (i %% 2L) + 0.4 * sin(i / 17) + rnorm(n),
  x1 = rnorm(n),
  x2 = runif(n, -1, 1),
  x3 = as.numeric(i %% 3L == 0L),
  group = rep(0:1, length.out = n),
  weight = 0.5 + (i %% 13L) / 10
)
benchmark <- benchmark_qr_solvers_repeated(
  y ~ x1 + x2 + x3,
  data = benchmark_data,
  group = "group",
  weights = "weight",
  solvers = c("br", "fn", "qfnb", "pfnb", "onestep"),
  reference_solver = "br",
  control = cf_control(
    nreg = if (quick) 19L else 100L,
    reported_quantiles = seq(0.1, 0.9, by = 0.1),
    crossing_diagnostics = FALSE,
    bootstrap_progress = FALSE
  ),
  repetitions = if (quick) 1L else 3L,
  warmup = if (quick) 0L else 1L,
  randomize_order = TRUE,
  seed = 20260810L
)
fwrite(benchmark$raw, file.path(output_dir, "solver_benchmark_raw.csv"))
fwrite(benchmark$summary, file.path(output_dir, "solver_benchmark_summary.csv"))

metadata <- data.frame(
  item = c(
    "timestamp", "package_version", "output_schema_version", "quick",
    "workers", "R_version", "platform"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    as.character(packageVersion("scalableCounterfactual")),
    "1.0", as.character(quick), as.character(workers),
    R.version.string, R.version$platform
  ),
  stringsAsFactors = FALSE
)
fwrite(metadata, file.path(output_dir, "release_metadata.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("Release validation written to ", normalizePath(
  output_dir, winslash = "/", mustWork = TRUE
))
