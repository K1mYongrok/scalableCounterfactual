#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) return(default)
  args[[position + 1L]]
}
output_dir <- value_after(
  "--output", file.path("output", "release_validation", "coverage")
)
workers <- as.integer(value_after(
  "--workers", max(1L, min(4L, parallel::detectCores(logical = FALSE)))
))
replications <- as.integer(value_after("--replications", 30L))
bootstrap_reps <- as.integer(value_after("--bootstrap-reps", 49L))
bootstrap_scheme <- value_after("--bootstrap-scheme", "counterfactual")
if (any(!is.finite(c(workers, replications, bootstrap_reps))) ||
    workers < 1L || replications < 2L || bootstrap_reps < 2L) {
  stop("invalid release coverage settings", call. = FALSE)
}
if (!bootstrap_scheme %in% c("counterfactual", "empirical", "multiplier")) {
  stop("invalid --bootstrap-scheme", call. = FALSE)
}

suppressPackageStartupMessages(library(scalableCounterfactual))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result <- simulate_counterfactual_validation(
  replications = replications,
  n_per_group = 500L,
  scenarios = c(
    "location_shift", "composition_shift", "combined_scale_shift"
  ),
  solver = "fn",
  nreg = 29L,
  trimming = 0.02,
  reported_quantiles = c(0.1, 0.5, 0.9),
  weighted = TRUE,
  bootstrap_reps = bootstrap_reps,
  bootstrap_scheme = bootstrap_scheme,
  workers = workers,
  seed = 20260811L,
  progress = TRUE,
  checkpoint_dir = file.path(output_dir, "checkpoints"),
  max_retries = 2L
)
write_simulation_validation(result, output_dir)
data.table::fwrite(data.frame(
  item = c(
    "timestamp", "package_version", "replications", "n_per_group",
    "bootstrap_reps", "bootstrap_scheme", "workers", "R_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    as.character(packageVersion("scalableCounterfactual")),
    replications, 500L, bootstrap_reps, bootstrap_scheme, workers,
    R.version.string
  )
), file.path(output_dir, "coverage_metadata.csv"))
message("Release coverage validation written to ", normalizePath(
  output_dir, winslash = "/", mustWork = TRUE
))
