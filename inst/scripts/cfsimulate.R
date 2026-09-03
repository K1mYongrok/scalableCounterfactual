#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scalableCounterfactual))

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_file)), "_cli_common.R"))
rm(script_file)

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript cfsimulate.R [--replications 100 --n-per-group 1000]",
    "    [--scenarios location_shift,composition_shift,combined_scale_shift]",
    "    [--solver fn --nreg 49 --quantiles 0.1,0.5,0.9]",
    "    [--weighted false --bootstrap-reps 0 --workers 1]",
    "    [--seed 20260809 --max-retries 1 --task-timeout-seconds Inf]",
    "    [--checkpoint-dir output/simulation/checkpoints]",
    "    --output output/simulation",
    sep = "\n"
  ))
}

args <- scalableCounterfactual:::parse_cli_args(commandArgs(trailingOnly = TRUE))
scalableCounterfactual:::validate_cli_args(args, c(
  "help", "replications", "n_per_group", "scenarios", "solver", "nreg",
  "trimming", "quantiles", "weighted", "bootstrap_reps",
  "bootstrap_scheme", "workers", "seed", "checkpoint_dir",
  "max_retries", "task_timeout_seconds", "output"
))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}
value <- scalableCounterfactual:::cli_value
scenarios <- strsplit(value(
  args, "scenarios",
  "location_shift,composition_shift,combined_scale_shift"
), ",", fixed = TRUE)[[1L]]
quantiles <- as.numeric(strsplit(value(
  args, "quantiles", "0.1,0.5,0.9"
), ",", fixed = TRUE)[[1L]])
output_dir <- value(args, "output", "output/simulation")
checkpoint_dir <- value(args, "checkpoint_dir", NULL)
simulation_output_paths <- stats::setNames(
  as.list(file.path(output_dir, c(
    "summary.csv", "raw.csv", "resources.csv", "curve_coverage.csv",
    "failures.csv", "simulation_validation.rds"
  ))), paste0("output_", seq_len(6L))
)
assert_distinct_cli_paths(c(
  list(checkpoint_directory = checkpoint_dir), simulation_output_paths
))
assert_cli_directories_not_below_files(
  list(checkpoint_directory = checkpoint_dir), simulation_output_paths
)
result <- simulate_counterfactual_validation(
  replications = value(args, "replications", 100L, as.integer),
  n_per_group = value(args, "n_per_group", 1000L, as.integer),
  scenarios = scenarios,
  solver = value(args, "solver", "fn"),
  nreg = value(args, "nreg", 49L, as.integer),
  trimming = value(args, "trimming", 0.02, as.numeric),
  reported_quantiles = quantiles,
  weighted = value(args, "weighted", FALSE, as_bool),
  bootstrap_reps = value(args, "bootstrap_reps", 0L, as.integer),
  bootstrap_scheme = value(
    args, "bootstrap_scheme", "counterfactual"
  ),
  workers = value(args, "workers", 1L, as.integer),
  seed = value(args, "seed", 20260809L, as.integer),
  progress = TRUE,
  checkpoint_dir = checkpoint_dir,
  max_retries = value(args, "max_retries", 1L, as.integer),
  task_timeout_seconds = value(
    args, "task_timeout_seconds", Inf, as.numeric
  )
)
write_simulation_validation(result, output_dir)
print(result)
message("Simulation outputs written to ", normalizePath(output_dir, winslash = "/"))
