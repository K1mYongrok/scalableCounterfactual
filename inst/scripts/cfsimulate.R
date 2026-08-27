#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scalableCounterfactual))

as_bool <- function(x) {
  normalized <- tolower(as.character(x))
  if (!normalized %in% c("1", "0", "true", "false", "yes", "no", "y", "n")) {
    stop("expected a true/false value, got: ", x, call. = FALSE)
  }
  normalized %in% c("1", "true", "yes", "y")
}

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
  checkpoint_dir = value(args, "checkpoint_dir", NULL),
  max_retries = value(args, "max_retries", 1L, as.integer),
  task_timeout_seconds = value(
    args, "task_timeout_seconds", Inf, as.numeric
  )
)
output_dir <- value(args, "output", "output/simulation")
write_simulation_validation(result, output_dir)
print(result)
message("Simulation outputs written to ", normalizePath(output_dir, winslash = "/"))
