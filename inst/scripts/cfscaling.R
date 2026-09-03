#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(scalableCounterfactual)
})

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_file)), "_cli_common.R"))
rm(script_file)

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript cfscaling.R --data FILE --formula \"y ~ x1 + x2\"",
    "    --group BINARY_COLUMN [--weights WEIGHT_COLUMN|none]",
    "    [--solvers fn,qfnb,pfnb,proqreg,profn] [--reference pfnb]",
    "    [--sample-sizes 5000,20000,100000] [--repetitions 3 --warmup 1]",
    "    [--nreg 100 --trimming 0.005 --point-workers 1]",
    "    [--rss-poll-ms 25 --timeout-seconds Inf]",
    "    [--checkpoint FILE.rds --resume true] --output FILE.csv",
    "",
    "Every recorded run uses a fresh R process. Output includes sampled peak RSS.",
    sep = "\n"
  ))
}

args <- scalableCounterfactual:::parse_cli_args(commandArgs(trailingOnly = TRUE))
scalableCounterfactual:::validate_cli_args(args, c(
  "help", "data", "formula", "group", "weights", "solvers", "reference",
  "sample_sizes", "repetitions", "warmup", "nreg", "trimming", "seed",
  "point_workers", "rss_poll_ms", "timeout_seconds", "output",
  "checkpoint", "resume"
))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}
value <- scalableCounterfactual:::cli_value
data_path <- value(args, "data")
formula_text <- value(args, "formula")
group_column <- value(args, "group")
if (is.null(data_path) || is.null(formula_text) || is.null(group_column)) {
  usage()
  stop("--data, --formula, and --group are required", call. = FALSE)
}
weights_column <- value(args, "weights", NULL)
if (!is.null(weights_column) &&
    tolower(weights_column) %in% c("none", "null", "na")) {
  weights_column <- NULL
}
solvers <- strsplit(value(
  args, "solvers", "fn,qfnb,pfnb,proqreg,profn"
), ",", fixed = TRUE)[[1L]]
sample_sizes <- as.integer(strsplit(value(
  args, "sample_sizes", "5000,20000,100000"
), ",", fixed = TRUE)[[1L]])
reference_solver <- value(args, "reference", "pfnb")
repetitions <- value(args, "repetitions", 3L, as.integer)
warmup <- value(args, "warmup", 1L, as.integer)
nreg <- value(args, "nreg", 100L, as.integer)
trimming <- value(args, "trimming", 0.005, as.numeric)
seed <- value(args, "seed", 20260719L, as.integer)
point_workers <- value(args, "point_workers", 1L, as.integer)
rss_poll_ms <- value(args, "rss_poll_ms", 25L, as.integer)
timeout_seconds <- value(args, "timeout_seconds", Inf, as.numeric)
checkpoint_path <- value(args, "checkpoint", NULL)
resume <- value(args, "resume", TRUE, as_bool)
output_path <- value(args, "output", "output/benchmark/qr_scaling.csv")
summary_path <- sub("\\.csv$", "_summary.csv", output_path, ignore.case = TRUE)
if (identical(summary_path, output_path)) summary_path <- paste0(output_path, "_summary.csv")
assert_distinct_cli_paths(list(
  input_data = data_path,
  benchmark_output = output_path,
  summary_output = summary_path,
  checkpoint = checkpoint_path
))

if (!file.exists(data_path)) stop("data file not found: ", data_path, call. = FALSE)
message("Reading ", data_path)
data <- if (tolower(tools::file_ext(data_path)) == "rds") {
  readRDS(data_path)
} else {
  data.table::fread(data_path, showProgress = TRUE)
}
benchmark <- benchmark_qr_scaling(
  formula = stats::as.formula(formula_text),
  data = data,
  group = group_column,
  weights = weights_column,
  solvers = solvers,
  reference_solver = reference_solver,
  control = cf_control(
    nreg = nreg,
    trimming = trimming,
    crossing_diagnostics = FALSE
  ),
  sample_sizes = sample_sizes,
  seed = seed,
  point_workers = point_workers,
  repetitions = repetitions,
  warmup = warmup,
  rss_poll_interval_ms = rss_poll_ms,
  timeout_seconds = timeout_seconds,
  checkpoint_path = checkpoint_path,
  resume = resume
)
scalableCounterfactual:::atomic_write_output_files(
  dirname(output_path),
  c(basename(output_path), basename(summary_path)),
  required_files = c(basename(output_path), basename(summary_path)),
  writer = function(stage) {
    data.table::fwrite(benchmark$raw, file.path(stage, basename(output_path)))
    data.table::fwrite(
      benchmark$summary, file.path(stage, basename(summary_path))
    )
  }
)
print(benchmark$summary)
message("Scaling benchmark written to ", normalizePath(output_path, winslash = "/"))
message("Summary written to ", normalizePath(summary_path, winslash = "/"))
