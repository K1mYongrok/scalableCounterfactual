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
    "  Rscript cfbenchmark.R --data FILE --formula \"y ~ x1 + x2\"",
    "    --group BINARY_COLUMN [--weights WEIGHT_COLUMN|none]",
    "    --solvers br,fn,pfn,qfnb,pfnb,proqreg,profn,onestep,auto",
    "    --reference br [--sample-n 5000] [--repetitions 5 --warmup 1]",
    "    [--nreg 100 --trimming 0.005 --point-workers 1]",
    "    --output FILE.csv",
    "",
    "Input FILE may be CSV or RDS. Group must be coded 0 and 1.",
    sep = "\n"
  ))
}

read_analysis_data <- function(path) {
  if (!file.exists(path)) stop("data file not found: ", path, call. = FALSE)
  if (tolower(tools::file_ext(path)) == "rds") {
    value <- readRDS(path)
    if (!is.data.frame(value)) stop("RDS input must contain a data frame")
    return(value)
  }
  data.table::fread(path, showProgress = TRUE)
}

args <- scalableCounterfactual:::parse_cli_args(commandArgs(trailingOnly = TRUE))
scalableCounterfactual:::validate_cli_args(args, c(
  "help", "data", "formula", "group", "weights", "solvers", "reference",
  "sample_n", "nreg", "trimming", "seed", "point_workers", "repetitions",
  "warmup", "qr_precondition", "output"
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
  args, "solvers", "br,fn,pfn,qfnb,pfnb,proqreg,profn,onestep,auto"
), ",", fixed = TRUE)[[1L]]
reference_solver <- value(args, "reference", "br")
sample_n <- value(args, "sample_n", NULL, as.integer)
nreg <- value(args, "nreg", 100L, as.integer)
trimming <- value(args, "trimming", 0.005, as.numeric)
seed <- value(args, "seed", 20260719L, as.integer)
point_workers <- value(args, "point_workers", 1L, as.integer)
repetitions <- value(args, "repetitions", 5L, as.integer)
warmup <- value(args, "warmup", 1L, as.integer)
qr_precondition <- value(args, "qr_precondition", TRUE, as_bool)
output_path <- value(args, "output", "output/qr_solvers.csv")
summary_path <- sub("\\.csv$", "_summary.csv", output_path, ignore.case = TRUE)
if (identical(summary_path, output_path)) summary_path <- paste0(output_path, "_summary.csv")
assert_distinct_cli_paths(list(
  input_data = data_path,
  benchmark_output = output_path,
  summary_output = summary_path
))

message("Reading ", data_path)
data <- read_analysis_data(data_path)
benchmark <- benchmark_qr_solvers_repeated(
  formula = stats::as.formula(formula_text),
  data = data,
  group = group_column,
  weights = weights_column,
  solvers = solvers,
  reference_solver = reference_solver,
  control = cf_control(
    nreg = nreg,
    trimming = trimming,
    qr_precondition = qr_precondition
  ),
  sample_n = sample_n,
  seed = seed,
  point_workers = point_workers,
  repetitions = repetitions,
  warmup = warmup,
  randomize_order = TRUE
)

scalableCounterfactual:::atomic_write_output_files(
  dirname(output_path),
  c(basename(output_path), basename(summary_path)),
  required_files = c(basename(output_path), basename(summary_path)),
  writer = function(stage) {
    fwrite(benchmark$raw, file.path(stage, basename(output_path)))
    fwrite(benchmark$summary, file.path(stage, basename(summary_path)))
  }
)
print(benchmark$summary)
message("Benchmark written to ", normalizePath(output_path, winslash = "/"))
message("Summary written to ", normalizePath(summary_path, winslash = "/"))
