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
    "  Rscript cfdecomp.R --data FILE --formula \"y ~ x1 + x2\"",
    "    --group BINARY_COLUMN [--weights WEIGHT_COLUMN|none]",
    "    --model qr|cqr|loc|locsca|cox|logit|probit|cloglog|lpm [--solver SOLVER]",
    "    [--censoring COLUMN|SCALAR --event COLUMN]",
    "    [--cox-boundary na|error|cap]",
    "    [--reps N --point-workers N --bootstrap-workers N]",
    "    [--bootstrap-max-retries N --bootstrap-progress true|false]",
    "    [--nreg N --trimming P --alpha P]",
    "    [--robust-se true|false]",
    "    [--reported-quantiles 0.1,0.5,0.9]",
    "    [--bootstrap-scheme counterfactual|empirical|multiplier]",
    "    [--qr-precondition true|false]",
    "    [--onestep-first-solver auto|br|fn|pfn]",
    "    [--onestep-bandwidth hall_sheather|bofinger]",
    "    [--qr-bootstrap-engine standard|xy_preprocess|onestep|auto]",
    "    [--marginal-method auto|matrix|chunked]",
    "    [--marginal-chunk-rows N --marginal-matrix-max-mb MB]",
    "    [--crossing-diagnostics true|false --crossing-tolerance X]",
    "    [--quantile-noncrossing none|rearrange]",
    "    [--dr-noncrossing cummax|rearrange|isotonic|none]",
    "    [--cqr-right true|false --cqr-nsteps N]",
    "    [--cqr-first-cut P --cqr-later-cut P]",
    "    [--linear-backend auto|qr|chol|fastglm]",
    "    [--dr-backend auto|glm|fastglm|speedglm|cuda]",
    "    [--dr-workers N --dr-warm-start true|false]",
    "    [--dr-maxit N --dr-tolerance X --dr-precondition true|false]",
    "    [--gpu-backend cpu|cuda --gpu-precision float64|float32]",
    "    [--gpu-python FILE --gpu-python-path DIR --gpu-module-path FILE]",
    "    [--gpu-block-columns N --gpu-qr-rho X --gpu-qr-maxit N]",
    "    [--gpu-qr-tolerance X --gpu-qr-allow-nonconvergence true|false]",
    "    --output DIR [--checkpoint-dir DIR]",
    "",
    "Input FILE may be CSV or RDS. Group must be coded 0 and 1.",
    "QR solvers: br, fn, pfn, qfnb, pfnb, proqreg, profn, onestep, cuda_admm, auto",
    sep = "\n"
  ))
}

as_number_list <- function(x) {
  as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
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
  "help", "data", "formula", "group", "weights", "model", "solver",
  "censoring", "event", "cox_boundary",
  "reps", "point_workers", "bootstrap_workers", "nreg", "trimming",
  "bootstrap_max_retries", "bootstrap_progress",
  "alpha", "robust_se", "reported_quantiles", "seed", "bootstrap_scheme",
  "weighted_bootstrap", "legacy_qr_shift", "legacy_weighted_quantile",
  "qr_precondition",
  "onestep_first_solver", "onestep_bandwidth", "qr_bootstrap_engine",
  "marginal_method", "marginal_chunk_rows", "marginal_matrix_max_mb",
  "marginal_histogram_bins", "marginal_candidate_max",
  "crossing_diagnostics", "crossing_tolerance",
  "quantile_noncrossing", "dr_noncrossing", "cqr_right", "cqr_nsteps",
  "cqr_first_cut", "cqr_later_cut",
  "linear_backend", "dr_backend", "dr_workers", "dr_warm_start",
  "dr_maxit", "dr_tolerance",
  "dr_precondition",
  "gpu_backend", "gpu_precision", "gpu_python", "gpu_python_path",
  "gpu_module_path", "gpu_block_columns", "gpu_qr_rho",
  "gpu_qr_maxit", "gpu_qr_tolerance", "gpu_qr_allow_nonconvergence",
  "output", "checkpoint_dir"
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
model <- value(args, "model", "qr")
solver <- value(args, "solver", if (model %in% c("qr", "cqr")) "auto" else NULL)
censoring <- value(args, "censoring", NULL)
event <- value(args, "event", NULL)
cox_boundary <- value(args, "cox_boundary", "na")
reps <- value(args, "reps", 0L, as.integer)
point_workers <- value(args, "point_workers", 1L, as.integer)
bootstrap_workers <- value(args, "bootstrap_workers", 1L, as.integer)
bootstrap_max_retries <- value(args, "bootstrap_max_retries", 2L, as.integer)
bootstrap_progress <- value(args, "bootstrap_progress", TRUE, as_bool)
nreg <- value(args, "nreg", 100L, as.integer)
trimming <- value(args, "trimming", 0.005, as.numeric)
alpha <- value(args, "alpha", 0.05, as.numeric)
robust_se <- value(args, "robust_se", FALSE, as_bool)
reported_quantiles <- value(
  args, "reported_quantiles", "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
  as_number_list
)
seed <- value(args, "seed", 20260719L, as.integer)
bootstrap_scheme <- value(args, "bootstrap_scheme", NULL)
weighted_bootstrap <- value(args, "weighted_bootstrap", NULL, as_bool)
legacy_shift <- value(args, "legacy_qr_shift", TRUE, as_bool)
legacy_weighted_quantile <- value(
  args, "legacy_weighted_quantile", FALSE, as_bool
)
qr_precondition <- value(args, "qr_precondition", TRUE, as_bool)
onestep_first_solver <- value(args, "onestep_first_solver", "auto")
onestep_bandwidth <- value(args, "onestep_bandwidth", "hall_sheather")
qr_bootstrap_engine <- value(args, "qr_bootstrap_engine", "standard")
marginal_method <- value(args, "marginal_method", "auto")
marginal_chunk_rows <- value(args, "marginal_chunk_rows", 50000L, as.integer)
marginal_matrix_max_mb <- value(
  args, "marginal_matrix_max_mb", 512, as.numeric
)
marginal_histogram_bins <- value(
  args, "marginal_histogram_bins", 262144L, as.integer
)
marginal_candidate_max <- value(
  args, "marginal_candidate_max", 2000000L, as.integer
)
crossing_diagnostics <- value(args, "crossing_diagnostics", TRUE, as_bool)
crossing_tolerance <- value(args, "crossing_tolerance", 1e-8, as.numeric)
quantile_noncrossing <- value(args, "quantile_noncrossing", "none")
dr_noncrossing <- value(args, "dr_noncrossing", "cummax")
cqr_right <- value(args, "cqr_right", FALSE, as_bool)
cqr_nsteps <- value(args, "cqr_nsteps", 3L, as.integer)
cqr_first_cut <- value(args, "cqr_first_cut", 0.1, as.numeric)
cqr_later_cut <- value(args, "cqr_later_cut", 0.05, as.numeric)
linear_backend <- value(args, "linear_backend", "auto")
dr_backend <- value(args, "dr_backend", "auto")
dr_workers <- value(args, "dr_workers", 1L, as.integer)
dr_warm_start <- value(args, "dr_warm_start", TRUE, as_bool)
dr_maxit <- value(args, "dr_maxit", 100L, as.integer)
dr_tolerance <- value(args, "dr_tolerance", 1e-8, as.numeric)
dr_precondition <- value(args, "dr_precondition", TRUE, as_bool)
gpu_backend <- value(args, "gpu_backend", "cpu")
gpu_precision <- value(args, "gpu_precision", "float64")
gpu_python <- value(args, "gpu_python", NULL)
gpu_python_path <- value(args, "gpu_python_path", NULL)
gpu_module_path <- value(args, "gpu_module_path", NULL)
gpu_block_columns <- value(args, "gpu_block_columns", 16L, as.integer)
gpu_qr_rho <- value(args, "gpu_qr_rho", 1, as.numeric)
gpu_qr_maxit <- value(args, "gpu_qr_maxit", 5000L, as.integer)
gpu_qr_tolerance <- value(args, "gpu_qr_tolerance", 1e-6, as.numeric)
gpu_qr_allow_nonconvergence <- value(
  args, "gpu_qr_allow_nonconvergence", FALSE, as_bool
)
output_dir <- value(args, "output", "output/counterfactual")
checkpoint_dir <- value(
  args, "checkpoint_dir", file.path(output_dir, "checkpoints")
)
managed_output_paths <- stats::setNames(
  file.path(output_dir, c(
    "decomposition.csv", "distribution_diagnostics.csv",
    "point_resources.csv", "marginalization_diagnostics.csv",
    "run_metadata.csv", "fit.rds", "bootstrap_resources.csv",
    "bootstrap_failures.csv", "functional_effect_tests.csv",
    "quantile_crossing_diagnostics.csv"
  )),
  paste0("output_", seq_len(10L))
)
assert_distinct_cli_paths(c(
  list(input_data = data_path),
  as.list(managed_output_paths)
))
assert_distinct_cli_paths(c(
  list(input_data = data_path, checkpoint_directory = checkpoint_dir),
  as.list(managed_output_paths)
))
assert_cli_directories_not_below_files(
  list(output_directory = output_dir, checkpoint_directory = checkpoint_dir),
  c(list(input_data = data_path), as.list(managed_output_paths))
)

message("Reading ", data_path)
data <- read_analysis_data(data_path)
if (!is.null(censoring) && !censoring %in% names(data)) {
  numeric_censoring <- suppressWarnings(as.numeric(censoring))
  if (!is.finite(numeric_censoring)) {
    stop("--censoring must be a column name or finite scalar", call. = FALSE)
  }
  censoring <- numeric_censoring
}
formula <- stats::as.formula(formula_text)
control_arguments <- list(
  nreg = nreg,
  trimming = trimming,
  reported_quantiles = reported_quantiles,
  alpha = alpha,
  robust_se = robust_se,
  weighted_bootstrap = weighted_bootstrap,
  legacy_qr_shift = legacy_shift,
  legacy_weighted_quantile = legacy_weighted_quantile,
  qr_precondition = qr_precondition,
  onestep_first_solver = onestep_first_solver,
  onestep_bandwidth = onestep_bandwidth,
  qr_bootstrap_engine = qr_bootstrap_engine,
  bootstrap_max_retries = bootstrap_max_retries,
  bootstrap_progress = bootstrap_progress,
  marginal_method = marginal_method,
  marginal_chunk_rows = marginal_chunk_rows,
  marginal_matrix_max_mb = marginal_matrix_max_mb,
  marginal_histogram_bins = marginal_histogram_bins,
  marginal_candidate_max = marginal_candidate_max,
  crossing_diagnostics = crossing_diagnostics,
  crossing_tolerance = crossing_tolerance,
  quantile_noncrossing = quantile_noncrossing,
  dr_noncrossing = dr_noncrossing,
  cqr_right = cqr_right,
  cqr_nsteps = cqr_nsteps,
  cqr_first_cut = cqr_first_cut,
  cqr_later_cut = cqr_later_cut,
  linear_backend = linear_backend,
  dr_backend = dr_backend,
  dr_workers = dr_workers,
  dr_warm_start = dr_warm_start,
  dr_maxit = dr_maxit,
  dr_tolerance = dr_tolerance,
  dr_precondition = dr_precondition,
  gpu_backend = gpu_backend,
  gpu_precision = gpu_precision,
  gpu_python = gpu_python,
  gpu_python_path = gpu_python_path,
  gpu_module_path = gpu_module_path,
  gpu_block_columns = gpu_block_columns,
  gpu_qr_rho = gpu_qr_rho,
  gpu_qr_maxit = gpu_qr_maxit,
  gpu_qr_tolerance = gpu_qr_tolerance,
  gpu_qr_allow_nonconvergence = gpu_qr_allow_nonconvergence,
  cox_boundary = cox_boundary
)
if (!is.null(bootstrap_scheme)) {
  control_arguments$bootstrap_scheme <- bootstrap_scheme
}
control <- do.call(cf_control, control_arguments)

fit <- counterfactual_decompose(
  formula = formula,
  data = data,
  group = group_column,
  weights = weights_column,
  model = model,
  solver = solver,
  censoring = censoring,
  event = event,
  control = control,
  bootstrap_reps = reps,
  point_workers = point_workers,
  bootstrap_workers = bootstrap_workers,
  checkpoint_dir = checkpoint_dir,
  seed = seed
)
write_cf_outputs(fit, output_dir)
print(fit)
message("Outputs written to ", normalizePath(output_dir, winslash = "/"))
