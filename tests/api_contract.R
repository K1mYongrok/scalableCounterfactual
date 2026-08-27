library(scalableCounterfactual)

required_exports <- c(
  "benchmark_conditional_backends", "benchmark_qr_scaling",
  "benchmark_qr_solvers", "benchmark_qr_solvers_repeated", "cf_control",
  "conditional_backend_registry", "counterfactual_decompose",
  "fit_weighted_qr", "functional_effect_tests", "gpu_backend_status",
  "qr_solver_registry", "register_conditional_backend",
  "register_qr_solver", "simulate_counterfactual_validation",
  "unregister_conditional_backend",
  "unregister_qr_solver", "write_cf_outputs",
  "write_simulation_validation"
)
actual_exports <- sort(getNamespaceExports("scalableCounterfactual"))
stopifnot(identical(sort(required_exports), actual_exports))
stopifnot(
  is.function(getS3method("summary", "cfdecomp")),
  is.function(getS3method("plot", "cfdecomp")),
  is.function(getS3method("print", "cfdecomp")),
  is.function(getS3method("as.data.frame", "cfdecomp"))
)

s3_contract <- list(
  "print.cfdecomp" = c("x", "..."),
  "as.data.frame.cfdecomp" = c("x", "row.names", "optional", "..."),
  "summary.cfdecomp" = c("object", "effects", "quantiles", "..."),
  "plot.cfdecomp" = c(
    "x", "effects", "interval", "col", "lwd", "pch", "xlab", "ylab",
    "main", "..."
  ),
  "print.summary.cfdecomp" = c("x", "...")
)
s3_generics <- c("print", "as.data.frame", "summary", "plot", "print")
s3_classes <- c(
  "cfdecomp", "cfdecomp", "cfdecomp", "cfdecomp", "summary.cfdecomp"
)
for (method_index in seq_along(s3_contract)) {
  actual_arguments <- names(formals(getS3method(
    s3_generics[[method_index]], s3_classes[[method_index]]
  )))
  stopifnot(identical(actual_arguments, s3_contract[[method_index]]))
}

api_contract <- list(
  benchmark_conditional_backends = c(
    "formula", "data", "group", "weights", "model", "backends",
    "reference_backend", "control", "sample_n", "seed", "point_workers"
  ),
  benchmark_qr_scaling = c(
    "formula", "data", "group", "weights", "solvers", "reference_solver",
    "control", "sample_sizes", "seed", "point_workers", "repetitions",
    "warmup", "randomize_order", "rss_poll_interval_ms", "timeout_seconds",
    "checkpoint_path", "resume"
  ),
  benchmark_qr_solvers = c(
    "formula", "data", "group", "weights", "solvers", "reference_solver",
    "control", "sample_n", "seed", "point_workers"
  ),
  benchmark_qr_solvers_repeated = c(
    "formula", "data", "group", "weights", "solvers", "reference_solver",
    "control", "sample_n", "seed", "point_workers", "repetitions",
    "warmup", "randomize_order"
  ),
  cf_control = c(
    "nreg", "trimming", "reported_quantiles", "alpha",
    "weighted_bootstrap", "bootstrap_scheme", "legacy_qr_shift",
    "legacy_weighted_quantile", "qr_precondition", "onestep_first_solver",
    "onestep_bandwidth", "qr_bootstrap_engine", "robust_se",
    "bootstrap_max_retries", "bootstrap_progress", "marginal_method",
    "marginal_chunk_rows", "marginal_matrix_max_mb",
    "marginal_histogram_bins", "marginal_candidate_max",
    "crossing_diagnostics", "crossing_tolerance", "quantile_noncrossing",
    "cqr_right", "cqr_nsteps", "cqr_first_cut", "cqr_later_cut",
    "cox_boundary", "linear_backend", "dr_backend", "dr_workers",
    "dr_warm_start", "dr_maxit", "dr_tolerance", "dr_precondition",
    "gpu_backend", "gpu_precision", "gpu_python", "gpu_python_path",
    "gpu_module_path", "gpu_block_columns", "gpu_qr_rho",
    "gpu_qr_maxit", "gpu_qr_tolerance", "gpu_qr_allow_nonconvergence"
  ),
  conditional_backend_registry = character(),
  counterfactual_decompose = c(
    "formula", "data", "group", "weights", "model", "solver", "control",
    "bootstrap_reps", "point_workers", "bootstrap_workers", "checkpoint_dir",
    "seed", "censoring", "event"
  ),
  fit_weighted_qr = c(
    "X", "y", "weights", "taus", "solver", "precondition",
    "onestep_first_solver", "onestep_bandwidth", "gpu_control"
  ),
  functional_effect_tests = c("object", "constants", "quantile_range"),
  gpu_backend_status = c("python", "python_path", "module_path"),
  qr_solver_registry = character(),
  register_conditional_backend = c(
    "name", "type", "fit", "objective_preserving", "description",
    "overwrite", "dependencies", "version"
  ),
  register_qr_solver = c(
    "name", "fit", "exact", "process_aware", "description", "overwrite",
    "dependencies", "version"
  ),
  simulate_counterfactual_validation = c(
    "replications", "n_per_group", "scenarios", "solver", "nreg",
    "trimming", "reported_quantiles", "weighted", "bootstrap_reps",
    "bootstrap_scheme", "workers", "seed", "progress", "checkpoint_dir",
    "max_retries", "task_timeout_seconds"
  ),
  unregister_conditional_backend = c("name", "type", "quiet"),
  unregister_qr_solver = c("name", "quiet"),
  write_cf_outputs = c("object", "output_dir"),
  write_simulation_validation = c("object", "output_dir")
)
for (function_name in names(api_contract)) {
  actual_arguments <- names(formals(getExportedValue(
    "scalableCounterfactual", function_name
  )))
  if (is.null(actual_arguments)) actual_arguments <- character()
  stopifnot(identical(actual_arguments, api_contract[[function_name]]))
}

split_contract_arguments <- function(value) {
  if (is.na(value) || !nzchar(value)) return(character())
  strsplit(value, ";", fixed = TRUE)[[1L]]
}
schema_dir <- system.file("schema", package = "scalableCounterfactual")
stopifnot(nzchar(schema_dir))
public_schema <- data.table::fread(
  file.path(schema_dir, "public_api_1.0.csv"), na.strings = ""
)
stopifnot(identical(public_schema[["function"]], names(api_contract)))
for (row in seq_len(nrow(public_schema))) {
  stopifnot(identical(
    split_contract_arguments(public_schema$ordered_arguments[[row]]),
    api_contract[[public_schema[["function"]][[row]]]]
  ))
}

s3_schema <- data.table::fread(file.path(schema_dir, "s3_methods_1.0.csv"))
schema_method_names <- paste(s3_schema$generic, s3_schema$class, sep = ".")
stopifnot(identical(schema_method_names, names(s3_contract)))
for (row in seq_len(nrow(s3_schema))) {
  stopifnot(identical(
    split_contract_arguments(s3_schema$ordered_arguments[[row]]),
    s3_contract[[schema_method_names[[row]]]]
  ))
}

set.seed(150)
n <- 160L
contract_data <- data.frame(
  y = rnorm(n) + rep(c(0, 0.2), each = n / 2L),
  x = rnorm(n),
  group = rep(0:1, each = n / 2L),
  weight = runif(n, 0.5, 2)
)
fit <- counterfactual_decompose(
  y ~ x,
  data = contract_data,
  group = "group",
  weights = "weight",
  model = "qr",
  solver = "fn",
  control = cf_control(
    nreg = 9L,
    reported_quantiles = c(0.25, 0.5, 0.75),
    legacy_qr_shift = FALSE,
    crossing_diagnostics = FALSE,
    bootstrap_progress = FALSE
  )
)
stopifnot(
  inherits(fit, "cfdecomp"),
  identical(fit$metadata$output_schema_version, "1.0")
)

output_dir <- tempfile("cf_api_contract_")
write_cf_outputs(fit, output_dir)
required_files <- c(
  "decomposition.csv", "distribution_diagnostics.csv",
  "point_resources.csv", "marginalization_diagnostics.csv",
  "run_metadata.csv", "fit.rds"
)
file_schema <- data.table::fread(file.path(schema_dir, "output_files_1.0.csv"))
stopifnot(identical(
  sort(file_schema[requirement == "required", file]),
  sort(required_files)
))
stopifnot(
  all(file.exists(file.path(output_dir, required_files))),
  identical(sort(list.files(output_dir)), sort(required_files))
)

decomposition <- data.table::fread(file.path(output_dir, "decomposition.csv"))
required_columns <- c(
  "model", "solver", "quantile", "effect", "estimate", "identified",
  "std_error", "pointwise_lower", "pointwise_upper", "uniform_lower",
  "uniform_upper", "bootstrap_reps", "bootstrap_reps_effective", "alpha"
)
decomposition_schema <- data.table::fread(
  file.path(schema_dir, "decomposition_1.0.csv")
)
stopifnot(
  identical(decomposition_schema$position, seq_along(required_columns)),
  identical(decomposition_schema$column, required_columns),
  identical(names(decomposition), required_columns),
  identical(unique(decomposition$effect), c("structure", "composition", "total")),
  is.character(decomposition$model),
  is.character(decomposition$solver),
  is.numeric(decomposition$quantile),
  is.character(decomposition$effect),
  is.numeric(decomposition$estimate),
  is.logical(decomposition$identified),
  is.integer(decomposition$bootstrap_reps),
  is.integer(decomposition$bootstrap_reps_effective),
  is.numeric(decomposition$alpha)
)
metadata <- data.table::fread(file.path(output_dir, "run_metadata.csv"))
stopifnot(
  metadata[item == "output_schema_version", value][[1L]] == "1.0",
  inherits(readRDS(file.path(output_dir, "fit.rds")), "cfdecomp")
)
unlink(output_dir, recursive = TRUE)
