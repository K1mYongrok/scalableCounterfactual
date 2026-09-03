library(scalableCounterfactual)

stopifnot(utils::packageVersion("quantreg") >= "5.61")

common_path <- system.file(
  "scripts", "_cli_common.R", package = "scalableCounterfactual"
)
stopifnot(nzchar(common_path), file.exists(common_path))
common <- new.env(parent = baseenv())
sys.source(common_path, envir = common)
stopifnot(
  isTRUE(common$as_bool(" TRUE ")),
  isTRUE(common$as_bool("y")),
  identical(common$as_bool("0"), FALSE),
  inherits(try(common$as_bool("sometimes"), silent = TRUE), "try-error"),
  inherits(try(common$as_bool(NA_character_), silent = TRUE), "try-error"),
  inherits(try(common$as_bool(c("true", "false")), silent = TRUE), "try-error")
)
stopifnot(
  inherits(try(common$assert_distinct_cli_paths(list(
    input = "same.csv", output = "SAME.csv"
  )), silent = TRUE), "try-error") == (.Platform$OS.type == "windows"),
  inherits(try(common$assert_cli_directories_not_below_files(
    list(checkpoint = file.path("out", "summary.csv", "checkpoints")),
    list(summary = file.path("out", "summary.csv"))
  ), silent = TRUE), "try-error")
)

benchmark_path <- system.file(
  "scripts", "cfbenchmark.R", package = "scalableCounterfactual"
)
decomp_path <- system.file(
  "scripts", "cfdecomp.R", package = "scalableCounterfactual"
)
scaling_path <- system.file(
  "scripts", "cfscaling.R", package = "scalableCounterfactual"
)
simulation_path <- system.file(
  "scripts", "cfsimulate.R", package = "scalableCounterfactual"
)
benchmark_text <- paste(readLines(benchmark_path, warn = FALSE), collapse = "\n")
stopifnot(
  grepl('value\\(args, "repetitions", 5L, as.integer\\)', benchmark_text),
  grepl('value\\(args, "warmup", 1L, as.integer\\)', benchmark_text)
)

run_cli <- function(script, arguments) {
  executable <- file.path(
    R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  output <- suppressWarnings(system2(
    executable,
    c("--vanilla", shQuote(script), arguments),
    stdout = TRUE, stderr = TRUE
  ))
  list(
    output = output,
    status = if (is.null(attr(output, "status"))) 0L else attr(output, "status")
  )
}

# Destructive CLI path aliases are rejected before any input or checkpoint is
# overwritten. The hash checks make this an end-to-end data-loss regression.
cli_root <- tempfile("scalablecf-cli-collision-")
dir.create(cli_root)
cli_data <- file.path(cli_root, "decomposition.csv")
utils::write.csv(data.frame(
  y = 1:6, x = 1:6, group = rep(0:1, each = 3L), weight = 1
), cli_data, row.names = FALSE)
cli_hash <- unname(tools::md5sum(cli_data))
common_cli <- c(
  "--data", shQuote(cli_data), "--formula", shQuote("y ~ x"),
  "--group", "group", "--weights", "weight"
)
benchmark_collision <- run_cli(
  benchmark_path, c(common_cli, "--output", shQuote(cli_data))
)
scaling_collision <- run_cli(
  scaling_path, c(
    common_cli, "--sample-sizes", "6", "--output",
    shQuote(file.path(cli_root, "scaling.csv")),
    "--checkpoint", shQuote(cli_data), "--resume", "false"
  )
)
decomp_collision <- run_cli(
  decomp_path, c(common_cli, "--output", shQuote(cli_root))
)
simulation_collision <- run_cli(
  simulation_path, c(
    "--replications", "1", "--n-per-group", "20",
    "--output", shQuote(cli_root),
    "--checkpoint-dir", shQuote(file.path(cli_root, "summary.csv", "nested"))
  )
)
stopifnot(
  benchmark_collision$status != 0L,
  scaling_collision$status != 0L,
  decomp_collision$status != 0L,
  simulation_collision$status != 0L,
  identical(unname(tools::md5sum(cli_data)), cli_hash),
  !dir.exists(file.path(cli_root, "summary.csv"))
)
unlink(cli_root, recursive = TRUE, force = TRUE)

decomp_text <- paste(readLines(decomp_path, warn = FALSE), collapse = "\n")
required_cli_arguments <- c(
  "robust_se", "dr_noncrossing", "gpu_backend", "gpu_precision", "gpu_python",
  "gpu_python_path", "gpu_module_path", "gpu_block_columns", "gpu_qr_rho",
  "gpu_qr_maxit", "gpu_qr_tolerance", "gpu_qr_allow_nonconvergence"
)
stopifnot(all(vapply(
  required_cli_arguments,
  function(argument) grepl(argument, decomp_text, fixed = TRUE),
  logical(1L)
)))
stopifnot(
  grepl(
    'value\\(args, "legacy_qr_shift", TRUE, as_bool\\)',
    decomp_text
  ),
  grepl(
    'value\\(args, "dr_noncrossing", "cummax"\\)',
    decomp_text
  )
)

fastglm_namespace_before <- "fastglm" %in% loadedNamespaces()
invisible(scalableCounterfactual:::bootstrap_runtime_identity(
  "qr", "fn", cf_control(dr_backend = "auto")
))
stopifnot(identical(
  "fastglm" %in% loadedNamespaces(), fastglm_namespace_before
))

auto_dr_control <- cf_control(dr_backend = "auto")
auto_dr_runtime <- scalableCounterfactual:::bootstrap_runtime_identity(
  "logit", NULL, auto_dr_control
)
stopifnot(
  "fastglm" %in% names(auto_dr_runtime$packages),
  identical(
    auto_dr_runtime$backend_availability$fastglm,
    requireNamespace("fastglm", quietly = TRUE)
  )
)

review_backend <- "reviewrequested"
register_conditional_backend(
  review_backend, type = "distribution",
  fit = function(X, response, weights, model, start, maxit, tolerance) {
    list(coefficients = rep(0, ncol(X)))
  },
  version = "1"
)
review_control <- cf_control(dr_backend = review_backend)
review_point <- list(fits = list(
  group0 = list(solver = "fn", selection_backend = "glm"),
  group1 = list(solver = "fn", selection_backend = "glm")
))
review_fingerprint <- scalableCounterfactual:::active_extension_fingerprint(
  "cqr", "fn", review_control, review_point
)
expected_fingerprint <- scalableCounterfactual:::extension_registry_fingerprint(
  list(
    qr = "fn", linear = character(),
    distribution = c(review_backend, "glm")
  )
)
unregister_conditional_backend(review_backend, type = "distribution")
stopifnot(identical(review_fingerprint, expected_fingerprint))

cpu_namespace_before <- "reticulate" %in% loadedNamespaces()
cpu_point_data <- data.frame(
  y = seq_len(80L) / 10 + stats::rnorm(80L),
  x = stats::rnorm(80L), group = rep(0:1, each = 40L),
  weight = rep(1, 80L)
)
cpu_point <- counterfactual_decompose(
  y ~ x, cpu_point_data, "group", "weight", model = "loc",
  control = cf_control(
    nreg = 5L, reported_quantiles = c(0.25, 0.5, 0.75),
    crossing_diagnostics = FALSE, gpu_backend = "cpu"
  )
)

# Group fitting has at most two independent tasks. Benchmarks report both the
# request and the effective worker count instead of overstating parallelism.
worker_benchmark <- benchmark_qr_solvers(
  y ~ x, cpu_point_data, "group", "weight",
  solvers = "fn", reference_solver = "fn",
  control = cf_control(
    nreg = 3L, reported_quantiles = 0.5, crossing_diagnostics = FALSE
  ),
  point_workers = 7L
)
stopifnot(
  identical(worker_benchmark$point_workers_requested, 7L),
  identical(worker_benchmark$point_workers, 2L)
)
cuda_fit_error <- tryCatch(benchmark_qr_solvers(
  y ~ x, cpu_point_data, "group", "weight",
  solvers = "cuda_admm", reference_solver = "cuda_admm",
  control = cf_control(crossing_diagnostics = FALSE),
  point_workers = 2L
), error = identity)
cuda_dr_error <- tryCatch(benchmark_conditional_backends(
  y ~ x, cpu_point_data, "group", "weight", model = "logit",
  backends = "cuda", reference_backend = "cuda",
  control = cf_control(crossing_diagnostics = FALSE),
  point_workers = 2L
), error = identity)
prediction_only <- scalableCounterfactual:::validate_execution_parallelism(
  "qr", "fn", cf_control(gpu_backend = "cuda"), 2L
)
qr_dr_worker_scope <- scalableCounterfactual:::validate_execution_parallelism(
  "qr", "fn", cf_control(dr_workers = 2L), 2L
)
stopifnot(
  inherits(cuda_fit_error, "error"),
  grepl("CUDA conditional fitting", conditionMessage(cuda_fit_error)),
  inherits(cuda_dr_error, "error"),
  grepl("CUDA conditional fitting", conditionMessage(cuda_dr_error)),
  identical(prediction_only$point_workers_effective, 2L),
  identical(qr_dr_worker_scope$point_workers_effective, 2L)
)
cuda_simulation_error <- tryCatch(simulate_counterfactual_validation(
  replications = 2L, n_per_group = 20L, scenarios = "location_shift",
  solver = "cuda_admm", workers = 2L, progress = FALSE
), error = identity)
stopifnot(
  inherits(cuda_simulation_error, "error"),
  grepl("requires one effective worker", conditionMessage(cuda_simulation_error))
)
worker_simulation <- simulate_counterfactual_validation(
  replications = 1L, n_per_group = 20L, scenarios = "location_shift",
  solver = "fn", workers = 20L, progress = FALSE
)
stopifnot(
  identical(worker_simulation$settings$workers_requested, 20L),
  identical(worker_simulation$settings$workers, 1L)
)
stopifnot(
  identical("reticulate" %in% loadedNamespaces(), cpu_namespace_before),
  !any(startsWith(names(cpu_point$metadata), "gpu_runtime_")),
  !any(startsWith(names(cpu_point$metadata), "gpu_python_"))
)

configured_gpu_python <- Sys.getenv("SCALABLECF_GPU_PYTHON", unset = "")
configured_gpu_path <- Sys.getenv("SCALABLECF_GPU_PATH", unset = "")
if (nzchar(configured_gpu_python) && nzchar(configured_gpu_path)) {
  status <- gpu_backend_status(configured_gpu_python, configured_gpu_path)
  stopifnot(
    isTRUE(status$available),
    isTRUE(status$capability_matmul),
    isTRUE(status$capability_sort),
    isTRUE(status$capability_solve),
    nzchar(status$module_file),
    grepl("^[0-9a-f]{64}$", status$module_sha256)
  )
  float32_control <- cf_control(
    dr_backend = "cuda", gpu_backend = "cuda", gpu_precision = "float32",
    gpu_python = configured_gpu_python,
    gpu_python_path = configured_gpu_path
  )
  module <- scalableCounterfactual:::gpu_module_from_control(float32_control)
  predictor <- seq(-8, 8, length.out = 100L)
  float32_fit <- module$fit_dr_process(
    cbind(1, predictor), predictor, rep(1, length(predictor)),
    as.array(0), "logit", "float32", 10L, 1e-6, 1L
  )
  stopifnot(isTRUE(as.logical(float32_fit$boundary[[1L]])))

  # A length-one threshold can cross the reticulate boundary as a scalar.  The
  # CUDA process must preserve a one-column result instead of indexing a 0-D
  # CuPy array.
  single_threshold_fit <- module$fit_dr_process(
    cbind(1, predictor), predictor, rep(1, length(predictor)),
    0, "logit", "float64", 100L, 1e-8, 1L
  )
  stopifnot(
    identical(dim(single_threshold_fit$coefficients), c(2L, 1L)),
    length(single_threshold_fit$converged) == 1L,
    length(single_threshold_fit$boundary) == 1L
  )

  set.seed(20260904L)
  latent <- stats::rnorm(800L)
  response <- stats::rbinom(800L, 1L, stats::plogis(0.6 * latent))
  scaled_fit <- module$fit_dr_process(
    cbind(1, latent / 1000), ifelse(response == 1L, -1, 1),
    rep(1, length(latent)), as.array(0), "logit", "float64",
    100L, 1e-8, 1L
  )
  stopifnot(
    abs(scaled_fit$coefficients[2L, 1L]) > 100,
    !isTRUE(as.logical(scaled_fit$boundary[[1L]]))
  )

  # Extreme weights exercise the shared rank-boundary tolerance.  The CUDA
  # cumulative sum is float64 even when prediction storage may be float32.
  set.seed(128L)
  parity_n <- 37L
  parity_X <- cbind(1, stats::rnorm(parity_n), stats::runif(parity_n, -2, 2))
  parity_coefficients <- matrix(stats::rnorm(15L), nrow = 3L, ncol = 5L)
  parity_weights <- 10^stats::runif(parity_n, -16, 16)
  parity_probs <- c(0.1, 0.33, 0.5, 0.77, 0.9)
  parity_draws <- as.vector(parity_X %*% parity_coefficients)
  parity_cpu <- scalableCounterfactual:::weighted_quantile(
    parity_draws,
    rep(parity_weights, times = ncol(parity_coefficients)),
    parity_probs
  )
  parity_control <- cf_control(
    gpu_backend = "cuda", gpu_precision = "float64",
    gpu_python = configured_gpu_python,
    gpu_python_path = configured_gpu_path
  )
  parity_gpu <- scalableCounterfactual:::gpu_qr_marginal_quantiles(
    parity_X, parity_coefficients, parity_weights, parity_probs, 0,
    parity_control
  )
  stopifnot(max(abs(parity_cpu - parity_gpu)) < 1e-10)

  frequency_X <- cbind(1, c(0, 10))
  frequency_coefficients <- matrix(c(0, 1), nrow = 2L)
  frequency_cpu <- scalableCounterfactual:::weighted_quantile(
    c(0, 10), c(2, 1), 0.5, normalization_n = 3
  )
  frequency_gpu <- scalableCounterfactual:::gpu_qr_marginal_quantiles(
    frequency_X, frequency_coefficients, c(2, 1), 0.5, 0,
    parity_control, normalization_rows = 3
  )
  stopifnot(identical(frequency_gpu, frequency_cpu))

  # CUDA ADMM is an iterative row-level algorithm. Non-unit multiplicities are
  # expanded internally so a compressed bootstrap draw is exactly the same
  # numerical problem as its explicit representation.
  set.seed(831L)
  admm_X <- cbind(1, stats::rnorm(28L), stats::runif(28L, -1, 1))
  admm_y <- drop(admm_X %*% c(0.3, -0.5, 0.2) + stats::rnorm(28L, sd = 0.2))
  admm_w <- stats::runif(28L, 0.5, 2)
  admm_frequency <- rep(c(1L, 2L, 3L, 1L), length.out = 28L)
  admm_index <- rep.int(seq_len(nrow(admm_X)), admm_frequency)
  admm_control <- cf_control(
    gpu_backend = "cuda", gpu_precision = "float64",
    gpu_python = configured_gpu_python,
    gpu_python_path = configured_gpu_path,
    gpu_qr_maxit = 5000L, gpu_qr_tolerance = 1e-7,
    gpu_qr_allow_nonconvergence = TRUE
  )
  admm_collapsed <- suppressWarnings(fit_weighted_qr(
    admm_X, admm_y, admm_w, c(0.25, 0.5, 0.75),
    solver = "cuda_admm", gpu_control = admm_control,
    frequency = admm_frequency
  ))
  admm_expanded <- suppressWarnings(fit_weighted_qr(
    admm_X[admm_index, , drop = FALSE], admm_y[admm_index],
    admm_w[admm_index], c(0.25, 0.5, 0.75),
    solver = "cuda_admm", gpu_control = admm_control
  ))
  stopifnot(
    identical(admm_collapsed$frequency_expanded, TRUE),
    identical(admm_expanded$frequency_expanded, FALSE),
    identical(
      admm_collapsed$frequency_effective_n,
      admm_expanded$frequency_effective_n
    ),
    isTRUE(all.equal(
      admm_collapsed$coefficients, admm_expanded$coefficients,
      tolerance = 0, check.attributes = FALSE
    ))
  )

  gpu_point <- counterfactual_decompose(
    y ~ x, cpu_point_data, "group", "weight", model = "loc",
    control = cf_control(
      nreg = 5L, reported_quantiles = c(0.25, 0.5, 0.75),
      crossing_diagnostics = FALSE, gpu_backend = "cuda",
      gpu_precision = "float64", gpu_python = configured_gpu_python,
      gpu_python_path = configured_gpu_path
    )
  )
  required_runtime_fields <- c(
    "gpu_python_executable", "gpu_python_version", "gpu_numpy_version",
    "gpu_cupy_version", "gpu_device_id", "gpu_device_name",
    "gpu_cuda_runtime_version", "gpu_cuda_driver_version",
    "gpu_module_file", "gpu_module_sha256", "gpu_capability_matmul",
    "gpu_capability_sort", "gpu_capability_solve", "gpu_python_isolated",
    "gpu_python_no_user_site", "gpu_user_site_enabled",
    "gpu_runtime_warnings"
  )
  stopifnot(
    all(required_runtime_fields %in% names(gpu_point$metadata)),
    all(vapply(
      gpu_point$metadata[required_runtime_fields], length, integer(1L)
    ) == 1L),
    grepl("^[0-9a-f]{64}$", gpu_point$metadata$gpu_module_sha256),
    isTRUE(gpu_point$metadata$gpu_capability_matmul),
    isTRUE(gpu_point$metadata$gpu_capability_sort),
    isTRUE(gpu_point$metadata$gpu_capability_solve)
  )
  metadata_output <- tempfile("scalablecf-gpu-metadata-")
  write_cf_outputs(gpu_point, metadata_output)
  written_metadata <- data.table::fread(
    file.path(metadata_output, "run_metadata.csv")
  )
  stopifnot(all(required_runtime_fields %in% written_metadata$item))
  unlink(metadata_output, recursive = TRUE, force = TRUE)
}

find_test_python <- function() {
  configured <- Sys.getenv("SCALABLECF_GPU_PYTHON", unset = "")
  candidates <- c(configured, Sys.which("python3"), Sys.which("python"))
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  for (candidate in candidates) {
    status <- suppressWarnings(system2(
      candidate, c("-c", shQuote("import sys; print(sys.executable)")),
      stdout = TRUE, stderr = FALSE
    ))
    if (length(status) && is.null(attr(status, "status")) && file.exists(status[[1L]])) {
      return(status[[1L]])
    }
  }
  launcher <- Sys.which("py")
  if (nzchar(launcher)) {
    status <- suppressWarnings(system2(
      launcher, c("-3", "-c", shQuote("import sys; print(sys.executable)")),
      stdout = TRUE, stderr = FALSE
    ))
    if (length(status) && is.null(attr(status, "status")) && file.exists(status[[1L]])) {
      return(status[[1L]])
    }
  }
  NULL
}

test_python <- find_test_python()
if (requireNamespace("reticulate", quietly = TRUE) && !is.null(test_python) &&
    isTRUE(tryCatch(
      {
        reticulate::use_python(test_python, required = TRUE)
        reticulate::py_available(initialize = TRUE)
      },
      error = function(error) FALSE
    ))) {
  make_fake_module <- function(path, marker, available = TRUE) {
    available_python <- if (available) "True" else "False"
    error_python <- if (available) "None" else "'fake unavailable'"
    writeLines(c(
      sprintf("MARKER = %s", as.integer(marker)),
      "def marker():",
      "    return MARKER",
      "def cuda_status(refresh=False):",
      "    return {",
      sprintf("        'available': %s,", available_python),
      sprintf("        'error': %s,", error_python),
      "        'error_type': None,",
      "        'python_executable': 'fake-python',",
      "        'python_version': '3.test',",
      "        'numpy_version': '2.test',",
      "        'cupy_version': '14.test',",
      "        'device_id': 0,",
      "        'device': 'fake-device',",
      "        'module_file': __file__,",
      sprintf("        'module_sha256': 'fake-sha-%s',", as.integer(marker)),
      "        'capability_matmul': True,",
      "        'capability_sort': True,",
      "        'capability_solve': True,",
      "        'runtime_warnings': []",
      "    }"
    ), path, useBytes = TRUE)
  }
  root <- tempfile("scalablecf-module-test-")
  first_dir <- file.path(root, "first")
  second_dir <- file.path(root, "second")
  dir.create(first_dir, recursive = TRUE)
  dir.create(second_dir, recursive = TRUE)
  first <- file.path(first_dir, "same_name.py")
  second <- file.path(second_dir, "same_name.py")
  make_fake_module(first, 1L)
  make_fake_module(second, 2L)

  first_module <- scalableCounterfactual:::load_cuda_module(
    module_path = first
  )
  second_module <- scalableCounterfactual:::load_cuda_module(
    module_path = second
  )
  stopifnot(
    identical(as.integer(first_module$marker()), 1L),
    identical(as.integer(second_module$marker()), 2L),
    !identical(first_module, second_module)
  )

  make_fake_module(first, 3L)
  changed_module <- scalableCounterfactual:::load_cuda_module(
    module_path = first
  )
  changed_status <- scalableCounterfactual:::gpu_runtime_metadata(
    module_path = first
  )
  stopifnot(
    identical(as.integer(changed_module$marker()), 3L),
    isTRUE(changed_status$available),
    identical(changed_status$module_sha256, "fake-sha-3"),
    isTRUE(changed_status$capability_matmul),
    isTRUE(changed_status$capability_sort),
    isTRUE(changed_status$capability_solve),
    identical(
      tolower(normalizePath(changed_status$module_file, winslash = "/")),
      tolower(normalizePath(first, winslash = "/"))
    )
  )
  unlink(root, recursive = TRUE, force = TRUE)
} else {
  message("Skipping module-identity test: reticulate Python is unavailable")
}

cat("GPU, CLI, and packaging review-fix tests passed.\n")
