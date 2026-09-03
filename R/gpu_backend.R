.cf_gpu_state <- new.env(parent = emptyenv())
.cf_gpu_state$modules <- new.env(parent = emptyenv())
.cf_gpu_state$runtime_paths <- new.env(parent = emptyenv())

gpu_python_module_path <- function() {
  path <- system.file(
    "python", "scalablecf_cuda.py", package = "scalableCounterfactual"
  )
  if (nzchar(path)) return(path)
  path <- file.path("inst", "python", "scalablecf_cuda.py")
  if (file.exists(path)) return(normalizePath(path))
  stop("scalablecf_cuda.py is missing", call. = FALSE)
}

split_environment_paths <- function(value) {
  if (!nzchar(value)) return(character())
  paths <- strsplit(value, .Platform$path.sep, fixed = TRUE)[[1L]]
  paths[nzchar(paths)]
}

prepend_unique_environment_paths <- function(name, paths) {
  paths <- paths[nzchar(paths)]
  existing <- split_environment_paths(Sys.getenv(name))
  combined <- c(paths, existing)
  key <- gsub("\\\\", "/", combined)
  if (.Platform$OS.type == "windows") key <- tolower(key)
  combined <- combined[!duplicated(key)]
  do.call(Sys.setenv, stats::setNames(
    list(paste(combined, collapse = .Platform$path.sep)), name
  ))
  invisible(combined)
}

configure_cuda_runtime_paths <- function(python_path) {
  if (is.null(python_path)) return(invisible(NULL))
  python_path <- normalizePath(python_path, winslash = "/", mustWork = TRUE)
  runtime_key <- if (.Platform$OS.type == "windows") {
    tolower(python_path)
  } else {
    python_path
  }
  if (isTRUE(.cf_gpu_state$runtime_paths[[runtime_key]])) {
    return(invisible(NULL))
  }
  prepend_unique_environment_paths("PYTHONPATH", python_path)
  nvidia <- file.path(python_path, "nvidia")
  if (dir.exists(nvidia)) {
    packages <- list.dirs(nvidia, recursive = FALSE, full.names = TRUE)
    bins <- file.path(packages, "bin")
    bins <- bins[dir.exists(bins)]
    if (length(bins)) {
      prepend_unique_environment_paths("PATH", bins)
    }
    runtime <- file.path(nvidia, "cuda_runtime")
    if (dir.exists(runtime)) Sys.setenv(CUDA_PATH = runtime)
  }
  .cf_gpu_state$runtime_paths[[runtime_key]] <- TRUE
  invisible(NULL)
}

normalize_cuda_path <- function(path, directory = FALSE) {
  if (is.null(path)) return(NULL)
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (directory && !dir.exists(normalized)) {
    stop("CUDA Python path is not a directory: ", normalized, call. = FALSE)
  }
  normalized
}

cuda_path_key <- function(path) {
  key <- gsub("\\\\", "/", path)
  if (.Platform$OS.type == "windows") tolower(key) else key
}

assert_cuda_status <- function(status) {
  if (!is.list(status) || length(status$available) != 1L ||
      !is.logical(status$available)) {
    stop("CUDA module returned an invalid capability status", call. = FALSE)
  }
  if (!isTRUE(status$available)) {
    warning_text <- paste(as.character(status$runtime_warnings %||% character()),
                          collapse = "; ")
    suffix <- if (nzchar(warning_text)) paste0(" [", warning_text, "]") else ""
    stop(
      "CUDA backend is unavailable: ", status$error %||% "unknown error", suffix,
      call. = FALSE
    )
  }
  invisible(status)
}

import_cuda_module_file <- function(module_path, module_name) {
  builtins <- reticulate::import_builtins(convert = FALSE)
  pathlib <- reticulate::import("pathlib", convert = FALSE)
  sys <- reticulate::import("sys", convert = FALSE)
  types <- reticulate::import("types", convert = FALSE)
  module_raw <- types$ModuleType(module_name)
  reticulate::py_set_attr(module_raw, "__file__", module_path)
  reticulate::py_set_attr(module_raw, "__package__", "")
  reticulate::py_set_item(sys$modules, module_name, module_raw)
  source <- pathlib$Path(module_path)$read_bytes()
  code <- builtins$compile(source, module_path, "exec")
  builtins$exec(
    code, reticulate::py_get_attr(module_raw, "__dict__")
  )
  actual_path <- reticulate::py_to_r(
    reticulate::py_get_attr(module_raw, "__file__")
  )
  actual_path <- normalizePath(actual_path, winslash = "/", mustWork = TRUE)
  if (!identical(cuda_path_key(actual_path), cuda_path_key(module_path))) {
    stop(
      "CUDA module identity mismatch: requested ", module_path,
      " but Python loaded ", actual_path,
      call. = FALSE
    )
  }
  reticulate::import(module_name, convert = TRUE)
}

load_cuda_module <- function(
    python = NULL, python_path = NULL, module_path = NULL,
    require_available = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("CUDA requires the optional R package 'reticulate'", call. = FALSE)
  }
  python <- normalize_cuda_path(python)
  python_path <- normalize_cuda_path(python_path, directory = TRUE)
  if (is.null(module_path)) module_path <- gpu_python_module_path()
  module_path <- normalize_cuda_path(module_path)
  module_content_md5 <- unname(tools::md5sum(module_path))
  identity <- list(
    python = python %||% "<default>",
    python_path = python_path %||% "<default-path>",
    module_path = cuda_path_key(module_path),
    module_content_md5 = module_content_md5
  )
  key <- object_md5(identity)
  cached <- .cf_gpu_state$modules[[key]]
  if (!is.null(cached)) {
    if (isTRUE(require_available)) assert_cuda_status(cached$cuda_status())
    return(cached)
  }
  configure_cuda_runtime_paths(python_path)
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)
  if (!is.null(python_path)) {
    sys <- reticulate::import("sys", convert = FALSE)
    existing <- reticulate::py_to_r(sys$path)
    existing_key <- vapply(existing, function(value) {
      cuda_path_key(gsub("\\\\", "/", as.character(value)))
    }, character(1L))
    if (!cuda_path_key(python_path) %in% existing_key) {
      sys$path$insert(0L, python_path)
    }
  }
  module_name <- paste0(
    "scalablecf_cuda_", object_md5(list(
      path = cuda_path_key(module_path), content_md5 = module_content_md5
    ))
  )
  module <- import_cuda_module_file(module_path, module_name)
  status <- module$cuda_status()
  .cf_gpu_state$modules[[key]] <- module
  if (isTRUE(require_available)) assert_cuda_status(status)
  module
}

#' Inspect the optional GPU computation backend
#'
#' @param python Optional Python executable containing CuPy.
#' @param python_path Optional project-local Python package directory.
#' @param module_path Optional development override for the CUDA module.
#' @return CUDA availability and runtime metadata.
#' @export
gpu_backend_status <- function(
    python = NULL, python_path = NULL, module_path = NULL) {
  tryCatch({
    module <- load_cuda_module(
      python, python_path, module_path, require_available = FALSE
    )
    status <- module$cuda_status()
    status$python <- status$python_executable %||% reticulate::py_config()$python
    status
  }, error = function(condition) list(
    available = FALSE, device_id = NULL, device = NULL,
    cupy_version = NULL, numpy_version = NULL,
    python = python, python_version = NULL,
    module_file = module_path, module_sha256 = NULL,
    capability_matmul = FALSE, capability_sort = FALSE,
    capability_solve = FALSE, runtime_warnings = character(),
    error = conditionMessage(condition), error_type = class(condition)[[1L]]
  ))
}

gpu_runtime_metadata <- function(
    control = NULL, python = NULL, python_path = NULL, module_path = NULL) {
  backend <- "cuda"
  prediction_backend <- NA_character_
  dr_backend <- NA_character_
  precision <- NA_character_
  if (!is.null(control)) {
    prediction_backend <- control$gpu_backend %||% NA_character_
    dr_backend <- control$dr_backend %||% NA_character_
    precision <- control$gpu_precision %||% NA_character_
    if (identical(prediction_backend, "cuda") ||
        identical(dr_backend, "cuda")) backend <- "cuda"
    python <- python %||% control$gpu_python
    python_path <- python_path %||% control$gpu_python_path
    module_path <- module_path %||% control$gpu_module_path
  }
  status <- gpu_backend_status(python, python_path, module_path)
  list(
    backend = backend,
    prediction_backend = prediction_backend,
    dr_backend = dr_backend,
    precision = precision,
    available = isTRUE(status$available),
    python = status$python %||% python %||% NA_character_,
    python_version = status$python_version %||% NA_character_,
    cupy_version = status$cupy_version %||% NA_character_,
    numpy_version = status$numpy_version %||% NA_character_,
    device_id = status$device_id %||% NA_integer_,
    device = status$device %||% NA_character_,
    module_file = status$module_file %||% module_path %||% NA_character_,
    module_sha256 = status$module_sha256 %||% NA_character_,
    module_hash = status$module_sha256 %||% NA_character_,
    cupy_file = status$cupy_file %||% NA_character_,
    numpy_file = status$numpy_file %||% NA_character_,
    cuda_runtime_version = status$cuda_runtime_version %||% NA_integer_,
    cuda_driver_version = status$cuda_driver_version %||% NA_integer_,
    capability_matmul = isTRUE(status$capability_matmul),
    capability_sort = isTRUE(status$capability_sort),
    capability_solve = isTRUE(status$capability_solve),
    python_isolated = status$python_isolated %||% NA,
    python_no_user_site = status$python_no_user_site %||% NA,
    user_site_enabled = status$user_site_enabled %||% NA,
    cupy_distributions = as.character(
      status$cupy_distributions %||% character()
    ),
    runtime_warnings = as.character(status$runtime_warnings %||% character()),
    error = status$error %||% NA_character_,
    error_type = status$error_type %||% NA_character_,
    reticulate_version = if (requireNamespace("reticulate", quietly = TRUE)) {
      as.character(utils::packageVersion("reticulate"))
    } else {
      NA_character_
    }
  )
}

gpu_run_metadata_fields <- function(control) {
  metadata <- gpu_runtime_metadata(control = control)
  list(
    gpu_runtime_available = metadata$available,
    gpu_python_executable = metadata$python,
    gpu_python_version = metadata$python_version,
    gpu_numpy_version = metadata$numpy_version,
    gpu_numpy_file = metadata$numpy_file,
    gpu_cupy_version = metadata$cupy_version,
    gpu_cupy_file = metadata$cupy_file,
    gpu_device_id = metadata$device_id,
    gpu_device_name = metadata$device,
    gpu_cuda_runtime_version = metadata$cuda_runtime_version,
    gpu_cuda_driver_version = metadata$cuda_driver_version,
    gpu_module_file = metadata$module_file,
    gpu_module_sha256 = metadata$module_sha256,
    gpu_capability_matmul = metadata$capability_matmul,
    gpu_capability_sort = metadata$capability_sort,
    gpu_capability_solve = metadata$capability_solve,
    gpu_python_isolated = metadata$python_isolated,
    gpu_python_no_user_site = metadata$python_no_user_site,
    gpu_user_site_enabled = metadata$user_site_enabled,
    gpu_cupy_distributions = paste(
      metadata$cupy_distributions, collapse = " | "
    ),
    gpu_runtime_warnings = paste(metadata$runtime_warnings, collapse = " | "),
    gpu_runtime_error_type = metadata$error_type,
    gpu_runtime_error = metadata$error,
    gpu_reticulate_version = metadata$reticulate_version
  )
}

gpu_module_from_control <- function(control) {
  load_cuda_module(
    control$gpu_python, control$gpu_python_path, control$gpu_module_path
  )
}

gpu_matmul <- function(X, coefficients, control) {
  coefficients <- as.matrix(coefficients)
  result <- gpu_module_from_control(control)$matmul(
    unname(as.matrix(X)), unname(coefficients), control$gpu_precision
  )
  matrix(as.numeric(result), nrow = nrow(X), ncol = ncol(coefficients))
}

gpu_qr_marginal_quantiles <- function(
    X, coefficients, weights, probs, shift, control,
    normalization_rows = nrow(X)) {
  if (!is.numeric(normalization_rows) || length(normalization_rows) != 1L ||
      !is.finite(normalization_rows) || normalization_rows < 1) {
    stop("normalization_rows must be one finite value >= 1", call. = FALSE)
  }
  weights <- normalize_weights(weights) *
    (as.numeric(normalization_rows) / length(weights))
  result <- gpu_module_from_control(control)$qr_marginal_quantiles(
    unname(as.matrix(X)), unname(as.matrix(coefficients)),
    as.numeric(weights), as.numeric(probs), as.numeric(shift),
    control$gpu_precision, as.numeric(normalization_rows)
  )
  as.numeric(result)
}

gpu_dr_marginal_cdf <- function(X, coefficients, weights, link, control) {
  as.numeric(gpu_module_from_control(control)$dr_marginal_cdf(
    unname(as.matrix(X)), unname(as.matrix(coefficients)), as.numeric(weights),
    link, control$gpu_precision, as.integer(control$gpu_block_columns)
  ))
}

fit_dr_process_cuda <- function(
    X, y, weights, thresholds, model, maxit, tolerance,
    constant_direction, control) {
  constant_zero <- thresholds < min(y)
  constant_one <- thresholds >= max(y)
  regular <- !(constant_zero | constant_one)
  p <- ncol(X)
  coefficients <- matrix(NA_real_, p, length(thresholds))
  iterations <- integer(length(thresholds))
  converged <- rep(TRUE, length(thresholds))
  boundary <- rep(FALSE, length(thresholds))
  link_value <- function(probability) switch(
    model,
    logit = stats::qlogis(probability),
    probit = stats::qnorm(probability),
    cloglog = log(-log1p(-probability))
  )
  if (any(constant_zero)) coefficients[, constant_zero] <-
    constant_direction * link_value(1e-8)
  if (any(constant_one)) coefficients[, constant_one] <-
    constant_direction * link_value(1 - 1e-8)
  elapsed <- 0
  if (any(regular)) {
    returned <- gpu_module_from_control(control)$fit_dr_process(
      unname(as.matrix(X)), as.numeric(y), as.numeric(weights),
      as.numeric(thresholds[regular]), model, control$gpu_precision,
      as.integer(maxit), as.numeric(tolerance),
      as.integer(control$gpu_block_columns)
    )
    coefficients[, regular] <- matrix(
      as.numeric(returned$coefficients), p, sum(regular)
    )
    iterations[regular] <- as.integer(returned$iterations)
    converged[regular] <- as.logical(returned$converged)
    boundary[regular] <- as.logical(returned$boundary)
    elapsed <- as.numeric(returned$elapsed_seconds)
  }
  list(
    coefficients = coefficients, iterations = iterations,
    converged = converged, boundary = boundary, elapsed_seconds = elapsed,
    backend = ifelse(regular, "cuda", "analytic"),
    fallback_reason = rep(NA_character_, length(thresholds))
  )
}

fit_qr_cuda_admm <- function(X, y, weights, taus, options) {
  control <- options$control
  preconditioned <- precondition_qr_design(
    X * normalize_weights(weights), isTRUE(options$precondition)
  )
  transform <- preconditioned$transform
  fit_X <- X %*% transform
  returned <- gpu_module_from_control(control)$fit_qr_admm(
    unname(as.matrix(fit_X)), as.numeric(y), normalize_weights(weights),
    as.numeric(taus), control$gpu_precision, as.numeric(control$gpu_qr_rho),
    as.integer(control$gpu_qr_maxit), as.numeric(control$gpu_qr_tolerance),
    as.integer(control$gpu_block_columns)
  )
  coefficients <- transform %*% matrix(
    as.numeric(returned$coefficients), ncol(X), length(taus)
  )
  conditioned_coefficients <- matrix(
    as.numeric(returned$coefficients), ncol(X), length(taus)
  )
  converged <- as.logical(returned$converged)
  if (any(!converged) && !isTRUE(control$gpu_qr_allow_nonconvergence)) {
    stop(
      "cuda_admm reached gpu_qr_maxit without convergence at tau(s): ",
      paste(signif(taus[!converged], 5), collapse = ", "),
      "; increase gpu_qr_maxit or use an exact CPU solver",
      call. = FALSE
    )
  }
  list(
    coefficients = coefficients,
    conditioned_coefficients = conditioned_coefficients,
    flag = ifelse(converged, 0L, 1L),
    nit = as.integer(returned$iterations),
    cuda_converged = converged,
    primal_residual = as.numeric(returned$primal_residual),
    dual_residual = as.numeric(returned$dual_residual),
    precondition_requested = isTRUE(options$precondition),
    preconditioned = preconditioned$method != "none",
    preconditioning_method = preconditioned$method,
    preconditioning_gram_rcond = preconditioned$gram_rcond,
    preconditioning_condition_estimate = preconditioned$condition_estimate,
    preconditioning_fallback_reason = preconditioned$fallback_reason,
    preconditioning_matrix = transform,
    warnings = if (all(converged)) character() else
      "cuda_admm reached its iteration limit for at least one quantile"
  )
}
