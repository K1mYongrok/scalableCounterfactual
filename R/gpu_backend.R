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

load_cuda_module <- function(
    python = NULL, python_path = NULL, module_path = NULL) {
  key <- paste(
    python %||% "<default>", python_path %||% "<default-path>",
    module_path %||% "<bundled>", sep = "|"
  )
  if (!is.null(.cf_gpu_state$modules[[key]])) {
    return(.cf_gpu_state$modules[[key]])
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("CUDA requires the optional R package 'reticulate'", call. = FALSE)
  }
  configure_cuda_runtime_paths(python_path)
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)
  if (!is.null(python_path)) {
    reticulate::py_run_string(sprintf(
      "import sys; p=r'''%s'''; sys.path.insert(0,p) if p not in sys.path else None",
      normalizePath(python_path, winslash = "/", mustWork = TRUE)
    ))
  }
  if (is.null(module_path)) module_path <- gpu_python_module_path()
  module_path <- normalizePath(module_path, winslash = "/", mustWork = TRUE)
  module <- reticulate::import_from_path(
    tools::file_path_sans_ext(basename(module_path)),
    path = dirname(module_path), convert = TRUE
  )
  status <- module$cuda_status()
  if (!isTRUE(status$available)) {
    stop("CUDA backend is unavailable: ", status$error, call. = FALSE)
  }
  .cf_gpu_state$modules[[key]] <- module
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
    module <- load_cuda_module(python, python_path, module_path)
    c(module$cuda_status(), list(python = reticulate::py_config()$python))
  }, error = function(condition) list(
    available = FALSE, device = NULL, cupy_version = NULL,
    python = python, error = conditionMessage(condition)
  ))
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
    X, coefficients, weights, probs, shift, control) {
  result <- gpu_module_from_control(control)$qr_marginal_quantiles(
    unname(as.matrix(X)), unname(as.matrix(coefficients)),
    as.numeric(weights), as.numeric(probs), as.numeric(shift),
    control$gpu_precision
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
