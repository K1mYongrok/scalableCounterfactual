supported_cf_models <- function() {
  c("qr", "cqr", "loc", "locsca", "cox", "logit", "probit", "cloglog", "lpm")
}

is_quantile_process_model <- function(model) model %in% c("qr", "cqr")

is_distribution_regression_model <- function(model) {
  model %in% c("logit", "probit", "cloglog", "lpm")
}

.datatable.aware <- TRUE
.cf_output_schema_version <- "1.0"

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

utils::globalVariables(c(
  ".", ".N", "bin", "effect", "error", "estimate", "pointwise_covered",
  "quantile", "replication", "scenario", "squared_error", "truth",
  "uniform_curve_covered", "value", "weight"
))

builtin_qr_solvers <- function() {
  c(
    "br", "fn", "pfn", "qfnb", "pfnb",
    "proqreg", "profn", "onestep", "cuda_admm", "auto"
  )
}

supported_qr_solvers <- function() {
  c(builtin_qr_solvers(), names(.cf_extension_registry$qr))
}

builtin_linear_backends <- function() c("auto", "qr", "chol", "fastglm")
builtin_dr_backends <- function() {
  c("auto", "glm", "fastglm", "speedglm", "cuda")
}

supported_linear_backends <- function() {
  c(builtin_linear_backends(), names(.cf_extension_registry$linear))
}

supported_dr_backends <- function() {
  c(builtin_dr_backends(), names(.cf_extension_registry$distribution))
}

#' Describe conditional-model computation backends
#'
#' @return A data frame listing models, backend option names, implementations,
#'   and whether they preserve the corresponding unpenalized model objective.
#' @export
conditional_backend_registry <- function() {
  builtins <- data.frame(
    model = c(
      rep("loc/locsca/lpm", 4L),
      rep("logit/probit/cloglog", 5L)
    ),
    backend = c(
      "auto", "qr", "chol", "fastglm",
      "auto", "glm", "fastglm", "speedglm", "cuda"
    ),
    implementation = c(
      "base R QR (robust default)",
      "base R reusable QR factorization",
      "base R reusable Cholesky factorization",
      "fastglm::fastglmPure Gaussian WLS (RcppEigen)",
      "fastglm when installed, otherwise stats::glm.fit",
      "stats::glm.fit",
      "fastglm::fastglmPure binomial IRLS (RcppEigen)",
      "speedglm::speedglm.wfit",
      "CuPy batched unpenalized binomial IRLS"
    ),
    objective_preserving = TRUE,
    optional_package = c(
      NA, NA, NA, "fastglm", NA, NA, "fastglm", "speedglm", "reticulate/CuPy"
    ),
    stringsAsFactors = FALSE
  )
  custom_linear <- lapply(.cf_extension_registry$linear, function(entry) {
    data.frame(
      model = "loc/locsca/lpm", backend = entry$name,
      implementation = entry$description,
      objective_preserving = entry$objective_preserving,
      optional_package = NA_character_, stringsAsFactors = FALSE
    )
  })
  custom_dr <- lapply(.cf_extension_registry$distribution, function(entry) {
    data.frame(
      model = "logit/probit/cloglog", backend = entry$name,
      implementation = entry$description,
      objective_preserving = entry$objective_preserving,
      optional_package = NA_character_, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(builtins), custom_linear, custom_dr))
}

#' Describe the available quantile-regression solvers
#'
#' @return A data frame describing each solver, its implementation, and whether
#'   it exactly minimizes the quantile-regression objective.
#' @export
qr_solver_registry <- function() {
  builtins <- data.frame(
    solver = builtin_qr_solvers(),
    implementation = c(
      "quantreg::rq.fit.br",
      "quantreg::rq.fit.fnb",
      "quantreg::rq.fit.pfn",
      "quantreg::rq.fit.qfnb",
      "quantreg::rq.fit.pfnb",
      "quantreg::rq.fit.ppro(pmethod='br')",
      "quantreg::rq.fit.ppro(pmethod='fn')",
      "qrprocess 1.1.3-compatible one-step R implementation",
      "CuPy batched ADMM (experimental)",
      "sample-size and quantile-grid dispatcher"
    ),
    exact = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, NA),
    process_aware = c(
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, NA
    ),
    stringsAsFactors = FALSE
  )
  custom <- lapply(.cf_extension_registry$qr, function(entry) {
    data.frame(
      solver = entry$name,
      implementation = entry$description,
      exact = entry$exact,
      process_aware = entry$process_aware,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(builtins), custom))
}

resolve_qr_solver <- function(solver, n, taus) {
  solver <- match.arg(solver, supported_qr_solvers())
  if (solver != "auto") return(solver)
  n <- as.integer(n)
  ntau <- length(taus)
  if (ntau < 10L) {
    return(if (n < 10000L) "br" else "pfn")
  }
  if (n < 10000L) return("qfnb")
  "pfnb"
}

assert_scalar_integer <- function(x, name, minimum = 0L) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= minimum && x <= .Machine$integer.max && x == floor(x)
  if (!valid) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  as.integer(x)
}

assert_probability <- function(x, name, open = TRUE) {
  valid <- length(x) == 1L && is.finite(x)
  if (open) valid <- valid && x > 0 && x < 1
  else valid <- valid && x >= 0 && x <= 1
  if (!valid) stop(name, " must be between 0 and 1", call. = FALSE)
  as.numeric(x)
}

assert_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

normalize_weights <- function(weights) {
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must be finite and strictly positive", call. = FALSE)
  }
  as.numeric(weights) / mean(as.numeric(weights))
}

weighted_quantile <- function(x, weights, probs, legacy = FALSE) {
  if (length(x) != length(weights)) {
    stop("x and weights have incompatible sizes", call. = FALSE)
  }
  if (isTRUE(legacy)) {
    return(as.numeric(Hmisc::wtd.quantile(
      x,
      weights = weights,
      probs = probs,
      na.rm = TRUE,
      normwt = FALSE
    )))
  }
  weights <- normalize_weights(weights)
  as.numeric(Hmisc::wtd.quantile(
    x,
    weights = weights,
    probs = probs,
    na.rm = TRUE,
    normwt = TRUE
  ))
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (file.exists(path) && unlink(path) != 0L) {
    stop("Could not replace checkpoint: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not create checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

object_md5 <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

measure_resources <- function(fun) {
  stopifnot(is.function(fun))
  invisible(gc(reset = TRUE))
  started <- proc.time()[["elapsed"]]
  value <- fun()
  elapsed <- proc.time()[["elapsed"]] - started
  memory <- gc()
  list(
    value = value,
    elapsed_seconds = unname(elapsed),
    peak_r_heap_mb = unname(sum(memory[, 6L, drop = TRUE]))
  )
}

package_worker_init <- function(cluster) {
  paths <- .libPaths()
  extensions <- registry_snapshot()
  parallel::clusterCall(cluster, function(lib_paths, registry) {
    .libPaths(lib_paths)
    suppressPackageStartupMessages(library(scalableCounterfactual))
    restore <- get(
      "restore_registry_snapshot",
      envir = asNamespace("scalableCounterfactual")
    )
    restore(registry)
    NULL
  }, paths, extensions)
  invisible(cluster)
}

safe_max_abs_difference <- function(x, y) {
  if (is.null(x) || is.null(y) || !identical(dim(x), dim(y))) return(NA_real_)
  max(abs(x - y), na.rm = TRUE)
}
