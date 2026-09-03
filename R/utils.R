supported_cf_models <- function() {
  c("qr", "cqr", "loc", "locsca", "cox", "logit", "probit", "cloglog", "lpm")
}

is_quantile_process_model <- function(model) model %in% c("qr", "cqr")

is_distribution_regression_model <- function(model) {
  model %in% c("logit", "probit", "cloglog", "lpm")
}

.datatable.aware <- TRUE
.cf_output_schema_version <- "1.1"

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

weighted_rank_tolerance <- function(total_weight) {
  64 * .Machine$double.eps * max(1, abs(as.numeric(total_weight)))
}

weighted_rank_index <- function(cumulative, target, total_weight) {
  tolerance <- weighted_rank_tolerance(total_weight)
  hit <- which(cumulative >= target - tolerance)[1L]
  if (is.na(hit)) length(cumulative) else hit
}

weighted_type7_quantile <- function(
    x, weights, probs, normalization_n = length(x)) {
  if (length(x) != length(weights)) {
    stop("x and weights have incompatible sizes", call. = FALSE)
  }
  if (!length(x) || any(!is.finite(x))) {
    stop("x must contain finite values", call. = FALSE)
  }
  probs <- as.numeric(probs)
  if (!length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop("probs must lie between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(normalization_n) || length(normalization_n) != 1L ||
      !is.finite(normalization_n) || normalization_n < 1) {
    stop("normalization_n must be one finite value >= 1", call. = FALSE)
  }
  weights <- normalize_weights(weights) *
    (as.numeric(normalization_n) / length(weights))
  ordering <- order(x)
  ordered_x <- as.numeric(x[ordering])
  cumulative <- cumsum(weights[ordering])
  total_positions <- as.numeric(normalization_n)
  cumulative[[length(cumulative)]] <- total_positions
  positions <- 1 + (total_positions - 1) * probs
  low <- pmax(floor(positions), 1)
  high <- pmin(low + 1, total_positions)
  target_values <- vapply(c(low, high), function(target) {
    ordered_x[[weighted_rank_index(cumulative, target, total_positions)]]
  }, numeric(1L))
  interpolation <- positions %% 1
  k <- length(probs)
  (1 - interpolation) * target_values[seq_len(k)] +
    interpolation * target_values[k + seq_len(k)]
}

weighted_quantile <- function(
    x, weights, probs, legacy = FALSE, normalization_n = length(x)) {
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
  weighted_type7_quantile(x, weights, probs, normalization_n)
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  if (dir.exists(path)) {
    stop("checkpoint path is an existing directory: ", path, call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  backup <- paste0(path, ".bak_", Sys.getpid())
  committed <- FALSE
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  on.exit({
    if (committed && file.exists(backup)) unlink(backup)
  }, add = TRUE)
  saveRDS(object, temporary, compress = compress)
  had_existing <- file.exists(path)
  if (had_existing && !file.rename(path, backup)) {
    stop("Could not stage the existing checkpoint: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    restored <- !had_existing || !file.exists(backup) ||
      file.rename(backup, path)
    if (!restored) {
      stop(
        "Could not create checkpoint or restore its predecessor; preserved ",
        "the prior checkpoint at: ", backup,
        call. = FALSE
      )
    }
    stop("Could not create checkpoint: ", path, call. = FALSE)
  }
  committed <- TRUE
  if (had_existing && file.exists(backup)) unlink(backup)
  invisible(path)
}

atomic_write_output_files <- function(
    output_dir, managed_files, writer, required_files = character()) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty path", call. = FALSE)
  }
  if (!is.function(writer)) stop("writer must be a function", call. = FALSE)
  managed_files <- unique(as.character(managed_files))
  required_files <- unique(as.character(required_files))
  invalid <- function(x) {
    !nzchar(x) | x %in% c(".", "..") | dirname(x) != "." |
      grepl("[/\\\\]", x)
  }
  if (!length(managed_files) || any(invalid(managed_files)) ||
      any(invalid(required_files)) ||
      !all(required_files %in% managed_files)) {
    stop("managed and required files must be safe base names", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("could not create output_dir: ", output_dir, call. = FALSE)
  }
  managed_paths <- file.path(output_dir, managed_files)
  managed_directories <- managed_files[dir.exists(managed_paths)]
  if (length(managed_directories)) {
    stop(
      "managed output path is an existing directory: ",
      paste(managed_directories, collapse = ", "),
      call. = FALSE
    )
  }
  stage <- tempfile(".scalablecf-stage-", tmpdir = output_dir)
  backup <- tempfile(".scalablecf-backup-", tmpdir = output_dir)
  dir.create(stage)
  dir.create(backup)
  committed <- FALSE
  on.exit({
    if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
    backup_empty <- dir.exists(backup) &&
      !length(list.files(backup, all.files = TRUE, no.. = TRUE))
    if ((committed || backup_empty) && dir.exists(backup)) {
      unlink(backup, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  writer(stage)
  missing_required <- required_files[!file.exists(file.path(stage, required_files))]
  if (length(missing_required)) {
    stop(
      "staged output is missing required files: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
  staged_files <- managed_files[file.exists(file.path(stage, managed_files))]
  if (any(dir.exists(file.path(stage, staged_files)))) {
    stop("writer produced a directory where a file was required", call. = FALSE)
  }
  existing_files <- managed_files[file.exists(file.path(output_dir, managed_files))]
  moved_old <- character()
  moved_new <- character()
  rollback <- function() {
    if (length(moved_new)) {
      unlink(file.path(output_dir, moved_new), force = TRUE)
    }
    for (name in rev(moved_old)) {
      file.rename(file.path(backup, name), file.path(output_dir, name))
    }
  }
  for (name in existing_files) {
    if (!file.rename(file.path(output_dir, name), file.path(backup, name))) {
      rollback()
      stop("could not stage existing output file: ", name, call. = FALSE)
    }
    moved_old <- c(moved_old, name)
  }
  for (name in staged_files) {
    if (!file.rename(file.path(stage, name), file.path(output_dir, name))) {
      rollback()
      stop("could not commit staged output file: ", name, call. = FALSE)
    }
    moved_new <- c(moved_new, name)
  }
  committed <- TRUE
  if (dir.exists(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  invisible(output_dir)
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
