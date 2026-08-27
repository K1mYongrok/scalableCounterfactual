.cf_extension_registry <- new.env(parent = emptyenv())
.cf_extension_registry$qr <- list()
.cf_extension_registry$linear <- list()
.cf_extension_registry$distribution <- list()

validate_extension_name <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !grepl("^[A-Za-z][A-Za-z0-9._]*$", name)) {
    stop(
      "extension name must start with a letter and contain only letters, ",
      "numbers, '.', or '_'",
      call. = FALSE
    )
  }
  tolower(name)
}

validate_extension_dependencies <- function(dependencies) {
  if (is.null(dependencies) || (is.list(dependencies) && !length(dependencies))) {
    return(list())
  }
  if (!is.list(dependencies) || is.null(names(dependencies)) ||
      any(!nzchar(names(dependencies))) || anyDuplicated(names(dependencies))) {
    stop("dependencies must be a uniquely named list", call. = FALSE)
  }
  dependencies
}

capture_extension_function <- function(fun, dependencies = NULL) {
  dependencies <- validate_extension_dependencies(dependencies)
  captured <- new.env(parent = baseenv())
  visited <- new.env(parent = emptyenv())

  capture_value <- function(name, value) {
    if (is.function(value) && !is.primitive(value) &&
        !isNamespace(environment(value))) {
      copied <- value
      environment(copied) <- captured
      assign(name, copied, envir = captured)
      capture_globals(value)
    } else {
      assign(name, value, envir = captured)
    }
  }
  capture_globals <- function(candidate) {
    key <- object_md5(candidate)
    if (exists(key, envir = visited, inherits = FALSE)) return(invisible(NULL))
    assign(key, TRUE, envir = visited)
    globals <- codetools::findGlobals(candidate, merge = FALSE)
    symbols <- unique(c(globals$variables, globals$functions))
    symbols <- setdiff(symbols, c(names(formals(candidate)), "..."))
    for (name in symbols) {
      if (exists(name, envir = captured, inherits = FALSE)) next
      if (name %in% names(dependencies)) {
        capture_value(name, dependencies[[name]])
      } else if (!exists(name, envir = baseenv(), inherits = FALSE) &&
                 exists(name, envir = environment(candidate), inherits = TRUE)) {
        capture_value(
          name,
          get(name, envir = environment(candidate), inherits = TRUE)
        )
      }
    }
    invisible(NULL)
  }

  for (name in names(dependencies)) {
    capture_value(name, dependencies[[name]])
  }
  capture_globals(fun)
  captured_fun <- fun
  environment(captured_fun) <- captured
  captured_fun
}

function_fingerprint <- function(fun, version = NULL) {
  object_md5(list(function_object = fun, version = version))
}

registry_snapshot <- function() {
  list(
    qr = .cf_extension_registry$qr,
    linear = .cf_extension_registry$linear,
    distribution = .cf_extension_registry$distribution
  )
}

restore_registry_snapshot <- function(snapshot) {
  if (is.null(snapshot)) return(invisible(NULL))
  for (type in c("qr", "linear", "distribution")) {
    .cf_extension_registry[[type]] <- snapshot[[type]] %||% list()
  }
  invisible(NULL)
}

extension_registry_fingerprint <- function(selection = NULL) {
  snapshot <- registry_snapshot()
  if (!is.null(selection)) {
    for (type in names(snapshot)) {
      requested <- selection[[type]] %||% character()
      snapshot[[type]] <- snapshot[[type]][intersect(
        names(snapshot[[type]]), requested
      )]
    }
  }
  compact <- lapply(snapshot, function(entries) {
    lapply(entries, function(entry) {
      list(
        name = entry$name,
        type = entry$type %||% "qr",
        exact = entry$exact %||% NULL,
        process_aware = entry$process_aware %||% NULL,
        objective_preserving = entry$objective_preserving %||% NULL,
        version = entry$version,
        function_fingerprint = entry$function_fingerprint
      )
    })
  })
  object_md5(compact)
}

active_extension_fingerprint <- function(model, solver, control, point = NULL) {
  selection <- list(qr = character(), linear = character(), distribution = character())
  if (is_quantile_process_model(model)) {
    selection$qr <- if (!is.null(point) && length(point$fits)) {
      unique(vapply(point$fits, `[[`, character(1L), "solver"))
    } else {
      solver
    }
  } else if (model %in% c("loc", "locsca", "lpm")) {
    selection$linear <- if (!is.null(point) && length(point$fits)) {
      unique(vapply(point$fits, `[[`, character(1L), "backend"))
    } else {
      control$linear_backend
    }
  } else {
    selection$distribution <- if (!is.null(point) && length(point$fits)) {
      unique(vapply(point$fits, `[[`, character(1L), "backend"))
    } else {
      control$dr_backend
    }
  }
  extension_registry_fingerprint(selection)
}

custom_registry_entry <- function(type, name) {
  .cf_extension_registry[[type]][[name]]
}

#' Register a custom quantile-regression solver
#'
#' The supplied function is called as `fit(X, y, weights, taus, control)` and
#' must return either a coefficient matrix or a list containing
#' `coefficients`. Coefficients must have one row per design column and one
#' column per quantile and must be in the original design units. Optional list
#' members are `flag` (zero denotes convergence), `iterations`, and `warnings`.
#' The function is responsible for applying case weights and minimizing the
#' claimed QR objective.
#'
#' @param name Unique solver name.
#' @param fit Solver function following the contract above.
#' @param exact Whether the function exactly minimizes the linear QR objective.
#' @param process_aware Whether it estimates multiple quantiles jointly.
#' @param description Human-readable implementation description.
#' @param overwrite Replace an existing custom solver with the same name.
#' @param dependencies Optional named list of non-base objects used by `fit`.
#'   Direct lexical dependencies are captured automatically; this argument is
#'   useful for dependencies reached through dynamic lookup.
#' @param version Optional user-defined implementation version included in
#'   checkpoint identities.
#' @return Invisibly returns the normalized solver name.
#' @export
register_qr_solver <- function(
    name, fit, exact, process_aware = FALSE, description = NULL,
    overwrite = FALSE, dependencies = NULL, version = NULL) {
  name <- validate_extension_name(name)
  if (name %in% builtin_qr_solvers()) {
    stop("built-in QR solvers cannot be replaced", call. = FALSE)
  }
  if (!is.function(fit)) stop("fit must be a function", call. = FALSE)
  dependencies <- validate_extension_dependencies(dependencies)
  if (!is.null(version) && (!is.character(version) || length(version) != 1L ||
                           is.na(version) || !nzchar(version))) {
    stop("version must be NULL or one nonempty character string", call. = FALSE)
  }
  exact <- assert_scalar_logical(exact, "exact")
  process_aware <- assert_scalar_logical(process_aware, "process_aware")
  overwrite <- assert_scalar_logical(overwrite, "overwrite")
  if (!is.null(custom_registry_entry("qr", name)) && !overwrite) {
    stop("custom QR solver already registered: ", name, call. = FALSE)
  }
  if (is.null(description)) description <- paste("custom QR solver", name)
  if (!is.character(description) || length(description) != 1L) {
    stop("description must be one character string", call. = FALSE)
  }
  fit <- capture_extension_function(fit, dependencies)
  .cf_extension_registry$qr[[name]] <- list(
    name = name,
    fit = fit,
    exact = exact,
    process_aware = process_aware,
    description = description,
    version = version,
    function_fingerprint = function_fingerprint(fit, version)
  )
  invisible(name)
}

#' @rdname register_qr_solver
#' @param quiet Do not error when the requested custom extension is absent.
#' @export
unregister_qr_solver <- function(name, quiet = FALSE) {
  name <- validate_extension_name(name)
  quiet <- assert_scalar_logical(quiet, "quiet")
  if (is.null(custom_registry_entry("qr", name))) {
    if (!quiet) stop("custom QR solver is not registered: ", name, call. = FALSE)
    return(invisible(FALSE))
  }
  .cf_extension_registry$qr[[name]] <- NULL
  invisible(TRUE)
}

#' Register a custom conditional-model computation backend
#'
#' Linear backends support `loc`, `locsca`, and `lpm`. Their function is called
#' as `fit(X, y, weights)` and returns coefficients or a list containing
#' `coefficients`, with optional `residuals` and `warnings`. Distribution
#' backends support logit/probit/cloglog threshold regressions. Their function is called
#' as `fit(X, response, weights, model, start, maxit, tolerance)` and returns
#' coefficients or a list with `coefficients` and optional `converged`,
#' `boundary`, `iterations`, and `warnings`.
#'
#' @param name Unique backend name.
#' @param type Either `linear` or `distribution`.
#' @param fit Backend function following the corresponding contract.
#' @param objective_preserving Whether the backend preserves the package's
#'   unpenalized conditional-model objective.
#' @param description Human-readable implementation description.
#' @param overwrite Replace an existing custom backend with the same name/type.
#' @param dependencies Optional named list of non-base objects used by `fit`.
#' @param version Optional implementation version included in checkpoint
#'   identities.
#' @return Invisibly returns the normalized backend name.
#' @export
register_conditional_backend <- function(
    name, type = c("linear", "distribution"), fit,
    objective_preserving = TRUE, description = NULL, overwrite = FALSE,
    dependencies = NULL, version = NULL) {
  name <- validate_extension_name(name)
  type <- match.arg(type)
  builtins <- if (type == "linear") {
    builtin_linear_backends()
  } else {
    builtin_dr_backends()
  }
  if (name %in% builtins) {
    stop("built-in conditional backends cannot be replaced", call. = FALSE)
  }
  if (!is.function(fit)) stop("fit must be a function", call. = FALSE)
  dependencies <- validate_extension_dependencies(dependencies)
  if (!is.null(version) && (!is.character(version) || length(version) != 1L ||
                           is.na(version) || !nzchar(version))) {
    stop("version must be NULL or one nonempty character string", call. = FALSE)
  }
  objective_preserving <- assert_scalar_logical(
    objective_preserving, "objective_preserving"
  )
  overwrite <- assert_scalar_logical(overwrite, "overwrite")
  if (!is.null(custom_registry_entry(type, name)) && !overwrite) {
    stop("custom ", type, " backend already registered: ", name,
         call. = FALSE)
  }
  if (is.null(description)) description <- paste("custom", type, "backend", name)
  if (!is.character(description) || length(description) != 1L) {
    stop("description must be one character string", call. = FALSE)
  }
  fit <- capture_extension_function(fit, dependencies)
  .cf_extension_registry[[type]][[name]] <- list(
    name = name,
    type = type,
    fit = fit,
    objective_preserving = objective_preserving,
    description = description,
    version = version,
    function_fingerprint = function_fingerprint(fit, version)
  )
  invisible(name)
}

#' @rdname register_conditional_backend
#' @param quiet Do not error when the requested custom extension is absent.
#' @export
unregister_conditional_backend <- function(
    name, type = c("linear", "distribution"), quiet = FALSE) {
  name <- validate_extension_name(name)
  type <- match.arg(type)
  quiet <- assert_scalar_logical(quiet, "quiet")
  if (is.null(custom_registry_entry(type, name))) {
    if (!quiet) {
      stop("custom ", type, " backend is not registered: ", name,
           call. = FALSE)
    }
    return(invisible(FALSE))
  }
  .cf_extension_registry[[type]][[name]] <- NULL
  invisible(TRUE)
}
