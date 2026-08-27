resolve_linear_backend <- function(backend) {
  backend <- match.arg(backend, supported_linear_backends())
  if (backend == "auto") "qr" else backend
}

require_optional_backend <- function(package, backend) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "backend '", backend, "' requires the optional package '", package,
      "'",
      call. = FALSE
    )
  }
}

make_weighted_linear_solver <- function(X, weights, backend = "auto") {
  X <- as.matrix(X)
  weights <- normalize_weights(weights)
  resolved <- resolve_linear_backend(backend)
  sqrt_weights <- sqrt(weights)
  weighted_X <- X * sqrt_weights
  warnings <- character()

  if (resolved == "qr") {
    decomposition <- qr(weighted_X, LAPACK = FALSE)
    if (decomposition$rank < ncol(X)) {
      stop("weighted QR design is rank deficient", call. = FALSE)
    }
    solve_response <- function(y) {
      coefficients <- as.numeric(qr.coef(decomposition, y * sqrt_weights))
      list(
        coefficients = coefficients,
        residuals = as.numeric(y - X %*% coefficients)
      )
    }
  } else if (resolved == "chol") {
    crossproduct <- crossprod(weighted_X)
    factor <- tryCatch(chol(crossproduct), error = function(error) NULL)
    if (is.null(factor)) {
      stop("weighted Cholesky design is not positive definite", call. = FALSE)
    }
    solve_response <- function(y) {
      right_hand_side <- crossprod(weighted_X, y * sqrt_weights)
      coefficients <- as.numeric(backsolve(
        factor,
        forwardsolve(t(factor), right_hand_side)
      ))
      list(
        coefficients = coefficients,
        residuals = as.numeric(y - X %*% coefficients)
      )
    }
  } else if (resolved == "fastglm") {
    require_optional_backend("fastglm", "fastglm")
    solve_response <- function(y) {
      fit_warnings <- character()
      fit <- withCallingHandlers(
        fastglm::fastglmPure(
          X,
          as.numeric(y),
          family = stats::gaussian(),
          weights = weights,
          maxit = 100L,
          tol = 1e-10
        ),
        warning = function(condition) {
          fit_warnings <<- c(fit_warnings, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }
      )
      if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
        stop("fastglm weighted linear fit did not converge", call. = FALSE)
      }
      warnings <<- unique(c(warnings, fit_warnings))
      list(
        coefficients = as.numeric(fit$coefficients),
        residuals = as.numeric(y - X %*% fit$coefficients)
      )
    }
  } else {
    entry <- custom_registry_entry("linear", resolved)
    if (is.null(entry)) stop("unregistered linear backend: ", resolved,
                             call. = FALSE)
    solve_response <- function(y) {
      backend_warnings <- character()
      returned <- withCallingHandlers(
        entry$fit(X, as.numeric(y), weights),
        warning = function(condition) {
          backend_warnings <<- c(
            backend_warnings, conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
      if (is.numeric(returned)) returned <- list(coefficients = returned)
      if (!is.list(returned) || is.null(returned$coefficients)) {
        stop("custom linear backend must return coefficients", call. = FALSE)
      }
      coefficients <- as.numeric(returned$coefficients)
      if (length(coefficients) != ncol(X) || any(!is.finite(coefficients))) {
        stop("custom linear backend returned invalid coefficients",
             call. = FALSE)
      }
      warnings <<- unique(c(
        warnings, backend_warnings, as.character(returned$warnings %||% character())
      ))
      residuals <- returned$residuals %||%
        as.numeric(y - X %*% coefficients)
      if (length(residuals) != nrow(X) || any(!is.finite(residuals))) {
        stop("custom linear backend returned invalid residuals", call. = FALSE)
      }
      list(coefficients = coefficients, residuals = as.numeric(residuals))
    }
  }

  list(
    fit = solve_response,
    backend_requested = backend,
    backend = resolved,
    warnings = function() warnings
  )
}

fit_location_model <- function(
    X, y, weights, taus, linear_backend = "auto") {
  solver <- make_weighted_linear_solver(X, weights, linear_backend)
  fitted <- solver$fit(y)
  residual_quantiles <- weighted_quantile(fitted$residuals, weights, taus)
  structure(list(
    model = "loc",
    coefficients = fitted$coefficients,
    residual_quantiles = residual_quantiles,
    taus = taus,
    backend_requested = solver$backend_requested,
    backend = solver$backend,
    warnings = solver$warnings()
  ), class = c("cf_loc_fit", "cf_conditional_fit"))
}

fit_location_scale_model <- function(
    X, y, weights, taus, linear_backend = "auto") {
  solver <- make_weighted_linear_solver(X, weights, linear_backend)
  location_fit <- solver$fit(y)
  residuals <- location_fit$residuals
  floor_value <- .Machine$double.eps
  scale_outcome <- log(pmax(residuals^2, floor_value))
  scale_fit <- solver$fit(scale_outcome)
  fitted_scale <- sqrt(exp(drop(X %*% scale_fit$coefficients)))
  standardized <- residuals / pmax(fitted_scale, sqrt(floor_value))
  residual_quantiles <- weighted_quantile(standardized, weights, taus)
  structure(list(
    model = "locsca",
    location_coefficients = location_fit$coefficients,
    scale_coefficients = scale_fit$coefficients,
    residual_quantiles = residual_quantiles,
    taus = taus,
    backend_requested = solver$backend_requested,
    backend = solver$backend,
    warnings = solver$warnings()
  ), class = c("cf_locsca_fit", "cf_conditional_fit"))
}
