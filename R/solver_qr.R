# Exact QR fits are delegated to public quantreg routines. The one-step process
# is the sole translated solver; its frozen source mapping is documented in
# inst/provenance/qrprocess_onestep.yml.
normalize_coefficient_matrix <- function(coefficients, p, taus) {
  coefficients <- as.matrix(coefficients)
  if (nrow(coefficients) != p && ncol(coefficients) == p) {
    coefficients <- t(coefficients)
  }
  if (!identical(dim(coefficients), c(p, length(taus)))) {
    stop(
      "QR solver returned coefficients with dimension ",
      paste(dim(coefficients), collapse = " x "),
      "; expected ", p, " x ", length(taus),
      call. = FALSE
    )
  }
  coefficients
}

single_tau_qr <- function(X, y, tau, solver) {
  switch(
    solver,
    br = quantreg::rq.fit.br(X, y, tau = tau),
    fn = quantreg::rq.fit.fnb(X, y, tau = tau),
    pfn = quantreg::rq.fit.pfn(X, y, tau = tau),
    stop("single-tau solver not implemented: ", solver, call. = FALSE)
  )
}

precondition_qr_design <- function(X, enabled, tolerance = 1e-10) {
  p <- ncol(X)
  identity <- diag(p)
  if (!isTRUE(enabled)) {
    return(list(
      X = X,
      transform = identity,
      method = "none",
      gram_rcond = NA_real_,
      condition_estimate = NA_real_,
      fallback_reason = NA_character_
    ))
  }

  gram <- crossprod(X)
  gram_rcond <- suppressWarnings(rcond(gram))
  condition_estimate <- if (is.finite(gram_rcond) && gram_rcond > 0) {
    1 / sqrt(gram_rcond)
  } else {
    Inf
  }
  cholesky_threshold <- .Machine$double.eps^(2 / 3)
  gram_factor <- if (is.finite(gram_rcond) &&
                     gram_rcond >= cholesky_threshold) {
    tryCatch(chol(gram), error = function(e) NULL)
  } else {
    NULL
  }
  if (!is.null(gram_factor)) {
    transform <- backsolve(gram_factor, identity)
    return(list(
      X = X %*% transform,
      transform = transform,
      method = "cholesky",
      gram_rcond = gram_rcond,
      condition_estimate = condition_estimate,
      fallback_reason = NA_character_
    ))
  }

  qr_fit <- qr(X, tol = tolerance, LAPACK = FALSE)
  if (qr_fit$rank < p) {
    stop("QR preconditioning found a rank-deficient weighted design",
         call. = FALSE)
  }
  upper <- qr.R(qr_fit, complete = FALSE)
  inverse_upper <- backsolve(upper, identity)
  transform <- matrix(0, p, p)
  transform[qr_fit$pivot, ] <- inverse_upper
  list(
    X = X %*% transform,
    transform = transform,
    method = "pivoted_qr",
    gram_rcond = gram_rcond,
    condition_estimate = condition_estimate,
    fallback_reason = if (is.finite(gram_rcond)) {
      "ill_conditioned_gram_matrix"
    } else {
      "nonfinite_gram_condition"
    }
  )
}

stata_weighted_quantile_type2 <- function(x, weights, probs) {
  x <- as.numeric(x)
  weights <- normalize_weights(weights)
  probs <- as.numeric(probs)
  ordering <- order(x, weights)
  x <- x[ordering]
  weights <- weights[ordering]
  runs <- !duplicated(x, fromLast = TRUE)
  unique_x <- x[runs]
  cumulative <- cumsum(weights)[runs]
  total <- cumulative[[length(cumulative)]]
  vapply(probs, function(probability) {
    target <- probability * total
    position <- which(cumulative >= target)[[1L]]
    if (cumulative[[position]] == target && position < length(unique_x)) {
      (unique_x[[position]] + unique_x[[position + 1L]]) / 2
    } else {
      unique_x[[position]]
    }
  }, numeric(1L))
}

stata_weighted_scale <- function(residuals, weights) {
  weights <- normalize_weights(weights)
  n <- length(residuals)
  center <- sum(weights * residuals) / sum(weights)
  variance <- if (n > 1L) {
    sum(weights * (residuals - center)^2) / (n - 1L)
  } else {
    0
  }
  quartiles <- stata_weighted_quantile_type2(
    residuals, weights, c(0.25, 0.75)
  )
  min(sqrt(max(variance, 0)), diff(quartiles) / 1.34)
}

stata_onestep_bandwidth <- function(
    tau, n, residuals, weights,
    method = c("hall_sheather", "bofinger")) {
  method <- match.arg(method)
  x0 <- stats::qnorm(tau)
  f0 <- stats::dnorm(x0)
  probability_bandwidth <- if (method == "hall_sheather") {
    n^(-1 / 3) * stats::qnorm(0.975)^(2 / 3) *
      ((1.5 * f0^2) / (2 * x0^2 + 1))^(1 / 3)
  } else {
    n^(-0.2) * ((4.5 * f0^4) / (2 * x0^2 + 1)^2)^0.2
  }
  if (tau - probability_bandwidth < 0.001) {
    probability_bandwidth <- tau * 0.5
  }
  if (tau + probability_bandwidth > 0.999) {
    probability_bandwidth <- (1 - tau) * 0.5
  }
  bandwidth <- (
    stats::qnorm(tau + probability_bandwidth) -
      stats::qnorm(tau - probability_bandwidth)
  ) * stata_weighted_scale(residuals, weights)
  if (!is.finite(bandwidth) || bandwidth <= 0) NA_real_ else bandwidth
}

exact_initial_qr <- function(X, y, weights, tau, solver) {
  fit <- single_tau_qr(X * weights, y * weights, tau, solver)
  list(
    coefficients = as.numeric(fit$coefficients),
    iterations = fit$nit
  )
}

# Native R transcription of qrprocess.ado 1.1.3, rq_1step(), lines 2271-2359
# at commit ec56830ef9c84ce54411ad59c5ce94535847d9df. The ordering,
# bandwidth, weighted scale, Jacobian, score, and fallback sequence mirror the
# Mata routine. Lower-level exact QR fits use quantreg, so cross-language
# bitwise identity is not claimed.
fit_qr_onestep <- function(
    X, y, weights, taus,
    first_solver = c("auto", "br", "fn", "pfn"),
    bandwidth_method = c("hall_sheather", "bofinger")) {
  first_solver <- match.arg(first_solver)
  bandwidth_method <- match.arg(bandwidth_method)
  if (length(taus) > 1L) {
    largest_step <- max(diff(taus))
    if (largest_step > 0.0501) {
      stop(
        "onestep requires a conditional-quantile grid with maximum spacing <= 0.05",
        call. = FALSE
      )
    }
  }
  n <- nrow(X)
  p <- ncol(X)
  ntau <- length(taus)
  resolved_first <- if (first_solver == "auto") {
    if (n < 100000L) "br" else "pfn"
  } else {
    first_solver
  }
  coefficients <- matrix(NA_real_, p, ntau)
  inverse_jacobians <- vector("list", ntau)
  convergence <- rep(NA_integer_, ntau)
  start <- which.min(abs(taus - 0.5))
  traversal <- c(
    seq.int(start, ntau),
    if (start > 1L) seq.int(start - 1L, 1L) else integer()
  )
  directions <- c(
    rep(1L, ntau - start + 1L),
    rep(-1L, start - 1L)
  )
  fallback_taus <- numeric()
  for (position in seq_along(traversal)) {
    current <- traversal[[position]]
    tau <- taus[[current]]
    candidate <- NULL
    if (position > 1L) {
      previous <- current - directions[[position]]
      previous_beta <- coefficients[, previous]
      fitted <- drop(X %*% previous_beta)
      residuals <- y - fitted
      bandwidth <- stata_onestep_bandwidth(
        tau, n, residuals, weights, method = bandwidth_method
      )
      if (is.finite(bandwidth)) {
        density <- stats::dnorm(residuals / bandwidth) / bandwidth
        jacobian <- crossprod(
          X,
          X * as.numeric(density * weights)
        ) / n
        inverse_jacobian <- tryCatch(
          solve(jacobian),
          error = function(e) NULL
        )
        if (!is.null(inverse_jacobian)) {
          score <- colSums(
            X * as.numeric(weights * (tau - (y <= fitted)))
          ) / sum(weights)
          update <- previous_beta + drop(inverse_jacobian %*% score)
          if (all(is.finite(update))) {
            candidate <- update
            inverse_jacobians[[current]] <- inverse_jacobian
            convergence[[current]] <- 1L
          }
        }
      }
    }
    if (is.null(candidate)) {
      exact <- exact_initial_qr(
        X, y, weights, tau, resolved_first
      )
      candidate <- exact$coefficients
      convergence[[current]] <- 0L
      if (position > 1L) fallback_taus <- c(fallback_taus, tau)
    }
    coefficients[, current] <- candidate
  }
  list(
    coefficients = coefficients,
    flag = rep(0L, ntau),
    nit = convergence,
    inverse_jacobians = inverse_jacobians,
    traversal = traversal,
    directions = directions,
    first_solver_requested = first_solver,
    first_solver_resolved = resolved_first,
    bandwidth_method = bandwidth_method,
    fallback_taus = fallback_taus,
    warnings = c(
      if (length(taus) > 1L && max(diff(taus)) > 0.0101) {
        "onestep mirrors qrprocess's warning for quantile spacing above 0.01"
      },
      if (length(fallback_taus)) {
        paste0(
          "onestep used exact ", resolved_first,
          " fallback at tau=", paste(signif(fallback_taus, 5), collapse = ",")
        )
      }
    )
  )
}

#' Fit a weighted conditional quantile process
#'
#' This is the solver-neutral QR backend. Positive case weights are applied by
#' premultiplying both the design rows and outcomes before calling lower-level
#' `quantreg` routines. This also preserves weights for multi-tau PFNB/QFNB.
#'
#' @param X Numeric design matrix including an intercept.
#' @param y Numeric outcome.
#' @param weights Strictly positive case weights.
#' @param taus Sorted conditional quantiles.
#' @param solver One of `br`, `fn`, `pfn`, `qfnb`, `pfnb`, `proqreg`,
#'   `profn`, `onestep`, or `auto`.
#' @param precondition Apply an invertible design preconditioner before fitting
#'   and transform coefficients back to the original units.
#' @param onestep_first_solver Exact solver used to initialize the one-step
#'   process. `auto` uses `br` below 100,000 observations and `pfn` otherwise.
#' @param onestep_bandwidth Bandwidth rule for the one-step process.
#' @param gpu_control Internal `cf_control` object required by the experimental
#'   `cuda_admm` solver.
#' @return A standardized QR fit list.
#' @export
fit_weighted_qr <- function(
    X,
    y,
    weights,
    taus,
    solver = "auto",
    precondition = TRUE,
    onestep_first_solver = c("auto", "br", "fn", "pfn"),
    onestep_bandwidth = c("hall_sheather", "bofinger"),
    gpu_control = NULL) {
  solver <- match.arg(solver, supported_qr_solvers())
  solver_requested <- solver
  onestep_first_solver <- match.arg(onestep_first_solver)
  onestep_bandwidth <- match.arg(onestep_bandwidth)
  X <- as.matrix(X)
  y <- as.numeric(y)
  weights <- normalize_weights(weights)
  taus <- sort(unique(as.numeric(taus)))
  if (nrow(X) != length(y) || length(y) != length(weights)) {
    stop("X, y, and weights have incompatible sizes", call. = FALSE)
  }
  if (!length(taus) || any(taus <= 0 | taus >= 1)) {
    stop("taus must lie strictly between 0 and 1", call. = FALSE)
  }
  solver <- resolve_qr_solver(solver, nrow(X), taus)
  custom_solver <- custom_registry_entry("qr", solver)
  weighted_X <- X * weights
  weighted_y <- y * weights
  precondition_effective <- isTRUE(precondition) &&
    !solver %in% c("onestep", "cuda_admm") &&
    is.null(custom_solver)
  preconditioned <- precondition_qr_design(weighted_X, precondition_effective)
  preconditioner <- preconditioned$transform
  weighted_X <- preconditioned$X
  conditioned_X <- X %*% preconditioner
  fit <- NULL
  single_fits <- NULL
  solver_warnings <- character()
  process_fallback_taus <- numeric()
  process_fallback_solver <- NA_character_
  process_fallback_reason <- NA_character_
  capture_solver_warnings <- function(expression) {
    withCallingHandlers(
      expression,
      warning = function(warning_condition) {
        solver_warnings <<- c(
          solver_warnings,
          conditionMessage(warning_condition)
        )
        invokeRestart("muffleWarning")
      }
    )
  }
  pfnb_m0 <- min(
    nrow(weighted_X),
    max(
      ncol(weighted_X) + 1L,
      floor(nrow(weighted_X)^(2 / 3) * sqrt(ncol(weighted_X)))
    )
  )
  coefficients <- if (!is.null(custom_solver)) {
    returned <- capture_solver_warnings(custom_solver$fit(
      X,
      y,
      weights,
      taus,
      list(
        precondition = precondition,
        onestep_first_solver = onestep_first_solver,
        onestep_bandwidth = onestep_bandwidth
      )
    ))
    if (is.matrix(returned) || is.numeric(returned)) {
      returned <- list(coefficients = returned)
    }
    if (!is.list(returned) || is.null(returned$coefficients)) {
      stop("custom QR solver must return coefficients", call. = FALSE)
    }
    fit <- returned
    fit$flag <- fit$flag %||% fit$convergence_flag
    fit$nit <- fit$nit %||% fit$iterations
    if (length(fit$warnings %||% character())) {
      solver_warnings <- c(solver_warnings, as.character(fit$warnings))
    }
    fit$coefficients
  } else if (solver %in% c("br", "fn", "pfn")) {
    single_fits <- lapply(taus, function(tau) {
      capture_solver_warnings(
        single_tau_qr(weighted_X, weighted_y, tau, solver)
      )
    })
    vapply(single_fits, function(single_fit) {
      as.numeric(single_fit$coefficients)
    }, numeric(ncol(X)))
  } else if (solver == "qfnb") {
    fit <- capture_solver_warnings(
      quantreg::rq.fit.qfnb(weighted_X, weighted_y, tau = taus)
    )
    fit$coefficients
  } else if (solver == "pfnb") {
    fit <- capture_solver_warnings(
      quantreg::rq.fit.pfnb(
        weighted_X,
        weighted_y,
        tau = taus,
        m0 = pfnb_m0
      )
    )
    fit$coefficients
  } else if (solver %in% c("proqreg", "profn")) {
    process_method <- if (solver == "proqreg") "br" else "fn"
    process_attempt <- tryCatch(
      capture_solver_warnings(
        quantreg::rq.fit.ppro(
          conditioned_X,
          y,
          tau = taus,
          weights = weights,
          pmethod = process_method
        )
      ),
      error = identity
    )
    suspect_warning <- any(grepl(
      "singular|Too many fixups|Error info",
      solver_warnings,
      ignore.case = TRUE
    ))
    if (inherits(process_attempt, "error") || suspect_warning) {
      process_fallback_taus <- taus
      process_fallback_solver <- process_method
      process_fallback_reason <- if (inherits(process_attempt, "error")) {
        conditionMessage(process_attempt)
      } else {
        "rq.fit.ppro reported numerical fixups or a possibly singular reduced design"
      }
      solver_warnings <- c(
        solver_warnings,
        paste0(
          solver, " validation fallback: refitted all conditional quantiles with ",
          process_method, "; reason: ", process_fallback_reason
        )
      )
      single_fits <- lapply(taus, function(tau) {
        capture_solver_warnings(
          single_tau_qr(weighted_X, weighted_y, tau, process_method)
        )
      })
      vapply(single_fits, function(single_fit) {
        as.numeric(single_fit$coefficients)
      }, numeric(ncol(X)))
    } else {
      fit <- process_attempt
      fit$coefficients
    }
  } else if (solver == "onestep") {
    fit <- capture_solver_warnings(
      fit_qr_onestep(
        conditioned_X,
        y,
        weights,
        taus,
        first_solver = onestep_first_solver,
        bandwidth_method = onestep_bandwidth
      )
    )
    if (length(fit$warnings)) {
      solver_warnings <- c(solver_warnings, fit$warnings)
    }
    fit$coefficients
  } else if (solver == "cuda_admm") {
    if (is.null(gpu_control)) {
      stop("cuda_admm requires gpu_control from cf_control()", call. = FALSE)
    }
    fit <- capture_solver_warnings(fit_qr_cuda_admm(
      X, y, weights, taus,
      list(precondition = precondition, control = gpu_control)
    ))
    if (length(fit$warnings)) {
      solver_warnings <- c(solver_warnings, fit$warnings)
    }
    fit$coefficients
  } else {
    stop("unimplemented QR solver: ", solver, call. = FALSE)
  }
  coefficients <- normalize_coefficient_matrix(coefficients, ncol(X), taus)
  conditioned_coefficients <- if (solver == "cuda_admm") {
    fit$conditioned_coefficients
  } else {
    coefficients
  }
  coefficients <- preconditioner %*% coefficients
  if (anyNA(coefficients) || any(!is.finite(coefficients))) {
    stop(solver, " returned invalid coefficients", call. = FALSE)
  }
  flag <- if (is.null(fit)) rep(NA_integer_, length(taus)) else {
    raw_flag <- fit$flag
    if (is.null(raw_flag)) raw_flag <- fit$info
    if (is.null(raw_flag)) rep(NA_integer_, length(taus)) else {
      rep(as.integer(raw_flag), length.out = length(taus))
    }
  }
  allow_cuda_audit <- solver == "cuda_admm" &&
    isTRUE(gpu_control$gpu_qr_allow_nonconvergence)
  if (any(!is.na(flag) & flag != 0L) && !allow_cuda_audit) {
    stop(
      solver, " convergence failure: ",
      paste(unique(flag[!is.na(flag) & flag != 0L]), collapse = ","),
      call. = FALSE
    )
  }
  structure(list(
    model = "qr",
    solver = solver,
    solver_requested = solver_requested,
    solver_implementation = qr_solver_registry()$implementation[
      match(solver, qr_solver_registry()$solver)
    ],
    solver_exact = qr_solver_registry()$exact[
      match(solver, qr_solver_registry()$solver)
    ],
    solver_process_aware = qr_solver_registry()$process_aware[
      match(solver, qr_solver_registry()$solver)
    ],
    taus = taus,
    coefficients = coefficients,
    conditioned_coefficients = conditioned_coefficients,
    convergence_flag = flag,
    convergence_diagnostics_available = any(!is.na(flag)),
    cuda_admm_converged = if (solver == "cuda_admm") {
      fit$cuda_converged
    } else {
      NULL
    },
    cuda_admm_primal_residual = if (solver == "cuda_admm") {
      fit$primal_residual
    } else {
      NULL
    },
    cuda_admm_dual_residual = if (solver == "cuda_admm") {
      fit$dual_residual
    } else {
      NULL
    },
    iterations = if (!is.null(single_fits)) {
      lapply(single_fits, function(x) x$nit)
    } else if (is.null(fit)) NULL else fit$nit,
    warnings = unique(solver_warnings),
    preprocessing_sample_size = if (solver == "pfnb") {
      pfnb_m0
    } else if (solver %in% c("proqreg", "profn")) {
      if (length(taus) > 1L) {
        as.integer(round(nrow(X) * sqrt(ncol(X)) * max(diff(taus))))
      } else {
        NA_integer_
      }
    } else {
      NA_integer_
    },
    process_fallback_taus = process_fallback_taus,
    process_fallback_solver = process_fallback_solver,
    process_fallback_reason = process_fallback_reason,
    onestep_first_solver_requested = if (solver == "onestep") {
      fit$first_solver_requested
    } else {
      NA_character_
    },
    onestep_first_solver_resolved = if (solver == "onestep") {
      fit$first_solver_resolved
    } else {
      NA_character_
    },
    onestep_fallback_taus = if (solver == "onestep") {
      fit$fallback_taus
    } else {
      numeric()
    },
    onestep_inverse_jacobians = if (solver == "onestep") {
      fit$inverse_jacobians
    } else {
      NULL
    },
    onestep_bandwidth = if (solver == "onestep") {
      fit$bandwidth_method
    } else {
      NA_character_
    },
    onestep_traversal = if (solver == "onestep") {
      fit$traversal
    } else {
      integer()
    },
    onestep_implementation_version = if (solver == "onestep") {
      "qrprocess_1.1.3_native_r"
    } else {
      NA_character_
    },
    stata_source_version = if (solver == "onestep") {
      "qrprocess 1.1.3"
    } else {
      NA_character_
    },
    stata_source_commit = if (solver == "onestep") {
      "ec56830ef9c84ce54411ad59c5ce94535847d9df"
    } else {
      NA_character_
    },
    precondition_requested = if (solver == "cuda_admm") {
      fit$precondition_requested
    } else isTRUE(precondition),
    preconditioned = if (solver == "cuda_admm") {
      fit$preconditioned
    } else precondition_effective,
    preconditioning_method = if (solver == "cuda_admm") {
      fit$preconditioning_method
    } else preconditioned$method,
    preconditioning_gram_rcond = if (solver == "cuda_admm") {
      fit$preconditioning_gram_rcond
    } else preconditioned$gram_rcond,
    preconditioning_condition_estimate = if (solver == "cuda_admm") {
      fit$preconditioning_condition_estimate
    } else preconditioned$condition_estimate,
    preconditioning_fallback_reason = if (solver == "cuda_admm") {
      fit$preconditioning_fallback_reason
    } else preconditioned$fallback_reason,
    preconditioning_matrix = if (solver == "cuda_admm") {
      fit$preconditioning_matrix
    } else preconditioner
  ), class = c("cf_qr_fit", "cf_conditional_fit"))
}
