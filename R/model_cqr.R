# Multi-step CQR follows Chernozhukov and Hong (2002),
# doi:10.1198/016214502388618663; exact QR stages are delegated to quantreg.
validate_cqr_solver <- function(solver) {
  allowed <- c("br", "fn", "pfn", "qfnb", "pfnb", "auto")
  match.arg(solver, allowed)
}

cqr_design_identified <- function(X, selected, quantile_frequency, p) {
  if (sum(quantile_frequency[selected]) <= p || sum(selected) < p) {
    return(FALSE)
  }
  qr(X[selected, , drop = FALSE], LAPACK = FALSE)$rank == p
}

fit_cqr_left_process <- function(
    X, y, weights, censoring, taus, solver, nsteps, first_cut, later_cut,
    precondition, dr_backend, dr_maxit, dr_tolerance,
    quantile_frequency) {
  uncensored <- as.numeric(y > censoring)
  if (length(unique(uncensored)) < 2L) {
    stop("CQR requires both censored and uncensored observations", call. = FALSE)
  }
  initial_selection_fit <- fit_binary_threshold(
    X, uncensored, weights, "logit", backend = dr_backend,
    maxit = dr_maxit, tolerance = dr_tolerance
  )
  selection_fit <- initial_selection_fit
  selection_fallback <- FALSE
  selection_problem <- !isTRUE(selection_fit$converged) ||
    isTRUE(selection_fit$boundary)
  if (selection_problem && !identical(selection_fit$backend, "glm")) {
    selection_fit <- fit_binary_threshold(
      X, uncensored, weights, "logit", backend = "glm",
      maxit = dr_maxit, tolerance = dr_tolerance
    )
    selection_fallback <- TRUE
  }
  if (!isTRUE(selection_fit$converged)) {
    stop(
      "CQR first-stage selection model failed to converge",
      if (length(selection_fit$iterations) == 1L &&
          is.finite(selection_fit$iterations)) {
        paste0(" after ", selection_fit$iterations, " iteration(s)")
      } else {
        ""
      },
      call. = FALSE
    )
  }
  if (isTRUE(selection_fit$boundary)) {
    stop(
      "CQR first-stage selection model reached boundary probabilities or ",
      "separation; increase overlap or revise the selection specification",
      call. = FALSE
    )
  }
  probabilities <- stats::plogis(drop(X %*% selection_fit$coefficients))
  p <- ncol(X)
  coefficients <- matrix(NA_real_, p, length(taus))
  selection_sizes <- matrix(
    NA_integer_, nrow = length(taus), ncol = nsteps - 1L,
    dimnames = list(paste0("tau_", signif(taus, 6)), paste0("step_", 2:nsteps))
  )
  warnings <- unique(c(initial_selection_fit$warnings, selection_fit$warnings))
  if (selection_fallback) {
    warnings <- c(
      warnings,
      paste0(
        "CQR selection used glm fallback after ",
        initial_selection_fit$backend,
        " reported nonconvergence, unavailable diagnostics, or a boundary fit"
      )
    )
  }
  resolved_solvers <- character(length(taus))

  for (j in seq_along(taus)) {
    tau <- taus[[j]]
    eligible <- probabilities > 1 - tau
    if (sum(quantile_frequency[eligible]) <= p) {
      stop(
        "CQR first-stage selection has too few observations at tau=",
        signif(tau, 6), call. = FALSE
      )
    }
    cutoff <- weighted_quantile(
      probabilities[eligible], weights[eligible], first_cut,
      normalization_n = sum(quantile_frequency[eligible])
    )
    selected <- probabilities >= cutoff
    selection_sizes[j, 1L] <- as.integer(round(sum(
      quantile_frequency[selected]
    )))
    if (!cqr_design_identified(X, selected, quantile_frequency, p)) {
      stop(
        "CQR second-step design is unidentified at tau=", signif(tau, 6),
        call. = FALSE
      )
    }
    qr_fit <- fit_weighted_qr(
      X[selected, , drop = FALSE], y[selected],
      weights[selected] / quantile_frequency[selected], tau,
      solver = solver, precondition = precondition,
      frequency = quantile_frequency[selected]
    )
    beta <- qr_fit$coefficients[, 1L]
    resolved_solvers[[j]] <- qr_fit$solver
    warnings <- c(warnings, qr_fit$warnings)

    if (nsteps >= 3L) {
      for (step in 3:nsteps) {
        gap <- drop(X %*% beta) - censoring
        eligible <- gap > 0
        if (sum(quantile_frequency[eligible]) <= p) {
          stop(
            "CQR iterative selection has too few uncensored predictions at tau=",
            signif(tau, 6), ", step=", step, call. = FALSE
          )
        }
        margin <- weighted_quantile(
          gap[eligible], weights[eligible], later_cut,
          normalization_n = sum(quantile_frequency[eligible])
        )
        selected <- gap >= margin
        selection_sizes[j, step - 1L] <- as.integer(round(sum(
          quantile_frequency[selected]
        )))
        if (!cqr_design_identified(X, selected, quantile_frequency, p)) {
          stop(
            "CQR iterative design is unidentified at tau=", signif(tau, 6),
            ", step=", step, call. = FALSE
          )
        }
        qr_fit <- fit_weighted_qr(
          X[selected, , drop = FALSE], y[selected],
          weights[selected] / quantile_frequency[selected], tau,
          solver = solver, precondition = precondition,
          frequency = quantile_frequency[selected]
        )
        beta <- qr_fit$coefficients[, 1L]
        resolved_solvers[[j]] <- qr_fit$solver
        warnings <- c(warnings, qr_fit$warnings)
      }
    }
    coefficients[, j] <- beta
  }

  list(
    coefficients = coefficients,
    selection_sizes = selection_sizes,
    resolved_solvers = resolved_solvers,
    selection_fit = selection_fit,
    initial_selection_fit = initial_selection_fit,
    selection_fallback = selection_fallback,
    warnings = unique(warnings)
  )
}

fit_weighted_cqr <- function(
    X, y, weights, censoring, taus, solver = "auto", right = FALSE,
    nsteps = 3L, first_cut = 0.1, later_cut = 0.05,
    precondition = TRUE, dr_backend = "auto", dr_maxit = 100L,
    dr_tolerance = 1e-8, quantile_frequency = rep(1, length(y))) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  weights <- normalize_weights(weights)
  censoring <- as.numeric(censoring)
  taus <- sort(unique(as.numeric(taus)))
  solver <- validate_cqr_solver(solver)
  if (length(censoring) != length(y) || nrow(X) != length(y)) {
    stop("CQR inputs have incompatible sizes", call. = FALSE)
  }
  if (!is.numeric(quantile_frequency) ||
      length(quantile_frequency) != length(y) ||
      any(!is.finite(quantile_frequency)) || any(quantile_frequency <= 0)) {
    stop("quantile_frequency must contain one finite positive value per row",
         call. = FALSE)
  }
  if (any(!is.finite(censoring))) {
    stop("CQR censoring points must be finite", call. = FALSE)
  }
  comparison_scale <- pmax(1, abs(y), abs(censoring))
  comparison_tolerance <- sqrt(.Machine$double.eps) * comparison_scale
  inconsistent <- if (isTRUE(right)) {
    y > censoring + comparison_tolerance
  } else {
    y < censoring - comparison_tolerance
  }
  if (any(inconsistent)) {
    direction <- if (isTRUE(right)) "above" else "below"
    stop(
      "CQR found ", sum(inconsistent), " observed outcome(s) ", direction,
      " the declared ", if (isTRUE(right)) "right" else "left",
      " censoring point",
      call. = FALSE
    )
  }
  transformed <- if (isTRUE(right)) {
    fit_cqr_left_process(
      X, -y, weights, -censoring, 1 - taus, solver, nsteps,
      first_cut, later_cut, precondition, dr_backend, dr_maxit, dr_tolerance,
      quantile_frequency
    )
  } else {
    fit_cqr_left_process(
      X, y, weights, censoring, taus, solver, nsteps,
      first_cut, later_cut, precondition, dr_backend, dr_maxit, dr_tolerance,
      quantile_frequency
    )
  }
  coefficients <- if (isTRUE(right)) -transformed$coefficients else {
    transformed$coefficients
  }
  selection_sizes <- transformed$selection_sizes
  rownames(selection_sizes) <- paste0("tau_", signif(taus, 6))
  structure(list(
    model = "cqr",
    solver = if (length(unique(transformed$resolved_solvers)) == 1L) {
      unique(transformed$resolved_solvers)
    } else {
      paste(unique(transformed$resolved_solvers), collapse = ",")
    },
    solver_requested = solver,
    solver_exact = TRUE,
    solver_process_aware = FALSE,
    taus = taus,
    coefficients = coefficients,
    right = isTRUE(right),
    nsteps = nsteps,
    first_cut = first_cut,
    later_cut = later_cut,
    selection_sizes = selection_sizes,
    selection_backend = transformed$selection_fit$backend,
    selection_fit = transformed$selection_fit,
    selection_diagnostics = list(
      converged = transformed$selection_fit$converged,
      boundary = transformed$selection_fit$boundary,
      iterations = transformed$selection_fit$iterations,
      backend = transformed$selection_fit$backend,
      fallback_used = transformed$selection_fallback,
      initial_backend = transformed$initial_selection_fit$backend,
      initial_converged = transformed$initial_selection_fit$converged,
      initial_boundary = transformed$initial_selection_fit$boundary,
      warnings = transformed$warnings
    ),
    warnings = transformed$warnings
  ), class = c("cf_cqr_fit", "cf_conditional_fit"))
}
