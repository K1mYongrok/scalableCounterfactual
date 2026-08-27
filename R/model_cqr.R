# Multi-step CQR follows Chernozhukov and Hong (2002),
# doi:10.1198/016214502388618663; exact QR stages are delegated to quantreg.
validate_cqr_solver <- function(solver) {
  allowed <- c("br", "fn", "pfn", "qfnb", "pfnb", "auto")
  match.arg(solver, allowed)
}

fit_cqr_left_process <- function(
    X, y, weights, censoring, taus, solver, nsteps, first_cut, later_cut,
    precondition, dr_backend, dr_maxit, dr_tolerance) {
  uncensored <- as.numeric(y > censoring)
  if (length(unique(uncensored)) < 2L) {
    stop("CQR requires both censored and uncensored observations", call. = FALSE)
  }
  selection_fit <- fit_binary_threshold(
    X, uncensored, weights, "logit", backend = dr_backend,
    maxit = dr_maxit, tolerance = dr_tolerance
  )
  probabilities <- stats::plogis(drop(X %*% selection_fit$coefficients))
  p <- ncol(X)
  coefficients <- matrix(NA_real_, p, length(taus))
  selection_sizes <- matrix(
    NA_integer_, nrow = length(taus), ncol = nsteps - 1L,
    dimnames = list(paste0("tau_", signif(taus, 6)), paste0("step_", 2:nsteps))
  )
  warnings <- selection_fit$warnings
  resolved_solvers <- character(length(taus))

  for (j in seq_along(taus)) {
    tau <- taus[[j]]
    eligible <- probabilities > 1 - tau
    if (sum(eligible) <= p) {
      stop(
        "CQR first-stage selection has too few observations at tau=",
        signif(tau, 6), call. = FALSE
      )
    }
    cutoff <- weighted_quantile(
      probabilities[eligible], weights[eligible], first_cut
    )
    selected <- probabilities >= cutoff
    selection_sizes[j, 1L] <- sum(selected)
    if (sum(selected) <= p) {
      stop(
        "CQR second-step design is unidentified at tau=", signif(tau, 6),
        call. = FALSE
      )
    }
    qr_fit <- fit_weighted_qr(
      X[selected, , drop = FALSE], y[selected], weights[selected], tau,
      solver = solver, precondition = precondition
    )
    beta <- qr_fit$coefficients[, 1L]
    resolved_solvers[[j]] <- qr_fit$solver
    warnings <- c(warnings, qr_fit$warnings)

    if (nsteps >= 3L) {
      for (step in 3:nsteps) {
        gap <- drop(X %*% beta) - censoring
        eligible <- gap > 0
        if (sum(eligible) <= p) {
          stop(
            "CQR iterative selection has too few uncensored predictions at tau=",
            signif(tau, 6), ", step=", step, call. = FALSE
          )
        }
        margin <- weighted_quantile(gap[eligible], weights[eligible], later_cut)
        selected <- gap >= margin
        selection_sizes[j, step - 1L] <- sum(selected)
        if (sum(selected) <= p) {
          stop(
            "CQR iterative design is unidentified at tau=", signif(tau, 6),
            ", step=", step, call. = FALSE
          )
        }
        qr_fit <- fit_weighted_qr(
          X[selected, , drop = FALSE], y[selected], weights[selected], tau,
          solver = solver, precondition = precondition
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
    warnings = unique(warnings)
  )
}

fit_weighted_cqr <- function(
    X, y, weights, censoring, taus, solver = "auto", right = FALSE,
    nsteps = 3L, first_cut = 0.1, later_cut = 0.05,
    precondition = TRUE, dr_backend = "auto", dr_maxit = 100L,
    dr_tolerance = 1e-8) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  weights <- normalize_weights(weights)
  censoring <- as.numeric(censoring)
  taus <- sort(unique(as.numeric(taus)))
  solver <- validate_cqr_solver(solver)
  if (length(censoring) != length(y) || nrow(X) != length(y)) {
    stop("CQR inputs have incompatible sizes", call. = FALSE)
  }
  if (any(!is.finite(censoring))) {
    stop("CQR censoring points must be finite", call. = FALSE)
  }
  transformed <- if (isTRUE(right)) {
    fit_cqr_left_process(
      X, -y, weights, -censoring, 1 - taus, solver, nsteps,
      first_cut, later_cut, precondition, dr_backend, dr_maxit, dr_tolerance
    )
  } else {
    fit_cqr_left_process(
      X, y, weights, censoring, taus, solver, nsteps,
      first_cut, later_cut, precondition, dr_backend, dr_maxit, dr_tolerance
    )
  }
  coefficients <- if (isTRUE(right)) -transformed$coefficients else {
    transformed$coefficients
  }
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
    selection_sizes = transformed$selection_sizes,
    selection_backend = transformed$selection_fit$backend,
    warnings = transformed$warnings
  ), class = c("cf_cqr_fit", "cf_conditional_fit"))
}
