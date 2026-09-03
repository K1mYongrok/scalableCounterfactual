# Distribution-regression counterfactuals follow Chernozhukov,
# Fernandez-Val, and Melly (2013), doi:10.3982/ECTA10582.
.cf_dr_worker_state <- new.env(parent = emptyenv())

select_distribution_thresholds <- function(y, nreg) {
  thresholds <- sort(unique(as.numeric(y)))
  if (length(thresholds) <= nreg || nreg == -1L) return(thresholds)
  interior_count <- nreg - 2L
  grid <- (seq_len(interior_count) - 0.5) / interior_count
  interior <- floor(grid * length(thresholds))
  interior <- pmax(1L, pmin(length(thresholds), interior))
  unique(thresholds[c(1L, interior, length(thresholds))])
}

resolve_dr_backend <- function(backend) {
  backend <- match.arg(backend, supported_dr_backends())
  if (backend == "auto") {
    return(if (requireNamespace("fastglm", quietly = TRUE)) "fastglm" else "glm")
  }
  backend
}

dr_link <- function(model) {
  switch(
    model,
    logit = "logit",
    probit = "probit",
    cloglog = "cloglog",
    stop("unsupported distribution-regression model: ", model, call. = FALSE)
  )
}

capture_backend_warnings <- function(expression) {
  messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      message <- conditionMessage(condition)
      if (!grepl("non-integer #successes", message, fixed = TRUE)) {
        messages <<- c(messages, message)
      }
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(messages))
}

fit_binary_threshold <- function(
    X, response, weights, model, backend = "auto", start = NULL,
    maxit = 100L, tolerance = 1e-8, constant_direction = NULL) {
  backend_resolved <- resolve_dr_backend(backend)
  if (backend_resolved == "cuda") {
    stop("cuda is a process-level DR backend", call. = FALSE)
  }
  family <- stats::binomial(link = dr_link(model))
  if (length(unique(response)) == 1L) {
    probability <- if (response[[1L]] == 1) 1 - 1e-8 else 1e-8
    link_value <- family$linkfun(probability)
    if (is.null(constant_direction)) {
      stop(
        "a constant binary response requires an identified intercept direction",
        call. = FALSE
      )
    } else {
      coefficients <- as.numeric(constant_direction) * link_value
    }
    return(list(
      coefficients = coefficients,
      converged = TRUE,
      boundary = FALSE,
      iterations = 0L,
      backend = backend_resolved,
      warnings = character()
    ))
  }

  fitted <- if (backend_resolved == "glm") {
    capture_backend_warnings(stats::glm.fit(
      x = X,
      y = response,
      weights = weights,
      family = family,
      intercept = FALSE,
      start = start,
      control = stats::glm.control(epsilon = tolerance, maxit = maxit)
    ))
  } else if (backend_resolved == "fastglm") {
    require_optional_backend("fastglm", "fastglm")
    capture_backend_warnings(fastglm::fastglmPure(
      X,
      response,
      family = family,
      weights = weights,
      start = start,
      maxit = maxit,
      tol = tolerance
    ))
  } else if (backend_resolved == "speedglm") {
    require_optional_backend("speedglm", "speedglm")
    capture_backend_warnings(speedglm::speedglm.wfit(
      y = response,
      X = X,
      intercept = FALSE,
      weights = weights,
      family = family,
      start = start,
      acc = tolerance,
      maxit = maxit,
      method = "eigen"
    ))
  } else {
    entry <- custom_registry_entry("distribution", backend_resolved)
    if (is.null(entry)) {
      stop("unregistered distribution backend: ", backend_resolved,
           call. = FALSE)
    }
    custom_warnings <- character()
    returned <- withCallingHandlers(
      entry$fit(
        X, response, weights, model, start, maxit, tolerance
      ),
      warning = function(condition) {
        custom_warnings <<- c(custom_warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    if (is.numeric(returned)) returned <- list(coefficients = returned)
    if (!is.list(returned) || is.null(returned$coefficients)) {
      stop("custom distribution backend must return coefficients",
           call. = FALSE)
    }
    returned$coefficients <- as.numeric(returned$coefficients)
    if (length(returned$coefficients) != ncol(X) ||
        any(!is.finite(returned$coefficients))) {
      stop("custom distribution backend returned invalid coefficients",
           call. = FALSE)
    }
    if (!is.null(returned$converged) &&
        (!is.logical(returned$converged) || length(returned$converged) != 1L ||
         is.na(returned$converged))) {
      stop("custom distribution backend converged must be one TRUE/FALSE value",
           call. = FALSE)
    }
    if (!is.null(returned$boundary) &&
        (!is.logical(returned$boundary) || length(returned$boundary) != 1L ||
         is.na(returned$boundary))) {
      stop("custom distribution backend boundary must be one TRUE/FALSE value",
           call. = FALSE)
    }
    if (!is.null(returned$iterations) &&
        (length(returned$iterations) != 1L ||
         !is.finite(returned$iterations) || returned$iterations < 0 ||
         returned$iterations != floor(returned$iterations))) {
      stop("custom distribution backend iterations must be one nonnegative integer",
           call. = FALSE)
    }
    if (!is.null(returned$fitted.values)) {
      returned$fitted.values <- as.numeric(returned$fitted.values)
      if (length(returned$fitted.values) != nrow(X) ||
          any(!is.finite(returned$fitted.values)) ||
          any(returned$fitted.values < 0 | returned$fitted.values > 1)) {
        stop("custom distribution backend returned invalid fitted.values",
             call. = FALSE)
      }
    }
    fit <- list(
      coefficients = returned$coefficients,
      fitted.values = returned$fitted.values %||%
        family$linkinv(drop(X %*% returned$coefficients)),
      converged = returned$converged %||% NA,
      boundary = returned$boundary %||% FALSE,
      iter = returned$iterations %||% NA_integer_,
      custom_warnings = unique(c(
        custom_warnings, as.character(returned$warnings %||% character())
      ))
    )
    fitted <- list(value = fit, warnings = fit$custom_warnings)
  }

  fit <- fitted$value
  coefficients <- as.numeric(fit$coefficients)
  linear_predictor <- drop(X %*% coefficients)
  probabilities <- if (!is.null(fit$fitted.values)) {
    as.numeric(fit$fitted.values)
  } else {
    as.numeric(family$linkinv(linear_predictor))
  }
  if (length(probabilities) != nrow(X) || any(!is.finite(probabilities)) ||
      any(probabilities < 0 | probabilities > 1)) {
    stop(backend_resolved, " backend returned invalid fitted probabilities",
         call. = FALSE)
  }
  converged <- if (backend_resolved == "speedglm") {
    isTRUE(fit$convergence)
  } else if (is.na(fit$converged)) {
    NA
  } else {
    isTRUE(fit$converged)
  }
  # Raw coefficient magnitudes and isolated extreme fitted probabilities are
  # not separation diagnostics.  Complete and quasi-complete separation do
  # imply that every fitted signed margin is nonnegative, with at least one
  # strictly positive margin.  A scale-aware tolerance keeps this test stable
  # under harmless covariate rescaling and floating-point roundoff.
  signed_margin <- (2 * response - 1) * linear_predictor
  margin_tolerance <- max(tolerance, 64 * .Machine$double.eps) *
    max(1, max(abs(linear_predictor)))
  separated <- all(signed_margin >= -margin_tolerance) &&
    any(signed_margin > margin_tolerance)
  boundary <- isTRUE(fit$boundary) || separated
  list(
    coefficients = coefficients,
    converged = converged,
    boundary = boundary,
    iterations = if (is.null(fit$iter)) NA_integer_ else as.integer(fit$iter),
    backend = backend_resolved,
    warnings = fitted$warnings
  )
}

initialize_dr_worker <- function(
    X, y, weights, thresholds, model, backend, maxit, tolerance,
    constant_direction) {
  .cf_dr_worker_state$X <- X
  .cf_dr_worker_state$y <- y
  .cf_dr_worker_state$weights <- weights
  .cf_dr_worker_state$thresholds <- thresholds
  .cf_dr_worker_state$model <- model
  .cf_dr_worker_state$backend <- backend
  .cf_dr_worker_state$maxit <- maxit
  .cf_dr_worker_state$tolerance <- tolerance
  .cf_dr_worker_state$constant_direction <- constant_direction
  invisible(NULL)
}

fit_dr_worker_threshold <- function(index) {
  state <- .cf_dr_worker_state
  fit_binary_threshold(
    state$X,
    as.numeric(state$y <= state$thresholds[[index]]),
    state$weights,
    state$model,
    state$backend,
    start = NULL,
    maxit = state$maxit,
    tolerance = state$tolerance,
    constant_direction = state$constant_direction
  )
}

fit_dr_thresholds_sequential <- function(
    X, y, weights, thresholds, model, backend, warm_start, maxit, tolerance,
    constant_direction) {
  fits <- vector("list", length(thresholds))
  if (!warm_start) {
    return(lapply(thresholds, function(threshold) {
      fit_binary_threshold(
        X, as.numeric(y <= threshold), weights, model, backend,
        maxit = maxit, tolerance = tolerance,
        constant_direction = constant_direction
      )
    }))
  }

  center <- which.min(abs(thresholds - stats::median(y)))
  traversal <- c(seq.int(center, length(thresholds)),
                 if (center > 1L) seq.int(center - 1L, 1L) else integer())
  previous <- NULL
  previous_direction <- 1L
  for (position in seq_along(traversal)) {
    index <- traversal[[position]]
    direction <- if (index >= center) 1L else -1L
    if (direction != previous_direction) {
      previous <- fits[[center]]$coefficients
    }
    fit <- fit_binary_threshold(
      X,
      as.numeric(y <= thresholds[[index]]),
      weights,
      model,
      backend,
      start = previous,
      maxit = maxit,
      tolerance = tolerance,
      constant_direction = constant_direction
    )
    fits[[index]] <- fit
    previous <- if (isTRUE(fit$converged)) fit$coefficients else NULL
    previous_direction <- direction
  }
  fits
}

fit_lpm_process <- function(X, y, weights, thresholds, linear_backend) {
  solver <- make_weighted_linear_solver(X, weights, linear_backend)
  fits <- lapply(thresholds, function(threshold) {
    fitted <- solver$fit(as.numeric(y <= threshold))
    list(
      coefficients = fitted$coefficients,
      converged = TRUE,
      boundary = FALSE,
      iterations = 1L,
      backend = solver$backend,
      warnings = character()
    )
  })
  list(fits = fits, warnings = solver$warnings())
}

refit_cuda_dr_failures <- function(
    process, X, y, weights, thresholds, model, maxit, tolerance,
    constant_direction) {
  fallback <- which(!process$converged)
  if (!length(fallback)) {
    return(list(process = process, fallback = integer()))
  }
  for (index in fallback) {
    cpu_fit <- fit_binary_threshold(
      X, as.numeric(y <= thresholds[[index]]), weights, model,
      backend = "glm", maxit = maxit, tolerance = tolerance,
      constant_direction = constant_direction
    )
    process$coefficients[, index] <- cpu_fit$coefficients
    process$converged[[index]] <- cpu_fit$converged
    process$boundary[[index]] <- cpu_fit$boundary
    process$iterations[[index]] <- cpu_fit$iterations
    process$backend[[index]] <- "glm"
    process$fallback_reason[[index]] <- "cuda_nonconvergence"
  }
  list(process = process, fallback = fallback)
}

fit_distribution_regression <- function(
    X, y, weights, model, nreg, dr_backend = "auto",
    linear_backend = "auto", dr_workers = 1L, warm_start = TRUE,
    maxit = 100L, tolerance = 1e-8, precondition = TRUE,
    control = NULL) {
  thresholds <- select_distribution_thresholds(y, nreg)
  execution_warnings <- character()
  transform <- diag(ncol(X))
  preconditioned <- list(
    method = "none",
    gram_rcond = NA_real_,
    condition_estimate = NA_real_,
    fallback_reason = NA_character_
  )

  if (model == "lpm") {
    process <- fit_lpm_process(
      X, y, weights, thresholds, linear_backend
    )
    fits <- process$fits
    execution_warnings <- process$warnings
    backend_requested <- linear_backend
  } else {
    preconditioned <- precondition_qr_design(
      X * sqrt(normalize_weights(weights)), precondition
    )
    transform <- preconditioned$transform
    fit_X <- X %*% transform
    intercept_index <- match("(Intercept)", colnames(X))
    if (is.na(intercept_index)) {
      stop(
        model, " distribution regression requires an intercept because its ",
        "threshold grid includes a constant-response endpoint",
        call. = FALSE
      )
    }
    constant_original <- numeric(ncol(X))
    constant_original[[intercept_index]] <- 1
    constant_direction <- as.numeric(solve(transform, constant_original))
    backend_requested <- dr_backend
    backend <- resolve_dr_backend(dr_backend)
    if (backend == "cuda") {
      if (is.null(control)) {
        stop("dr_backend='cuda' requires a cf_control object", call. = FALSE)
      }
      if (dr_workers > 1L) {
        stop("dr_backend='cuda' requires dr_workers=1", call. = FALSE)
      }
      if (warm_start) {
        execution_warnings <- c(
          execution_warnings,
          "DR warm starts are disabled by the batched CUDA backend"
        )
      }
      process <- fit_dr_process_cuda(
        fit_X, y, weights, thresholds, model, maxit, tolerance,
        constant_direction, control
      )
      fallback_result <- refit_cuda_dr_failures(
        process, fit_X, y, weights, thresholds, model, maxit, tolerance,
        constant_direction
      )
      process <- fallback_result$process
      fallback <- fallback_result$fallback
      if (length(fallback)) {
        execution_warnings <- c(
          execution_warnings,
          paste0(
            "CUDA DR used stats::glm.fit fallback at threshold(s): ",
            paste(signif(thresholds[fallback], 6), collapse = ", ")
          )
        )
      }
      fits <- lapply(seq_along(thresholds), function(index) list(
        coefficients = process$coefficients[, index],
        converged = process$converged[[index]],
        boundary = process$boundary[[index]],
        iterations = process$iterations[[index]],
        backend = process$backend[[index]],
        fallback_reason = process$fallback_reason[[index]],
        warnings = character()
      ))
    } else if (dr_workers > 1L) {
      if (warm_start) {
        execution_warnings <- c(
          execution_warnings,
          "DR warm starts are disabled when dr_workers > 1"
        )
      }
      cluster <- parallel::makeCluster(min(dr_workers, length(thresholds)))
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      package_worker_init(cluster)
      parallel::clusterCall(
        cluster,
        function(X, y, weights, thresholds, model, backend, maxit, tolerance,
                 constant_direction) {
          initialize <- get(
            "initialize_dr_worker",
            envir = asNamespace("scalableCounterfactual")
          )
          initialize(
            X, y, weights, thresholds, model, backend, maxit, tolerance,
            constant_direction
          )
        },
        fit_X, y, weights, thresholds, model, backend, maxit, tolerance,
        constant_direction
      )
      fits <- parallel::parLapplyLB(cluster, seq_along(thresholds), function(i) {
        worker <- get(
          "fit_dr_worker_threshold",
          envir = asNamespace("scalableCounterfactual")
        )
        worker(i)
      })
    } else {
      fits <- fit_dr_thresholds_sequential(
        fit_X, y, weights, thresholds, model, backend, warm_start,
        maxit, tolerance, constant_direction
      )
    }
  }

  coefficients <- vapply(fits, `[[`, numeric(ncol(X)), "coefficients")
  coefficients <- matrix(
    as.numeric(coefficients), nrow = ncol(X), ncol = length(fits)
  )
  if (model != "lpm") coefficients <- transform %*% coefficients
  if (anyNA(coefficients) || any(!is.finite(coefficients))) {
    stop(model, " distribution regression returned invalid coefficients",
         call. = FALSE)
  }
  converged <- vapply(fits, `[[`, logical(1L), "converged")
  if (any(!converged, na.rm = TRUE)) {
    parallel_hint <- if (model != "lpm" && dr_workers > 1L) {
      "; parallel threshold fits do not use warm starts; try dr_workers=1 or a different backend"
    } else {
      ""
    }
    stop(
      model, " distribution regression failed to converge at threshold(s): ",
      paste(signif(thresholds[!converged], 6), collapse = ", "),
      parallel_hint,
      call. = FALSE
    )
  }
  boundary <- vapply(fits, `[[`, logical(1L), "boundary")
  fit_warnings <- unique(c(execution_warnings, unlist(Map(
    function(fit, threshold) {
      messages <- fit$warnings
      if (fit$boundary) {
        messages <- c(messages, "boundary probabilities or separation detected")
      }
      if (!length(messages)) return(character())
      paste0("threshold=", signif(threshold, 6), ": ", messages)
    },
    fits,
    thresholds
  ))))
  threshold_backend <- vapply(fits, `[[`, character(1L), "backend")
  threshold_fallback_reason <- vapply(fits, function(fit) {
    fit$fallback_reason %||% NA_character_
  }, character(1L))
  structure(list(
    model = model,
    thresholds = thresholds,
    coefficients = coefficients,
    convergence_flag = ifelse(is.na(converged), NA_integer_,
                              ifelse(boundary, 2L, 0L)),
    convergence_diagnostics_available = any(!is.na(converged)),
    iterations = vapply(fits, `[[`, integer(1L), "iterations"),
    backend_requested = backend_requested,
    backend = paste(unique(threshold_backend), collapse = "+"),
    threshold_backend = threshold_backend,
    threshold_fallback_reason = threshold_fallback_reason,
    threshold_workers = if (model == "lpm" || backend_requested == "cuda") {
      1L
    } else {
      min(dr_workers, length(thresholds))
    },
    warm_start = model != "lpm" && backend_requested != "cuda" &&
      dr_workers == 1L && warm_start,
    preconditioned = model != "lpm" && isTRUE(precondition),
    preconditioning_method = preconditioned$method,
    preconditioning_gram_rcond = preconditioned$gram_rcond,
    preconditioning_condition_estimate = preconditioned$condition_estimate,
    preconditioning_fallback_reason = preconditioned$fallback_reason,
    warnings = fit_warnings
  ), class = c("cf_dr_fit", "cf_conditional_fit"))
}
