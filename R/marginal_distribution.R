# Marginalization implements the counterfactual integration step; optional
# row-wise rearrangement follows doi:10.3982/ECTA7880.
conditional_matmul <- function(X, coefficients, control) {
  if (identical(control$gpu_backend, "cuda")) {
    return(gpu_matmul(X, coefficients, control))
  }
  X %*% coefficients
}

predict_conditional_draws <- function(fit, X, control) {
  if (inherits(fit, "cf_qr_fit") || inherits(fit, "cf_cqr_fit")) {
    shift <- if (inherits(fit, "cf_qr_fit") &&
                 isTRUE(control$legacy_qr_shift)) control$trimming else 0
    draws <- conditional_matmul(X, fit$coefficients, control) + shift
    if (control$quantile_noncrossing == "rearrange") {
      draws <- rearrange_quantile_rows(draws)
    }
    return(draws)
  }
  if (inherits(fit, "cf_loc_fit")) {
    location <- drop(conditional_matmul(X, fit$coefficients, control))
    return(outer(location, rep(1, length(fit$residual_quantiles))) +
      rep(fit$residual_quantiles, each = length(location)))
  }
  if (inherits(fit, "cf_locsca_fit")) {
    location <- drop(conditional_matmul(
      X, fit$location_coefficients, control
    ))
    scale <- sqrt(exp(drop(conditional_matmul(
      X, fit$scale_coefficients, control
    ))))
    return(outer(location, rep(1, length(fit$residual_quantiles))) +
      outer(scale, fit$residual_quantiles))
  }
  stop("draw prediction is not defined for this model", call. = FALSE)
}

conditional_draw_count <- function(fit) {
  if (inherits(fit, "cf_qr_fit") || inherits(fit, "cf_cqr_fit")) {
    return(ncol(fit$coefficients))
  }
  if (inherits(fit, "cf_loc_fit") || inherits(fit, "cf_locsca_fit")) {
    return(length(fit$residual_quantiles))
  }
  stop("draw count is not defined for this model", call. = FALSE)
}

rearrange_quantile_rows <- function(draws) {
  draws <- as.matrix(draws)
  if (nrow(draws) == 1L) {
    return(matrix(sort(draws[1L, ]), nrow = 1L,
                  dimnames = list(rownames(draws), colnames(draws))))
  }
  rearranged <- t(apply(draws, 1L, sort, method = "quick"))
  dimnames(rearranged) <- dimnames(draws)
  rearranged
}

qr_crossing_diagnostics_single <- function(
    fit, X, weights, distribution, chunk_rows, tolerance,
    noncrossing = "none", stage = "raw", control) {
  pairs_per_row <- max(0L, length(fit$taus) - 1L)
  if (!pairs_per_row) {
    return(data.frame(
      distribution = distribution, rows = nrow(X), quantile_pairs = 0,
      crossing_rows = 0, crossing_row_share = 0,
      weighted_crossing_row_share = 0, crossing_pairs = 0,
      crossing_pair_share = 0, max_crossing = 0,
      mean_crossing = 0, tolerance = tolerance,
      correction = noncrossing, stage = stage,
      stringsAsFactors = FALSE
    ))
  }
  crossing_rows <- 0
  weighted_crossing_rows <- 0
  crossing_pairs <- 0
  crossing_sum <- 0
  crossing_max <- 0
  chunks <- draw_chunk_rows(nrow(X), chunk_rows)
  for (rows in chunks) {
    predictions <- conditional_matmul(
      X[rows, , drop = FALSE], fit$coefficients, control
    )
    if (noncrossing == "rearrange") {
      predictions <- rearrange_quantile_rows(predictions)
    }
    differences <- predictions[, -1L, drop = FALSE] -
      predictions[, -ncol(predictions), drop = FALSE]
    crossed <- differences < -tolerance
    row_crossed <- rowSums(crossed) > 0L
    magnitudes <- -differences[crossed]
    crossing_rows <- crossing_rows + sum(row_crossed)
    weighted_crossing_rows <- weighted_crossing_rows +
      sum(weights[rows] * row_crossed)
    crossing_pairs <- crossing_pairs + length(magnitudes)
    if (length(magnitudes)) {
      crossing_sum <- crossing_sum + sum(magnitudes)
      crossing_max <- max(crossing_max, magnitudes)
    }
  }
  data.frame(
    distribution = distribution,
    rows = nrow(X),
    quantile_pairs = as.numeric(nrow(X)) * pairs_per_row,
    crossing_rows = crossing_rows,
    crossing_row_share = crossing_rows / nrow(X),
    weighted_crossing_row_share = weighted_crossing_rows / sum(weights),
    crossing_pairs = crossing_pairs,
    crossing_pair_share = crossing_pairs /
      (as.numeric(nrow(X)) * pairs_per_row),
    max_crossing = crossing_max,
    mean_crossing = if (crossing_pairs) crossing_sum / crossing_pairs else 0,
    tolerance = tolerance,
    correction = noncrossing,
    stage = stage,
    stringsAsFactors = FALSE
  )
}

evaluate_qr_crossing_diagnostics <- function(prepared, fit0, fit1, control) {
  designs <- list(
    list(fit0, prepared$X0, prepared$w0, "reference_model_X0"),
    list(fit0, prepared$X1, prepared$w1, "reference_model_X1_counterfactual"),
    list(fit1, prepared$X1, prepared$w1, "comparison_model_X1")
  )
  raw <- lapply(designs, function(item) {
    qr_crossing_diagnostics_single(
      item[[1L]], item[[2L]], item[[3L]], item[[4L]],
      control$marginal_chunk_rows, control$crossing_tolerance,
      noncrossing = "none", stage = "raw", control = control
    )
  })
  if (control$quantile_noncrossing == "none") return(do.call(rbind, raw))
  corrected <- lapply(designs, function(item) {
    qr_crossing_diagnostics_single(
      item[[1L]], item[[2L]], item[[3L]], item[[4L]],
      control$marginal_chunk_rows, control$crossing_tolerance,
      noncrossing = control$quantile_noncrossing, stage = "corrected",
      control = control
    )
  })
  do.call(rbind, c(raw, corrected))
}

resolve_marginal_method <- function(fit, X, control) {
  requested <- control$marginal_method
  draws <- conditional_draw_count(fit)
  estimated_mb <- as.numeric(nrow(X)) * draws * 8 / 1024^2
  if (isTRUE(control$legacy_weighted_quantile)) {
    if (requested == "chunked") {
      stop(
        "legacy_weighted_quantile is incompatible with chunked marginalization",
        call. = FALSE
      )
    }
    return(list(
      requested = requested,
      resolved = "matrix",
      estimated_matrix_mb = estimated_mb,
      draws_per_row = draws
    ))
  }
  resolved <- if (requested == "auto") {
    if (estimated_mb <= control$marginal_matrix_max_mb) "matrix" else "chunked"
  } else {
    requested
  }
  list(
    requested = requested,
    resolved = resolved,
    estimated_matrix_mb = estimated_mb,
    draws_per_row = draws
  )
}

draw_chunk_rows <- function(n, chunk_rows) {
  starts <- seq.int(1L, n, by = chunk_rows)
  Map(seq.int, starts, pmin(n, starts + chunk_rows - 1L))
}

draw_bin_index <- function(values, lower, upper, bins) {
  scaled <- (values - lower) / (upper - lower)
  pmax(1L, pmin(bins, as.integer(floor(scaled * bins)) + 1L))
}

chunked_draw_range <- function(fit, X, control, chunks) {
  lower <- Inf
  upper <- -Inf
  for (rows in chunks) {
    draws <- predict_conditional_draws(fit, X[rows, , drop = FALSE], control)
    lower <- min(lower, draws)
    upper <- max(upper, draws)
  }
  c(lower = lower, upper = upper)
}

chunked_draw_histogram <- function(
    fit, X, weights, control, chunks, lower, upper, bins) {
  weight_histogram <- numeric(bins)
  count_histogram <- numeric(bins)
  for (rows in chunks) {
    draws <- predict_conditional_draws(fit, X[rows, , drop = FALSE], control)
    values <- as.vector(draws)
    bin <- draw_bin_index(values, lower, upper, bins)
    draw_weights <- rep(weights[rows], times = ncol(draws))
    aggregated <- rowsum(draw_weights, bin, reorder = FALSE)
    occupied <- as.integer(rownames(aggregated))
    weight_histogram[occupied] <- weight_histogram[occupied] +
      as.numeric(aggregated[, 1L])
    count_histogram <- count_histogram + tabulate(bin, nbins = bins)
  }
  list(weights = weight_histogram, counts = count_histogram)
}

collect_target_bin_draws <- function(
    fit, X, weights, control, chunks, lower, upper, bins, target_bins) {
  pieces <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    rows <- chunks[[i]]
    draws <- predict_conditional_draws(fit, X[rows, , drop = FALSE], control)
    values <- as.vector(draws)
    bin <- draw_bin_index(values, lower, upper, bins)
    keep <- bin %in% target_bins
    if (!any(keep)) next
    draw_weights <- rep(weights[rows], times = ncol(draws))
    piece <- data.table::data.table(
      bin = bin[keep], value = values[keep], weight = draw_weights[keep]
    )
    pieces[[i]] <- piece[, list(weight = sum(weight)), by = c("bin", "value")]
  }
  pieces <- pieces[vapply(pieces, function(x) !is.null(x) && nrow(x), logical(1L))]
  if (!length(pieces)) {
    stop("chunked marginalization did not retain target-bin draws", call. = FALSE)
  }
  candidates <- data.table::rbindlist(pieces)
  candidates <- candidates[
    , list(weight = sum(weight)), by = c("bin", "value")
  ]
  data.table::setorder(candidates, bin, value)
  candidates
}

chunked_weighted_draw_quantile <- function(fit, X, weights, probs, control) {
  n <- nrow(X)
  draws_per_row <- conditional_draw_count(fit)
  chunks <- draw_chunk_rows(n, control$marginal_chunk_rows)
  range <- chunked_draw_range(fit, X, control, chunks)
  if (!all(is.finite(range))) {
    stop("conditional draws contain non-finite values", call. = FALSE)
  }
  if (range[["lower"]] == range[["upper"]]) {
    result <- rep(range[["lower"]], length(probs))
    attr(result, "marginal_diagnostics") <- list(
      method = "chunked", passes = 1L, histogram_bins = 0L,
      candidate_draws = 1, estimated_matrix_mb = n * draws_per_row * 8 / 1024^2
    )
    return(result)
  }

  total_positions <- as.numeric(n) * draws_per_row
  order <- 1 + (total_positions - 1) * probs
  low_position <- pmax(floor(order), 1)
  high_position <- pmin(low_position + 1, total_positions)
  interpolation <- order %% 1
  target_positions <- c(low_position, high_position)
  total_weight <- sum(weights) * draws_per_row
  target_weights <- target_positions / total_positions * total_weight

  bins <- control$marginal_histogram_bins
  max_bins <- 4194304L
  repeat {
    histogram <- chunked_draw_histogram(
      fit, X, weights, control, chunks,
      range[["lower"]], range[["upper"]], bins
    )
    cumulative <- cumsum(histogram$weights)
    cumulative[[length(cumulative)]] <- total_weight
    target_bins <- vapply(target_weights, function(target) {
      which(cumulative >= target)[[1L]]
    }, integer(1L))
    candidate_draws <- sum(histogram$counts[unique(target_bins)])
    if (candidate_draws <= control$marginal_candidate_max || bins >= max_bins) {
      break
    }
    bins <- min(max_bins, bins * 4L)
  }
  if (candidate_draws > control$marginal_candidate_max) {
    stop(
      "chunked marginalization target bins contain ", candidate_draws,
      " draws; increase marginal_candidate_max or marginal_histogram_bins",
      call. = FALSE
    )
  }

  candidates <- collect_target_bin_draws(
    fit, X, weights, control, chunks,
    range[["lower"]], range[["upper"]], bins, unique(target_bins)
  )
  selected <- vapply(seq_along(target_weights), function(i) {
    target_bin <- target_bins[[i]]
    before <- if (target_bin == 1L) 0 else cumulative[[target_bin - 1L]]
    local <- candidates[candidates$bin == target_bin]
    local_cumulative <- before + cumsum(local$weight)
    local_cumulative[[length(local_cumulative)]] <- cumulative[[target_bin]]
    hit <- which(local_cumulative >= target_weights[[i]])[[1L]]
    local$value[[hit]]
  }, numeric(1L))
  k <- length(probs)
  result <- (1 - interpolation) * selected[seq_len(k)] +
    interpolation * selected[k + seq_len(k)]
  attr(result, "marginal_diagnostics") <- list(
    method = "chunked",
    passes = 2L + if (bins == control$marginal_histogram_bins) 1L else {
      1L + round(log(bins / control$marginal_histogram_bins, base = 4))
    },
    histogram_bins = bins,
    candidate_draws = candidate_draws,
    estimated_matrix_mb = total_positions * 8 / 1024^2
  )
  result
}

inverse_step_cdf <- function(thresholds, cdf, probs) {
  ordering <- order(thresholds)
  thresholds <- thresholds[ordering]
  cdf <- cdf[ordering]
  cdf <- cummax(pmin(1, pmax(0, cdf)))
  vapply(probs, function(probability) {
    hit <- which(cdf >= probability)[1L]
    if (is.na(hit)) thresholds[[length(thresholds)]] else thresholds[[hit]]
  }, numeric(1L))
}

marginal_quantiles <- function(fit, X, weights, probs, control) {
  if (inherits(fit, "cf_cox_fit")) {
    return(cox_marginal_quantiles(
      fit, X, weights, probs, control$cox_boundary
    ))
  }
  if (inherits(fit, "cf_dr_fit")) {
    if (identical(control$gpu_backend, "cuda")) {
      marginal_cdf <- gpu_dr_marginal_cdf(
        X, fit$coefficients, weights, fit$model, control
      )
    } else {
      linear_predictor <- X %*% fit$coefficients
      probabilities <- switch(
        fit$model,
        logit = stats::plogis(linear_predictor),
        probit = stats::pnorm(linear_predictor),
        cloglog = {
          values <- -expm1(-exp(pmin(700, as.numeric(linear_predictor))))
          matrix(values, nrow = nrow(linear_predictor),
                 ncol = ncol(linear_predictor))
        },
        lpm = linear_predictor
      )
      marginal_cdf <- colSums(probabilities * weights) / sum(weights)
    }
    result <- inverse_step_cdf(fit$thresholds, marginal_cdf, probs)
    attr(result, "marginal_diagnostics") <- list(
      method = if (identical(control$gpu_backend, "cuda")) {
        "cdf_cuda"
      } else {
        "cdf"
      }, passes = 1L, histogram_bins = NA_integer_,
      candidate_draws = NA_real_, estimated_matrix_mb =
        as.numeric(nrow(X)) * ncol(fit$coefficients) * 8 / 1024^2
    )
    return(result)
  }
  if (identical(control$gpu_backend, "cuda") &&
      (inherits(fit, "cf_qr_fit") || inherits(fit, "cf_cqr_fit")) &&
      !isTRUE(control$legacy_weighted_quantile) &&
      !identical(control$marginal_method, "chunked")) {
    shift <- if (inherits(fit, "cf_qr_fit") &&
                 isTRUE(control$legacy_qr_shift)) control$trimming else 0
    result <- gpu_qr_marginal_quantiles(
      X, fit$coefficients, weights, probs, shift, control
    )
    attr(result, "marginal_diagnostics") <- list(
      method = "cuda_weighted_quantile", passes = 1L,
      histogram_bins = NA_integer_, candidate_draws =
        as.numeric(nrow(X)) * ncol(fit$coefficients),
      estimated_matrix_mb =
        as.numeric(nrow(X)) * ncol(fit$coefficients) * 8 / 1024^2
    )
    return(result)
  }
  method <- resolve_marginal_method(fit, X, control)
  if (method$resolved == "chunked") {
    return(chunked_weighted_draw_quantile(fit, X, weights, probs, control))
  }
  draws <- predict_conditional_draws(fit, X, control)
  result <- weighted_quantile(
    as.vector(draws),
    rep(weights, times = ncol(draws)),
    probs,
    legacy = control$legacy_weighted_quantile
  )
  rm(draws)
  invisible(gc())
  attr(result, "marginal_diagnostics") <- list(
    method = if (identical(control$gpu_backend, "cuda")) {
      "matrix_cuda_matmul_cpu_quantile"
    } else {
      "matrix"
    }, passes = 1L, histogram_bins = NA_integer_,
    candidate_draws = as.numeric(nrow(X)) * method$draws_per_row,
    estimated_matrix_mb = method$estimated_matrix_mb
  )
  result
}
