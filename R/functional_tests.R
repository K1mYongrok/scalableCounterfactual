effect_test_statistics <- function(
    bootstrap_numerator, observed_numerator, variance_draws, robust_se) {
  variance <- if (isTRUE(robust_se)) {
    apply(variance_draws, 2L, function(x) {
      (diff(stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)) / 1.34)^2
    })
  } else {
    apply(variance_draws, 2L, stats::var, na.rm = TRUE)
  }
  variance <- variance + 1e-9
  standardized_boot <- sweep(
    bootstrap_numerator^2, 2L, variance, "/"
  )
  standardized_observed <- observed_numerator^2 / variance
  bootstrap_ks <- apply(sqrt(standardized_boot), 1L, max, na.rm = TRUE)
  bootstrap_cms <- rowMeans(standardized_boot, na.rm = TRUE)
  observed_ks <- max(sqrt(standardized_observed), na.rm = TRUE)
  observed_cms <- mean(standardized_observed, na.rm = TRUE)
  c(
    ks_statistic = observed_ks,
    ks_p_value = mean(bootstrap_ks > observed_ks, na.rm = TRUE),
    cms_statistic = observed_cms,
    cms_p_value = mean(bootstrap_cms > observed_cms, na.rm = TRUE)
  )
}

functional_test_row <- function(
    effect, test, null_value, bootstrap_numerator, observed_numerator,
    variance_draws, robust_se, quantile_min, quantile_max) {
  statistics <- effect_test_statistics(
    bootstrap_numerator, observed_numerator, variance_draws, robust_se
  )
  data.frame(
    effect = effect,
    test = test,
    null_value = null_value,
    quantile_min = quantile_min,
    quantile_max = quantile_max,
    quantiles_tested = ncol(bootstrap_numerator),
    bootstrap_reps = nrow(bootstrap_numerator),
    ks_statistic = unname(statistics[["ks_statistic"]]),
    ks_p_value = unname(statistics[["ks_p_value"]]),
    cms_statistic = unname(statistics[["cms_statistic"]]),
    cms_p_value = unname(statistics[["cms_p_value"]]),
    stringsAsFactors = FALSE
  )
}

#' Bootstrap KS and CMS tests for decomposition-effect curves
#'
#' Implements the studentized Kolmogorov--Smirnov (KS) and Cramer--von Mises
#' (CMS) effect-curve tests used by `Counterfactual` 1.2. Tests are available
#' for a zero or user-specified constant effect, equality to the median effect,
#' and one-sided nonnegative or nonpositive effects. The function requires a
#' fitted object with bootstrap replications.
#'
#' @param object A `cfdecomp` object estimated with `bootstrap_reps > 0`.
#' @param constants Optional nonzero constant-effect null values.
#' @param quantile_range Length-two range of reported quantiles used by the
#'   tests.
#' @return A data frame containing KS and CMS statistics and bootstrap
#'   p-values for each decomposition effect and null hypothesis.
#' @export
functional_effect_tests <- function(
    object, constants = numeric(), quantile_range = c(0.1, 0.9)) {
  if (!inherits(object, "cfdecomp")) {
    stop("object must be a cfdecomp fit", call. = FALSE)
  }
  if (is.null(object$bootstrap) || !length(object$bootstrap$effects)) {
    stop("functional tests require bootstrap replications", call. = FALSE)
  }
  if (nrow(object$bootstrap$effects[[1L]]) < 2L) {
    stop("functional tests require at least two bootstrap replications",
         call. = FALSE)
  }
  if (!is.numeric(constants) || any(!is.finite(constants))) {
    stop("constants must contain finite numeric values", call. = FALSE)
  }
  constants <- sort(unique(constants[constants != 0]))
  if (!is.numeric(quantile_range) || length(quantile_range) != 2L ||
      any(!is.finite(quantile_range)) || quantile_range[[1L]] < 0 ||
      quantile_range[[2L]] > 1 || quantile_range[[1L]] >= quantile_range[[2L]]) {
    stop("quantile_range must be an increasing length-two vector in [0, 1]",
         call. = FALSE)
  }
  quantiles <- object$point$quantiles
  selected_requested <- which(
    quantiles >= quantile_range[[1L]] & quantiles <= quantile_range[[2L]]
  )
  if (length(selected_requested) < 2L) {
    stop("quantile_range must retain at least two reported quantiles",
         call. = FALSE)
  }
  rows <- list()
  append_test <- function(...) {
    rows[[length(rows) + 1L]] <<- functional_test_row(...)
  }
  for (effect in c("structure", "composition", "total")) {
    draws_all <- object$bootstrap$effects[[effect]]
    point_all <- as.numeric(object$point$effects[effect, ])
    selected <- selected_requested[
      is.finite(point_all[selected_requested]) &
        vapply(selected_requested, function(index) {
          all(is.finite(draws_all[, index]))
        }, logical(1L))
    ]
    if (length(selected) < 2L) {
      warning(
        "Functional tests skipped for ", effect,
        ": fewer than two identified quantiles have complete bootstrap draws",
        call. = FALSE
      )
      next
    }
    median_index <- selected[[which.min(abs(quantiles[selected] - 0.5))]]
    if (abs(quantiles[[median_index]] - 0.5) > sqrt(.Machine$double.eps)) {
      warning(
        "0.5 is not identified for ", effect,
        "; the nearest identified reported quantile is used as median",
        call. = FALSE
      )
    }
    draws <- draws_all[, selected, drop = FALSE]
    point <- point_all[selected]
    centered <- sweep(draws, 2L, point, "-")
    append_test(
      effect, "zero_effect", 0, centered, point, draws,
      object$control$robust_se, min(quantiles[selected]), max(quantiles[selected])
    )
    for (constant in constants) {
      append_test(
        effect, "constant_effect", constant, centered, point - constant, draws,
        object$control$robust_se, min(quantiles[selected]),
        max(quantiles[selected])
      )
    }

    nonmedian <- setdiff(selected, median_index)
    if (length(nonmedian)) {
      median_draws <- draws_all[, median_index]
      transformed_draws <- sweep(
        draws_all[, nonmedian, drop = FALSE], 1L, median_draws, "-"
      )
      transformed_point <- point_all[nonmedian] - point_all[[median_index]]
      transformed_centered <- sweep(
        transformed_draws, 2L, transformed_point, "-"
      )
      append_test(
        effect, "constant_across_quantiles", point_all[[median_index]],
        transformed_centered, transformed_point, transformed_draws,
        object$control$robust_se, min(quantiles[nonmedian]),
        max(quantiles[nonmedian])
      )
    }

    negative_boot <- centered * (centered <= 0)
    positive_boot <- centered * (centered >= 0)
    append_test(
      effect, "nonnegative_effect", 0, negative_boot,
      point * (point <= 0), draws, object$control$robust_se,
      min(quantiles[selected]), max(quantiles[selected])
    )
    append_test(
      effect, "nonpositive_effect", 0, positive_boot,
      point * (point >= 0), draws, object$control$robust_se,
      min(quantiles[selected]), max(quantiles[selected])
    )
  }
  if (!length(rows)) return(data.frame())
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}
