#' @export
print.cfdecomp <- function(x, ...) {
  cat("Scalable counterfactual decomposition\n")
  cat("  model:", x$metadata$model, "\n")
  if (!is.na(x$metadata$solver)) cat("  QR solver:", x$metadata$solver, "\n")
  cat("  observations:", format(x$metadata$analysis_rows, big.mark = ","), "\n")
  cat("  point elapsed seconds:", round(x$metadata$point_elapsed_seconds, 3), "\n")
  cat("  bootstrap replications:", x$metadata$bootstrap_reps, "\n")
  if (x$metadata$bootstrap_reps > 0L) {
    cat(
      "  bootstrap retries used:",
      x$metadata$bootstrap_retries_used %||% 0L, "\n"
    )
  }
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.cfdecomp <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$results
}

#' Summarize a counterfactual decomposition
#'
#' @param object A `cfdecomp` object.
#' @param effects Effects to retain: any of `structure`, `composition`, and
#'   `total`.
#' @param quantiles Optional reported quantiles to retain. The nearest stored
#'   quantile is used for each requested value.
#' @param ... Unused.
#' @return An object of class `summary.cfdecomp`.
#' @export
summary.cfdecomp <- function(
    object,
    effects = c("structure", "composition", "total"),
    quantiles = NULL,
    ...) {
  if (!inherits(object, "cfdecomp")) {
    stop("object must be a cfdecomp fit", call. = FALSE)
  }
  effects <- unique(as.character(effects))
  invalid <- setdiff(effects, c("structure", "composition", "total"))
  if (!length(effects) || length(invalid)) {
    stop(
      "effects must contain structure, composition, or total",
      call. = FALSE
    )
  }
  table <- object$results[object$results$effect %in% effects, , drop = FALSE]
  if (!is.null(quantiles)) {
    quantiles <- as.numeric(quantiles)
    if (!length(quantiles) || any(!is.finite(quantiles)) ||
        any(quantiles <= 0 | quantiles >= 1)) {
      stop("quantiles must lie strictly between 0 and 1", call. = FALSE)
    }
    available <- sort(unique(table$quantile))
    selected <- unique(vapply(quantiles, function(probability) {
      available[[which.min(abs(available - probability))]]
    }, numeric(1L)))
    table <- table[table$quantile %in% selected, , drop = FALSE]
  }
  table <- table[order(match(table$effect, effects), table$quantile), , drop = FALSE]
  crossing <- object$point$crossing_diagnostics
  finite_max <- function(values) {
    values <- values[is.finite(values)]
    if (length(values)) max(values) else NA_real_
  }
  diagnostics <- list(
    maximum_identity_residual = finite_max(abs(object$point$identity_residual)),
    maximum_crossing_row_share = if (!is.null(crossing) && nrow(crossing)) {
      finite_max(crossing$crossing_row_share)
    } else {
      NA_real_
    },
    marginalization_methods = paste(
      unique(object$point$marginal_diagnostics$method), collapse = ","
    ),
    bootstrap_retries = object$metadata$bootstrap_retries_used %||% 0L,
    warning_count = length(object$point$warnings),
    warnings = paste(object$point$warnings, collapse = " | ")
  )
  structure(list(
    call = object$call,
    model = object$metadata$model,
    solver = object$metadata$solver,
    observations = object$metadata$analysis_rows,
    group0_rows = object$metadata$group0_rows,
    group1_rows = object$metadata$group1_rows,
    bootstrap_reps = object$metadata$bootstrap_reps,
    alpha = object$control$alpha,
    effects = table,
    diagnostics = diagnostics,
    functional_tests = object$functional_tests %||% data.frame()
  ), class = "summary.cfdecomp")
}

#' @export
print.summary.cfdecomp <- function(x, ...) {
  cat("Scalable counterfactual decomposition summary\n")
  cat("  model:", x$model, "\n")
  if (!is.na(x$solver)) cat("  QR solver:", x$solver, "\n")
  cat("  observations:", format(x$observations, big.mark = ","), "\n")
  cat("  groups:", format(x$group0_rows, big.mark = ","), "reference;",
      format(x$group1_rows, big.mark = ","), "comparison\n")
  cat("  bootstrap replications:", x$bootstrap_reps, "\n")
  cat("  marginalization:", x$diagnostics$marginalization_methods, "\n")
  cat("  max decomposition identity residual:",
      format(x$diagnostics$maximum_identity_residual, scientific = TRUE), "\n")
  print(x$effects, row.names = FALSE)
  if (nrow(x$functional_tests)) {
    cat("\nFunctional tests are available in $functional_tests.\n")
  }
  invisible(x)
}

#' Plot counterfactual decomposition effects
#'
#' Draws one base-R panel per requested effect with a zero reference line and
#' optional pointwise or uniform confidence band.
#'
#' @param x A `cfdecomp` object.
#' @param effects Effects to plot.
#' @param interval `pointwise`, `uniform`, or `none`.
#' @param col Line and band colors, recycled across effects.
#' @param lwd Line width.
#' @param pch Point symbol.
#' @param xlab,ylab Axis labels.
#' @param main Optional overall title. Effect names remain panel titles.
#' @param ... Additional arguments passed to [graphics::plot.default()].
#' @return Invisibly returns the plotted effect data.
#' @export
plot.cfdecomp <- function(
    x,
    effects = c("structure", "composition", "total"),
    interval = c("pointwise", "uniform", "none"),
    col = c("#1B7837", "#D95F02", "#2166AC"),
    lwd = 2,
    pch = 16,
    xlab = "Unconditional quantile",
    ylab = "Effect",
    main = NULL,
    ...) {
  if (!inherits(x, "cfdecomp")) stop("x must be a cfdecomp fit", call. = FALSE)
  interval <- match.arg(interval)
  effects <- unique(as.character(effects))
  invalid <- setdiff(effects, c("structure", "composition", "total"))
  if (!length(effects) || length(invalid)) {
    stop("invalid effects requested", call. = FALSE)
  }
  if (!length(col) || anyNA(col)) stop("col must contain at least one color",
                                      call. = FALSE)
  colors <- rep(col, length.out = length(effects))
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  if (length(effects) > 1L) {
    graphics::par(mfrow = grDevices::n2mfrow(length(effects)))
  }
  if (!is.null(main)) graphics::par(oma = c(0, 0, 2, 0))
  plotted <- vector("list", length(effects))
  names(plotted) <- effects
  plot_arguments <- list(...)
  supplied_ylim <- plot_arguments$ylim
  plot_arguments$type <- NULL
  plot_arguments$ylim <- NULL
  for (i in seq_along(effects)) {
    effect <- effects[[i]]
    selected <- x$results[x$results$effect == effect, , drop = FALSE]
    selected <- selected[order(selected$quantile), , drop = FALSE]
    lower_name <- paste0(interval, "_lower")
    upper_name <- paste0(interval, "_upper")
    interval_rows <- rep(FALSE, nrow(selected))
    if (interval != "none") {
      interval_rows <- is.finite(selected[[lower_name]]) &
        is.finite(selected[[upper_name]])
    }
    has_interval <- any(interval_rows)
    y_values <- selected$estimate
    if (has_interval) {
      y_values <- c(
        y_values,
        selected[[lower_name]][interval_rows],
        selected[[upper_name]][interval_rows]
      )
    }
    panel_ylim <- supplied_ylim %||% range(c(0, y_values), finite = TRUE)
    do.call(graphics::plot.default, c(list(
      x = selected$quantile,
      y = selected$estimate,
      type = "n",
      xlab = xlab,
      ylab = ylab,
      main = effect,
      ylim = panel_ylim
    ), plot_arguments))
    graphics::abline(h = 0, lty = 2, col = "grey45")
    if (has_interval) {
      indices <- which(interval_rows)
      interval_runs <- split(indices, cumsum(c(TRUE, diff(indices) > 1L)))
      for (run in interval_runs) {
        if (length(run) == 1L) {
          graphics::segments(
            selected$quantile[run], selected[[lower_name]][run],
            selected$quantile[run], selected[[upper_name]][run],
            col = grDevices::adjustcolor(colors[[i]], alpha.f = 0.35)
          )
        } else {
          graphics::polygon(
            c(selected$quantile[run], rev(selected$quantile[run])),
            c(
              selected[[lower_name]][run],
              rev(selected[[upper_name]][run])
            ),
            border = NA,
            col = grDevices::adjustcolor(colors[[i]], alpha.f = 0.18)
          )
        }
      }
    }
    graphics::lines(
      selected$quantile, selected$estimate,
      col = colors[[i]], lwd = lwd
    )
    graphics::points(
      selected$quantile, selected$estimate,
      col = colors[[i]], pch = pch
    )
    plotted[[i]] <- selected
  }
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, line = 0.5)
  invisible(plotted)
}

#' Write decomposition outputs
#'
#' @param object A `cfdecomp` object.
#' @param output_dir Destination directory.
#' @return Invisibly returns `output_dir`.
#' @export
write_cf_outputs <- function(object, output_dir) {
  if (!inherits(object, "cfdecomp")) stop("object must be cfdecomp", call. = FALSE)
  managed_files <- c(
    "decomposition.csv", "distribution_diagnostics.csv",
    "point_resources.csv", "marginalization_diagnostics.csv",
    "run_metadata.csv", "fit.rds", "bootstrap_resources.csv",
    "bootstrap_failures.csv", "functional_effect_tests.csv",
    "quantile_crossing_diagnostics.csv"
  )
  atomic_write_output_files(
    output_dir, managed_files,
    required_files = c(
      "decomposition.csv", "distribution_diagnostics.csv",
      "point_resources.csv", "marginalization_diagnostics.csv",
      "run_metadata.csv", "fit.rds"
    ),
    writer = function(output_dir) {
  optional_outputs <- file.path(output_dir, c(
    "bootstrap_resources.csv", "bootstrap_failures.csv",
    "functional_effect_tests.csv",
    "quantile_crossing_diagnostics.csv"
  ))
  existing_optional <- optional_outputs[file.exists(optional_outputs)]
  if (length(existing_optional)) unlink(existing_optional)
  data.table::fwrite(object$results, file.path(output_dir, "decomposition.csv"))
  diagnostics <- data.frame(
    quantile = object$point$quantiles,
    t(object$point$diagnostics),
    check.names = FALSE
  )
  data.table::fwrite(
    diagnostics, file.path(output_dir, "distribution_diagnostics.csv")
  )
  data.table::fwrite(
    object$point$resources, file.path(output_dir, "point_resources.csv")
  )
  if (!is.null(object$point$marginal_diagnostics)) {
    data.table::fwrite(
      object$point$marginal_diagnostics,
      file.path(output_dir, "marginalization_diagnostics.csv")
    )
  }
  if (!is.null(object$point$crossing_diagnostics) &&
      nrow(object$point$crossing_diagnostics)) {
    data.table::fwrite(
      object$point$crossing_diagnostics,
      file.path(output_dir, "quantile_crossing_diagnostics.csv")
    )
  }
  metadata <- data.frame(
    item = names(object$metadata),
    value = vapply(object$metadata, function(x) paste(x, collapse = ","), character(1L)),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(metadata, file.path(output_dir, "run_metadata.csv"))
  if (!is.null(object$bootstrap)) {
    data.table::fwrite(
      object$bootstrap$resources,
      file.path(output_dir, "bootstrap_resources.csv")
    )
    if (!is.null(object$bootstrap$failures) &&
        nrow(object$bootstrap$failures)) {
      data.table::fwrite(
        object$bootstrap$failures,
        file.path(output_dir, "bootstrap_failures.csv")
      )
    }
  }
  if (!is.null(object$functional_tests) && nrow(object$functional_tests)) {
    data.table::fwrite(
      object$functional_tests,
      file.path(output_dir, "functional_effect_tests.csv")
    )
  }
  saveRDS(object, file.path(output_dir, "fit.rds"), compress = "xz")
    }
  )
  invisible(output_dir)
}
