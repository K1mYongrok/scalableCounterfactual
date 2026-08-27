stata_exe <- Sys.getenv("STATA_EXE", unset = "")
if (!nzchar(stata_exe) || !file.exists(stata_exe)) {
  message("Skipping optional Stata parity test; set STATA_EXE to enable it.")
} else {
  library(scalableCounterfactual)

  script_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
  script_directory <- if (length(script_argument)) {
    dirname(normalizePath(sub("^--file=", "", script_argument[[1L]])))
  } else {
    getwd()
  }
  do_file <- file.path(script_directory, "stata_parity.do")
  output_file <- tempfile(fileext = ".csv")
  status <- system2(
    stata_exe,
    args = c("/e", "do", shQuote(do_file), shQuote(output_file)),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!file.exists(output_file)) {
    stop(
      "Stata parity run did not produce coefficients:\n",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }

  i <- seq_len(400L)
  x1 <- (i %% 17L - 8) / 5
  x2 <- i %% 2L
  weights <- 1 + (i %% 7L) / 10
  y <- 1 + 0.6 * x1 - 0.25 * x2 + 0.15 * sin(i / 3) +
    (0.35 + 0.08 * x2) * cos(i / 11)
  X <- cbind(x1 = x1, x2 = x2, `_cons` = 1)
  r_fit <- fit_weighted_qr(
    X,
    y,
    weights,
    taus = seq(0.1, 0.9, by = 0.01),
    solver = "onestep",
    precondition = FALSE,
    onestep_first_solver = "br"
  )
  stata_fit <- utils::read.csv(output_file, check.names = FALSE)
  stata_coefficients <- as.matrix(stata_fit[, -1L, drop = FALSE])
  rownames(stata_coefficients) <- stata_fit[[1L]]
  stata_coefficients <- stata_coefficients[c("x1", "x2", "_cons"), ]
  difference <- max(abs(r_fit$coefficients - stata_coefficients))
  tolerance <- as.numeric(Sys.getenv("STATA_PARITY_TOL", unset = "1e-6"))
  if (!is.finite(difference) || difference > tolerance) {
    stop(
      "Stata onestep parity difference ", signif(difference, 8),
      " exceeds tolerance ", tolerance,
      call. = FALSE
    )
  }
  evidence_path <- Sys.getenv("STATA_PARITY_OUTPUT", unset = "")
  if (nzchar(evidence_path)) {
    dir.create(dirname(evidence_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(data.frame(
      timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      upstream = "bmelly/Stata qrprocess 1.1.3",
      upstream_commit = "ec56830ef9c84ce54411ad59c5ce94535847d9df",
      quantile_low = 0.1,
      quantile_high = 0.9,
      quantile_step = 0.01,
      max_abs_coefficient_difference = difference,
      tolerance = tolerance,
      R_version = R.version.string,
      stringsAsFactors = FALSE
    ), evidence_path, row.names = FALSE)
  }
  message("Stata onestep parity max difference: ", signif(difference, 8))
}
