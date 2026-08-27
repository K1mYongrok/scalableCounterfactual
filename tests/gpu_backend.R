library(scalableCounterfactual)

gpu_python <- Sys.getenv("SCALABLECF_GPU_PYTHON", unset = NA_character_)
gpu_path <- Sys.getenv("SCALABLECF_GPU_PATH", unset = NA_character_)
if (is.na(gpu_python) || is.na(gpu_path)) {
  message("Skipping CUDA tests: set SCALABLECF_GPU_PYTHON and SCALABLECF_GPU_PATH")
} else {
  status <- gpu_backend_status(gpu_python, gpu_path)
  stopifnot(isTRUE(status$available))
  path_after_first <- Sys.getenv("PATH")
  pythonpath_after_first <- Sys.getenv("PYTHONPATH")
  repeated_status <- replicate(
    4L, gpu_backend_status(gpu_python, gpu_path), simplify = FALSE
  )
  stopifnot(
    all(vapply(repeated_status, `[[`, logical(1L), "available")),
    identical(Sys.getenv("PATH"), path_after_first),
    identical(Sys.getenv("PYTHONPATH"), pythonpath_after_first)
  )

  set.seed(42)
  n <- 1000L
  example <- data.frame(
    y = 1 + 0.3 * rep(0:1, each = n / 2) +
      0.7 * rnorm(n) + rep(c(0.8, 1), each = n / 2) * rnorm(n),
    x = rnorm(n),
    group = rep(0:1, each = n / 2),
    weight = runif(n, 0.5, 2)
  )
  common <- list(
    nreg = 9L, reported_quantiles = c(0.25, 0.5, 0.75),
    legacy_qr_shift = FALSE, crossing_diagnostics = TRUE,
    marginal_method = "matrix"
  )
  cpu_control <- do.call(cf_control, c(common, list(gpu_backend = "cpu")))
  cuda_control <- do.call(cf_control, c(common, list(
    gpu_backend = "cuda", gpu_precision = "float64",
    gpu_python = gpu_python, gpu_python_path = gpu_path
  )))
  cpu <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "qr", solver = "fn",
    control = cpu_control
  )
  cuda <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "qr", solver = "fn",
    control = cuda_control
  )
  stopifnot(max(abs(cpu$point$effects - cuda$point$effects)) < 1e-10)
  stopifnot(max(abs(
    cpu$point$crossing_diagnostics$max_crossing -
      cuda$point$crossing_diagnostics$max_crossing
  )) < 1e-10)
  cuda_two_fit_workers <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "qr", solver = "fn",
    control = cuda_control, point_workers = 2L
  )
  stopifnot(
    max(abs(cpu$point$effects -
      cuda_two_fit_workers$point$effects)) < 1e-10,
    identical(cuda_two_fit_workers$metadata$fit_device, "cpu"),
    identical(cuda_two_fit_workers$metadata$prediction_device, "cuda"),
    identical(cuda_two_fit_workers$metadata$marginalization_device, "cuda")
  )

  admm_control <- cf_control(
    nreg = 3L, trimming = 0.2, reported_quantiles = 0.5,
    crossing_diagnostics = FALSE,
    gpu_backend = "cuda", gpu_precision = "float64",
    gpu_python = gpu_python, gpu_python_path = gpu_path,
    gpu_qr_maxit = 10L, gpu_qr_tolerance = 1e-12,
    gpu_qr_allow_nonconvergence = TRUE
  )
  admm_design <- stats::model.matrix(~ x, example[seq_len(200L), ])
  admm_fit <- suppressWarnings(fit_weighted_qr(
    admm_design, example$y[seq_len(200L)],
    example$weight[seq_len(200L)], 0.5, solver = "cuda_admm",
    precondition = TRUE, gpu_control = admm_control
  ))
  final_admm_convergence <-
    admm_fit$cuda_admm_primal_residual <= admm_control$gpu_qr_tolerance &
    admm_fit$cuda_admm_dual_residual <= admm_control$gpu_qr_tolerance
  stopifnot(
    isTRUE(admm_fit$preconditioned),
    admm_fit$preconditioning_method != "none",
    isTRUE(admm_fit$convergence_diagnostics_available),
    all(admm_fit$convergence_flag %in% c(0L, 1L)),
    identical(admm_fit$cuda_admm_converged, final_admm_convergence)
  )

  # Ties and unequal weights exercise the Hmisc-compatible rank rule used by
  # GPU-resident QR marginalization.
  tied <- example
  tied$y <- round(tied$y, 1)
  tied$weight <- rep(c(0.5, 1, 2, 3), length.out = nrow(tied))
  tied_cpu <- counterfactual_decompose(
    y ~ x, tied, "group", "weight", model = "qr", solver = "fn",
    control = cpu_control
  )
  tied_cuda <- counterfactual_decompose(
    y ~ x, tied, "group", "weight", model = "qr", solver = "fn",
    control = cuda_control
  )
  stopifnot(max(abs(
    tied_cpu$point$effects - tied_cuda$point$effects
  )) < 1e-10)

  cpu_dr <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "logit",
    control = do.call(cf_control, utils::modifyList(common, list(
      crossing_diagnostics = FALSE, dr_backend = "glm"
    )))
  )
  cuda_dr <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "logit",
    control = do.call(cf_control, utils::modifyList(common, list(
      crossing_diagnostics = FALSE, dr_backend = "cuda",
      gpu_backend = "cuda", gpu_precision = "float64",
      gpu_python = gpu_python, gpu_python_path = gpu_path
    )))
  )
  stopifnot(max(abs(
    cpu_dr$point$fits$group0$coefficients -
      cuda_dr$point$fits$group0$coefficients
  )) < 1e-5)
  stopifnot(max(abs(cpu_dr$point$effects - cuda_dr$point$effects)) < 1e-5)

  bootstrap_common <- list(
    nreg = 5L, reported_quantiles = c(0.25, 0.5, 0.75),
    legacy_qr_shift = FALSE, crossing_diagnostics = FALSE,
    bootstrap_progress = FALSE
  )
  cpu_boot_dir <- tempfile("gpu_test_cpu_boot_")
  cuda_boot_dir <- tempfile("gpu_test_cuda_boot_")
  cpu_boot <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "qr", solver = "fn",
    control = do.call(cf_control, c(bootstrap_common, list(
      gpu_backend = "cpu"
    ))),
    bootstrap_reps = 2L, bootstrap_workers = 1L,
    checkpoint_dir = cpu_boot_dir, seed = 20260809L
  )
  cuda_boot <- counterfactual_decompose(
    y ~ x, example, "group", "weight", model = "qr", solver = "fn",
    control = do.call(cf_control, c(bootstrap_common, list(
      gpu_backend = "cuda", gpu_precision = "float64",
      gpu_python = gpu_python, gpu_python_path = gpu_path
    ))),
    bootstrap_reps = 2L, point_workers = 2L, bootstrap_workers = 1L,
    checkpoint_dir = cuda_boot_dir, seed = 20260809L
  )
  stopifnot(
    max(abs(cpu_boot$point$effects - cuda_boot$point$effects)) < 1e-10,
    max(abs(cpu_boot$results$std_error -
      cuda_boot$results$std_error)) < 1e-10,
    identical(cuda_boot$metadata$bootstrap_point_workers, 2L),
    identical(cuda_boot$bootstrap$draw_point_workers, 2L)
  )
  unlink(cpu_boot_dir, recursive = TRUE)
  unlink(cuda_boot_dir, recursive = TRUE)
}
