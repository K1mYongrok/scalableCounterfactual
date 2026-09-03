library(scalableCounterfactual)

ns <- asNamespace("scalableCounterfactual")
internal <- function(name) get(name, envir = ns, inherits = FALSE)

# Pointwise inference may use all finite draws at a quantile, but a simultaneous
# band is reported only on quantiles sharing the complete bootstrap support.
control <- cf_control(
  alpha = 0.05,
  robust_se = FALSE,
  reported_quantiles = c(0.25, 0.5, 0.75)
)
draws <- cbind(
  complete = c(0.1, 0.2, 0.3),
  partial = c(0.1, NA_real_, 0.3),
  constant = c(0.2, 0.2, 0.2)
)
inference <- internal("effect_inference")(
  effect = "total",
  estimate = c(0.2, 0.2, 0.2),
  draws = draws,
  control = control,
  model = "qr",
  solver = "br"
)
stopifnot(
  identical(inference$bootstrap_reps_effective, c(3L, 2L, 3L)),
  is.finite(inference$pointwise_lower[[2L]]),
  is.na(inference$uniform_lower[[2L]]),
  is.na(inference$uniform_upper[[2L]]),
  is.finite(inference$uniform_lower[[1L]]),
  identical(inference$uniform_lower[[3L]], 0.2),
  identical(inference$uniform_upper[[3L]], 0.2)
)

# Functional tests are invariant to a change in outcome units.
set.seed(811L)
scale_draws <- cbind(
  stats::rnorm(999L, sd = 0.001),
  stats::rnorm(999L, sd = 1),
  stats::rnorm(999L, sd = 10)
)
scale_observed <- c(0.003, 0.6, 3)
scale_centered <- sweep(scale_draws, 2L, colMeans(scale_draws), "-")
original_statistics <- internal("effect_test_statistics")(
  scale_centered, scale_observed, scale_draws, FALSE
)
rescaled_statistics <- internal("effect_test_statistics")(
  scale_centered * 1e-6, scale_observed * 1e-6,
  scale_draws * 1e-6, FALSE
)
stopifnot(isTRUE(all.equal(
  original_statistics, rescaled_statistics,
  tolerance = 1e-12, check.attributes = FALSE
)))

# Failed staged writes leave the previous managed output set intact.
transaction_dir <- tempfile("review_output_transaction_")
dir.create(transaction_dir)
old_path <- file.path(transaction_dir, "result.txt")
writeLines("old-valid-result", old_path)
old_hash <- unname(tools::md5sum(old_path))
transaction_error <- tryCatch(
  internal("atomic_write_output_files")(
    transaction_dir, "result.txt", required_files = "result.txt",
    writer = function(stage) {
      writeLines("incomplete-new-result", file.path(stage, "result.txt"))
      stop("synthetic write failure")
    }
  ),
  error = identity
)
stopifnot(
  inherits(transaction_error, "error"),
  identical(unname(tools::md5sum(old_path)), old_hash),
  identical(readLines(old_path), "old-valid-result")
)
directory_checkpoint <- file.path(transaction_dir, "checkpoint.rds")
dir.create(directory_checkpoint)
writeLines("sentinel", file.path(directory_checkpoint, "sentinel.txt"))
checkpoint_directory_error <- tryCatch(
  internal("atomic_save_rds")(list(ok = TRUE), directory_checkpoint),
  error = identity
)
stopifnot(
  inherits(checkpoint_directory_error, "error"),
  file.exists(file.path(directory_checkpoint, "sentinel.txt"))
)
unlink(transaction_dir, recursive = TRUE, force = TRUE)

# ylim and type supplied through ... must not collide with the panel setup.
plot_object <- structure(list(results = inference), class = "cfdecomp")
plot_path <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_path)
plot_error <- tryCatch({
  plot(
    plot_object,
    effects = "total",
    interval = "uniform",
    ylim = c(-1, 1),
    type = "l"
  )
  NULL
}, error = identity)
grDevices::dev.off()
unlink(plot_path)
stopifnot(is.null(plot_error))

# A one-quantile fit remains a valid decomposition result; functional curve
# tests are simply unavailable for it.
one_quantile_object <- structure(list(
  point = list(quantiles = 0.5),
  bootstrap = list(effects = stats::setNames(
    replicate(3L, matrix(c(0.1, 0.2), ncol = 1L), simplify = FALSE),
    c("structure", "composition", "total")
  )),
  control = control
), class = "cfdecomp")
functional <- suppressWarnings(functional_effect_tests(one_quantile_object))
stopifnot(is.data.frame(functional), nrow(functional) == 0L)

# The legacy Hmisc path uses the original bootstrap sample size as its weight
# denominator. Collapsing repeated rows to frequency weights must therefore
# agree with explicitly expanding the same bootstrap sample.
stopifnot(identical(
  internal("weighted_quantile")(
    c(0, 10), c(2, 1), 0.5, normalization_n = 3
  ),
  internal("weighted_quantile")(c(0, 0, 10), rep(1, 3), 0.5)
))
set.seed(4102)
n_group <- 45L
legacy_data <- data.frame(
  y = c(rnorm(n_group), 0.25 + rnorm(n_group)),
  x = rnorm(2L * n_group),
  group = rep(0:1, each = n_group),
  weight = runif(2L * n_group, 0.5, 2)
)
legacy_control <- cf_control(
  nreg = 5L,
  reported_quantiles = c(0.25, 0.5, 0.75),
  bootstrap_scheme = "empirical",
  legacy_qr_shift = FALSE,
  legacy_weighted_quantile = TRUE,
  bootstrap_progress = FALSE,
  crossing_diagnostics = FALSE,
  marginal_method = "matrix"
)
prepared <- internal("prepare_cf_data")(
  y ~ x, legacy_data, "group", "weight"
)
point <- internal("estimate_point_prepared")(
  prepared,
  model = "qr",
  solver = "br",
  control = legacy_control,
  point_workers = 1L,
  point_seed = 902L,
  keep_fits = TRUE
)
common <- list(
  prepared = prepared,
  point = point,
  model = "qr",
  solver = "br",
  control = legacy_control,
  seed = 902L,
  signature = "review-legacy-denominator",
  data_fingerprint = "review-data",
  bootstrap_engine = "standard",
  draw_point_workers = 1L,
  run_dir = tempfile("review_bootstrap_")
)
collapsed <- internal("bootstrap_replication_attempt")(1L, 0L, common)
seed_used <- internal("bootstrap_attempt_seed")(902L, 1L, 0L)
set.seed(seed_used)
multiplier0 <- internal("bootstrap_multipliers")(
  prepared$n0, legacy_control$bootstrap_scheme
)
multiplier1 <- internal("bootstrap_multipliers")(
  prepared$n1, legacy_control$bootstrap_scheme
)
index0 <- rep.int(seq_len(prepared$n0), multiplier0)
index1 <- rep.int(seq_len(prepared$n1), multiplier1)
expanded <- prepared
expanded$X0 <- prepared$X0[index0, , drop = FALSE]
expanded$y0 <- prepared$y0[index0]
expanded$w0 <- prepared$w0[index0]
expanded$X1 <- prepared$X1[index1, , drop = FALSE]
expanded$y1 <- prepared$y1[index1]
expanded$w1 <- prepared$w1[index1]
expanded$n0 <- length(index0)
expanded$n1 <- length(index1)
expanded$n <- expanded$n0 + expanded$n1
expanded <- internal("reduce_prepared_design")(expanded)
expanded_point <- internal("estimate_point_prepared")(
  expanded,
  model = "qr",
  solver = "br",
  control = legacy_control,
  point_workers = 1L,
  point_seed = internal("bootstrap_attempt_seed")(902L, 1L, 0L, stream = 1L),
  keep_fits = FALSE
)
stopifnot(
  isTRUE(all.equal(
    collapsed$effects,
    expanded_point$effects,
    tolerance = 1e-10,
    check.attributes = FALSE
  )),
  identical(collapsed$group0_active_rows, as.integer(sum(multiplier0 > 0))),
  identical(collapsed$group1_active_rows, as.integer(sum(multiplier1 > 0))),
  identical(
    collapsed$active_rows,
    as.integer(sum(multiplier0 > 0) + sum(multiplier1 > 0))
  )
)

# The default stable type-7 path also treats integer frequency weights as the
# collapsed representation of the explicit bootstrap sample. This applies to
# fitted residual distributions as well as final marginalization.
default_control <- cf_control(
  nreg = 5L,
  reported_quantiles = c(0.25, 0.5, 0.75),
  bootstrap_scheme = "empirical",
  legacy_qr_shift = FALSE,
  legacy_weighted_quantile = FALSE,
  bootstrap_progress = FALSE,
  crossing_diagnostics = FALSE,
  marginal_method = "matrix"
)
for (bootstrap_model in c("qr", "loc", "locsca")) {
  bootstrap_solver <- if (bootstrap_model == "qr") "br" else NA_character_
  model_point <- internal("estimate_point_prepared")(
    prepared, model = bootstrap_model, solver = bootstrap_solver,
    control = default_control, point_workers = 1L, point_seed = 902L,
    keep_fits = TRUE
  )
  model_common <- common
  model_common$point <- model_point
  model_common$model <- bootstrap_model
  model_common$solver <- bootstrap_solver
  model_common$control <- default_control
  model_common$signature <- paste0("review-frequency-", bootstrap_model)
  model_collapsed <- internal("bootstrap_replication_attempt")(
    1L, 0L, model_common
  )
  model_expanded <- internal("estimate_point_prepared")(
    expanded, model = bootstrap_model, solver = bootstrap_solver,
    control = default_control, point_workers = 1L,
    point_seed = internal("bootstrap_attempt_seed")(
      902L, 1L, 0L, stream = 1L
    ),
    keep_fits = FALSE
  )
  stopifnot(isTRUE(all.equal(
    model_collapsed$effects, model_expanded$effects,
    tolerance = 1e-9, check.attributes = FALSE
  )))
}

# The standard bootstrap may use the process-aware one-step solver explicitly.
# Its compressed frequencies must follow the same numerical path as literal
# row replication, including bandwidth and Jacobian sample-size calculations.
onestep_control <- cf_control(
  nreg = 21L,
  reported_quantiles = c(0.25, 0.5, 0.75),
  bootstrap_scheme = "empirical",
  legacy_qr_shift = FALSE,
  bootstrap_progress = FALSE,
  crossing_diagnostics = FALSE,
  marginal_method = "matrix",
  onestep_first_solver = "fn"
)
onestep_point <- suppressWarnings(internal("estimate_point_prepared")(
  prepared, model = "qr", solver = "onestep", control = onestep_control,
  point_workers = 1L, point_seed = 902L, keep_fits = TRUE
))
onestep_common <- common
onestep_common$point <- onestep_point
onestep_common$solver <- "onestep"
onestep_common$control <- onestep_control
onestep_common$signature <- "review-frequency-onestep-standard"
onestep_collapsed <- suppressWarnings(
  internal("bootstrap_replication_attempt")(1L, 0L, onestep_common)
)
onestep_expanded <- suppressWarnings(internal("estimate_point_prepared")(
  expanded, model = "qr", solver = "onestep", control = onestep_control,
  point_workers = 1L,
  point_seed = internal("bootstrap_attempt_seed")(
    902L, 1L, 0L, stream = 1L
  ),
  keep_fits = TRUE
))
stopifnot(
  isTRUE(all.equal(
    onestep_collapsed$effects, onestep_expanded$effects,
    tolerance = 1e-12, check.attributes = FALSE
  )),
  all(vapply(
    onestep_expanded$fits,
    function(fit) identical(fit$frequency_expanded, FALSE),
    logical(1L)
  ))
)

# The saved-Jacobian engine may retain exactly p distinct, full-rank rows when
# bootstrap multiplicities imply an effective sample larger than p.
tiny_onestep_point <- list(
  taus = 0.5,
  onestep_traversal = 1L,
  onestep_first_solver_resolved = "fn",
  onestep_inverse_jacobians = list(NULL),
  solver_requested = "onestep",
  onestep_implementation_version = "test",
  stata_source_version = "test",
  stata_source_commit = "test"
)
tiny_onestep_X <- cbind("(Intercept)" = 1, x = c(-1, 1))
tiny_onestep_fit <- internal("fit_bootstrap_onestep")(
  tiny_onestep_X, c(0, 2), c(1, 1), c(2, 2), tiny_onestep_point
)
tiny_expanded_index <- rep(1:2, each = 2L)
tiny_onestep_expanded <- internal("fit_bootstrap_onestep")(
  tiny_onestep_X[tiny_expanded_index, , drop = FALSE],
  c(0, 2)[tiny_expanded_index], rep(1, 4), rep(1, 4),
  tiny_onestep_point
)
stopifnot(
  all(is.finite(tiny_onestep_fit$coefficients)),
  isTRUE(all.equal(
    tiny_onestep_fit$coefficients, tiny_onestep_expanded$coefficients,
    tolerance = 1e-12, check.attributes = FALSE
  ))
)
rank_deficient_onestep <- try(
  internal("fit_bootstrap_onestep")(
    cbind("(Intercept)" = 1, x = c(1, 1)), c(0, 2), c(1, 1), c(2, 2),
    tiny_onestep_point
  ),
  silent = TRUE
)
too_small_onestep <- try(
  internal("fit_bootstrap_onestep")(
    tiny_onestep_X, c(0, 2), c(1, 1), c(1, 1), tiny_onestep_point
  ),
  silent = TRUE
)
stopifnot(
  inherits(rank_deficient_onestep, "try-error"),
  inherits(too_small_onestep, "try-error")
)

# The exported low-level API validates multiplicities and preserves its old
# behavior when frequency is omitted.
frequency_api_X <- cbind("(Intercept)" = 1, x = seq(-1, 1, length.out = 24L))
frequency_api_y <- 0.5 * frequency_api_X[, "x"] +
  sin(seq_len(nrow(frequency_api_X)))
frequency_api_w <- seq(0.5, 1.5, length.out = nrow(frequency_api_X))
frequency_api_tau <- seq(0.1, 0.9, by = 0.05)
frequency_api_fit <- suppressWarnings(fit_weighted_qr(
  frequency_api_X, frequency_api_y, frequency_api_w, frequency_api_tau,
  solver = "onestep", onestep_first_solver = "fn"
))
frequency_api_default <- suppressWarnings(fit_weighted_qr(
  frequency_api_X, frequency_api_y, frequency_api_w, frequency_api_tau,
  solver = "onestep", onestep_first_solver = "fn", frequency = rep(1, 24L)
))
stopifnot(
  identical(frequency_api_fit$coefficients,
            frequency_api_default$coefficients),
  inherits(try(fit_weighted_qr(
    frequency_api_X, frequency_api_y, frequency_api_w, 0.5,
    solver = "fn", frequency = rep(1, 23L)
  ), silent = TRUE), "try-error"),
  inherits(try(fit_weighted_qr(
    frequency_api_X, frequency_api_y, frequency_api_w, 0.5,
    solver = "fn", frequency = c(0, rep(1, 23L))
  ), silent = TRUE), "try-error"),
  inherits(try(fit_weighted_qr(
    frequency_api_X, frequency_api_y, frequency_api_w, 0.5,
    solver = "fn", frequency = c(1.5, rep(1, 23L))
  ), silent = TRUE), "try-error")
)

set.seed(4103)
cqr_group_n <- 100L
cqr_data <- data.frame(
  y = c(
    rep(0, 10L), exp(rnorm(90L)),
    rep(0, 10L), exp(0.15 + rnorm(90L))
  ),
  x = rnorm(2L * cqr_group_n),
  group = rep(0:1, each = cqr_group_n),
  weight = runif(2L * cqr_group_n, 0.5, 2)
)
cqr_prepared <- internal("prepare_cf_data")(
  y ~ x, cqr_data, "group", "weight", model = "cqr", censoring = 0
)
cqr_control <- cf_control(
  nreg = 3L, reported_quantiles = 0.5,
  bootstrap_scheme = "empirical", legacy_qr_shift = FALSE,
  bootstrap_progress = FALSE, crossing_diagnostics = FALSE,
  marginal_method = "matrix"
)
cqr_point <- suppressWarnings(internal("estimate_point_prepared")(
  cqr_prepared, model = "cqr", solver = "fn", control = cqr_control,
  point_workers = 1L, point_seed = 903L, keep_fits = TRUE
))
cqr_common <- list(
  prepared = cqr_prepared, point = cqr_point, model = "cqr", solver = "fn",
  control = cqr_control, seed = 903L, signature = "review-cqr-frequency",
  data_fingerprint = "review-cqr-data", bootstrap_engine = "standard",
  draw_point_workers = 1L, run_dir = tempfile("review_cqr_bootstrap_")
)
cqr_collapsed <- suppressWarnings(internal("bootstrap_replication_attempt")(
  1L, 0L, cqr_common
))
set.seed(internal("bootstrap_attempt_seed")(903L, 1L, 0L))
cqr_multiplier0 <- internal("bootstrap_multipliers")(
  cqr_prepared$n0, "empirical"
)
cqr_multiplier1 <- internal("bootstrap_multipliers")(
  cqr_prepared$n1, "empirical"
)
cqr_index0 <- rep.int(seq_len(cqr_prepared$n0), cqr_multiplier0)
cqr_index1 <- rep.int(seq_len(cqr_prepared$n1), cqr_multiplier1)
cqr_expanded <- cqr_prepared
for (field in c("X", "y", "w", "censoring")) {
  source0 <- cqr_prepared[[paste0(field, "0")]]
  source1 <- cqr_prepared[[paste0(field, "1")]]
  cqr_expanded[[paste0(field, "0")]] <- if (field == "X") {
    source0[cqr_index0, , drop = FALSE]
  } else source0[cqr_index0]
  cqr_expanded[[paste0(field, "1")]] <- if (field == "X") {
    source1[cqr_index1, , drop = FALSE]
  } else source1[cqr_index1]
}
cqr_expanded$n0 <- length(cqr_index0)
cqr_expanded$n1 <- length(cqr_index1)
cqr_expanded$n <- cqr_expanded$n0 + cqr_expanded$n1
cqr_expanded <- internal("reduce_prepared_design")(cqr_expanded)
cqr_expanded_point <- suppressWarnings(internal("estimate_point_prepared")(
  cqr_expanded, model = "cqr", solver = "fn", control = cqr_control,
  point_workers = 1L,
  point_seed = internal("bootstrap_attempt_seed")(
    903L, 1L, 0L, stream = 1L
  ),
  keep_fits = FALSE
))
cqr_frequency_difference <- max(abs(
  cqr_collapsed$effects - cqr_expanded_point$effects
))
if (cqr_frequency_difference >= 1e-5) {
  stop(
    "collapsed CQR bootstrap differs from explicit expansion by ",
    signif(cqr_frequency_difference, 6)
  )
}

# A matching-signature but malformed bootstrap draw is recomputed rather than
# admitted into inference.
dir.create(common$run_dir, recursive = TRUE)
first_draw <- internal("bootstrap_replication")(7L, common)
stopifnot(first_draw$status == "ok", identical(first_draw$cached, FALSE))
first_checkpoint <- readRDS(first_draw$path)
first_checkpoint$effects <- first_checkpoint$effects[-1L, , drop = FALSE]
saveRDS(first_checkpoint, first_draw$path, compress = FALSE)
repaired_draw <- internal("bootstrap_replication")(7L, common)
repaired_checkpoint <- readRDS(repaired_draw$path)
stopifnot(
  repaired_draw$status == "ok",
  identical(repaired_draw$cached, FALSE),
  identical(
    dim(repaired_checkpoint$effects),
    c(3L, length(common$control$reported_quantiles))
  ),
  identical(
    rownames(repaired_checkpoint$effects),
    c("structure", "composition", "total")
  )
)
checkpoint_mutations <- list(
  function(x) { x$effects[] <- NA_real_; x },
  function(x) { x$elapsed_seconds <- NA_real_; x },
  function(x) { x$peak_r_heap_mb <- -Inf; x },
  function(x) {
    x$group0_active_rows <- -1L
    x$active_rows <- x$group0_active_rows + x$group1_active_rows
    x
  },
  function(x) { x$active_rows <- x$active_rows + 1L; x },
  function(x) { x$draw_point_workers <- 0L; x }
)
for (mutate_checkpoint in checkpoint_mutations) {
  malformed_checkpoint <- mutate_checkpoint(unserialize(serialize(
    repaired_checkpoint, NULL
  )))
  saveRDS(malformed_checkpoint, first_draw$path, compress = FALSE)
  stopifnot(is.null(internal("valid_bootstrap_checkpoint")(
    first_draw$path, common, 7L
  )))
}
unlink(common$run_dir, recursive = TRUE, force = TRUE)

runtime <- internal("bootstrap_runtime_identity")(
  "qr", "br", legacy_control, point
)
stopifnot(
  identical(runtime$R, R.version.string),
  identical(runtime$platform, R.version$platform),
  all(c("Hmisc", "quantreg") %in% names(runtime$packages)),
  length(runtime$blas) == 1L,
  identical(runtime$lapack, as.character(base::La_version())),
  is.null(runtime$gpu)
)
default_runtime <- internal("bootstrap_runtime_identity")(
  "qr", "br", cf_control(), point
)
stopifnot(!"Hmisc" %in% names(default_runtime$packages))
gpu_module <- tempfile(fileext = ".py")
writeLines("def cuda_status(): return {'available': False}", gpu_module)
gpu_control <- legacy_control
gpu_control$gpu_backend <- "cuda"
gpu_control$gpu_module_path <- gpu_module
gpu_identity_1 <- internal("bootstrap_gpu_runtime_identity")(
  gpu_control, "qr", "br"
)
writeLines(c(
  "def cuda_status(): return {'available': False}",
  "MODULE_REVISION = 2"
), gpu_module)
gpu_identity_2 <- internal("bootstrap_gpu_runtime_identity")(
  gpu_control, "qr", "br"
)
unlink(gpu_module)
stopifnot(!identical(
  internal("object_md5")(gpu_identity_1),
  internal("object_md5")(gpu_identity_2)
))

# Simulation seeds and cached results are properties of scenario/replication,
# not of the order in which scenarios happen to be requested.
seed_task <- list(
  scenario = "location_shift", replication = 2L, master_seed = 7401L,
  task_id = 1L
)
same_seed_task <- seed_task
same_seed_task$task_id <- 999L
stopifnot(identical(
  internal("simulation_attempt_seed")(seed_task, 0L),
  internal("simulation_attempt_seed")(same_seed_task, 0L)
))
stopifnot(
  identical(internal("simulation_mean_elapsed")(c(1, NA, 3)), 2),
  is.na(internal("simulation_mean_elapsed")(c(NA_real_, Inf)))
)

shared_checkpoints <- tempfile("review_simulation_shared_")
fresh_checkpoints <- tempfile("review_simulation_fresh_")
simulation_arguments <- list(
  replications = 1L,
  n_per_group = 40L,
  solver = "fn",
  nreg = 5L,
  reported_quantiles = c(0.25, 0.5, 0.75),
  weighted = FALSE,
  bootstrap_reps = 0L,
  workers = 1L,
  seed = 7401L,
  progress = FALSE
)
forward <- do.call(simulate_counterfactual_validation, c(
  simulation_arguments,
  list(
    scenarios = c("location_shift", "composition_shift"),
    checkpoint_dir = shared_checkpoints
  )
))
reverse_cached <- do.call(simulate_counterfactual_validation, c(
  simulation_arguments,
  list(
    scenarios = c("composition_shift", "location_shift"),
    checkpoint_dir = shared_checkpoints
  )
))
reverse_fresh <- do.call(simulate_counterfactual_validation, c(
  simulation_arguments,
  list(
    scenarios = c("composition_shift", "location_shift"),
    checkpoint_dir = fresh_checkpoints
  )
))
canonical_raw <- function(object) {
  output <- object$raw[order(
    object$raw$scenario,
    object$raw$replication,
    object$raw$effect,
    object$raw$quantile
  ), c("scenario", "replication", "effect", "quantile", "estimate")]
  rownames(output) <- NULL
  output
}
stopifnot(
  reverse_cached$resources$cached_tasks == 2L,
  isTRUE(all.equal(canonical_raw(forward), canonical_raw(reverse_fresh))),
  isTRUE(all.equal(
    canonical_raw(reverse_cached), canonical_raw(reverse_fresh)
  ))
)
unlink(shared_checkpoints, recursive = TRUE, force = TRUE)
unlink(fresh_checkpoints, recursive = TRUE, force = TRUE)

# Matching error checkpoints record a prior attempt but must not permanently
# cache a transient failure. A resumed task recomputes and replaces them.
recovery_checkpoints <- tempfile("review_simulation_recovery_")
invisible(do.call(simulate_counterfactual_validation, c(
  simulation_arguments,
  list(scenarios = "location_shift", checkpoint_dir = recovery_checkpoints)
)))
recovery_file <- list.files(
  recovery_checkpoints, pattern = "[.]rds$", full.names = TRUE,
  recursive = TRUE
)[[1L]]
failed_checkpoint <- readRDS(recovery_file)
failed_checkpoint$status <- "error"
failed_checkpoint$error <- "synthetic transient failure"
saveRDS(failed_checkpoint, recovery_file, compress = FALSE)
recovered <- do.call(simulate_counterfactual_validation, c(
  simulation_arguments,
  list(scenarios = "location_shift", checkpoint_dir = recovery_checkpoints)
))
stopifnot(
  recovered$resources$successful_tasks == 1L,
  recovered$resources$cached_tasks == 0L,
  identical(readRDS(recovery_file)$status, "ok")
)
# Structural corruption with an otherwise matching signature is also repaired.
malformed_simulation <- readRDS(recovery_file)
malformed_simulation$rows <- malformed_simulation$rows[-1L, , drop = FALSE]
saveRDS(malformed_simulation, recovery_file, compress = FALSE)
repaired_simulation <- do.call(simulate_counterfactual_validation, c(
  simulation_arguments,
  list(scenarios = "location_shift", checkpoint_dir = recovery_checkpoints)
))
stopifnot(
  repaired_simulation$resources$cached_tasks == 0L,
  nrow(readRDS(recovery_file)$rows) ==
    3L * length(simulation_arguments$reported_quantiles)
)
unlink(recovery_checkpoints, recursive = TRUE, force = TRUE)

# A reference solver with no successful repetition must fail explicitly rather
# than returning a summary whose relative timing is entirely missing.
failing_solver <- "review_inference_failure_solver"
register_qr_solver(
  failing_solver,
  fit = function(X, y, weights, taus, control) {
    stop("intentional review failure", call. = FALSE)
  },
  exact = TRUE,
  process_aware = FALSE
)
benchmark_error <- tryCatch(
  benchmark_qr_solvers_repeated(
    y ~ x,
    data = legacy_data,
    group = "group",
    weights = "weight",
    solvers = failing_solver,
    reference_solver = failing_solver,
    control = cf_control(
      nreg = 3L,
      reported_quantiles = 0.5,
      crossing_diagnostics = FALSE
    ),
    repetitions = 1L,
    randomize_order = FALSE
  ),
  error = identity
)
unregister_qr_solver(failing_solver)
stopifnot(
  inherits(benchmark_error, "error"),
  grepl("failed in every recorded repetition", conditionMessage(benchmark_error))
)

# warmup affects the work represented by a resumable scaling checkpoint and is
# therefore part of its compatibility signature.
if (requireNamespace("processx", quietly = TRUE)) {
  scaling_checkpoint <- tempfile(fileext = ".rds")
  scaling_control <- cf_control(
    nreg = 3L,
    reported_quantiles = 0.5,
    crossing_diagnostics = FALSE
  )
  invisible(benchmark_qr_scaling(
    y ~ x,
    data = legacy_data,
    group = "group",
    weights = "weight",
    solvers = "fn",
    reference_solver = "fn",
    control = scaling_control,
    sample_sizes = nrow(legacy_data),
    repetitions = 1L,
    warmup = 0L,
    randomize_order = FALSE,
    checkpoint_path = scaling_checkpoint
  ))
  scaling_error <- tryCatch(
    benchmark_qr_scaling(
      y ~ x,
      data = legacy_data,
      group = "group",
      weights = "weight",
      solvers = "fn",
      reference_solver = "fn",
      control = scaling_control,
      sample_sizes = nrow(legacy_data),
      repetitions = 1L,
      warmup = 1L,
      randomize_order = FALSE,
      checkpoint_path = scaling_checkpoint,
      resume = TRUE
    ),
    error = identity
  )
  unlink(scaling_checkpoint)
  stopifnot(
    inherits(scaling_error, "error"),
    grepl("checkpoint is incompatible", conditionMessage(scaling_error))
  )

  # Resume uses the same absolute repetition seeds and solver order as an
  # uninterrupted run, even after warmup has already completed.
  full_checkpoint <- tempfile(fileext = ".rds")
  partial_checkpoint <- tempfile(fileext = ".rds")
  full_scaling <- benchmark_qr_scaling(
    y ~ x,
    data = legacy_data,
    group = "group",
    weights = "weight",
    solvers = c("fn", "qfnb"),
    reference_solver = "fn",
    control = scaling_control,
    sample_sizes = nrow(legacy_data),
    repetitions = 2L,
    warmup = 1L,
    randomize_order = TRUE,
    checkpoint_path = full_checkpoint,
    seed = 2L
  )
  partial <- readRDS(full_checkpoint)
  keep_rows <- vapply(
    partial$rows, function(row) identical(row$repetition, 1L), logical(1L)
  )
  partial$rows <- partial$rows[keep_rows]
  keep_keys <- vapply(strsplit(names(partial$fits), "::", fixed = TRUE), function(x) {
    identical(x[[2L]], "1")
  }, logical(1L))
  partial$fits <- partial$fits[keep_keys]
  saveRDS(partial, partial_checkpoint, compress = FALSE)
  resumed_scaling <- benchmark_qr_scaling(
    y ~ x,
    data = legacy_data,
    group = "group",
    weights = "weight",
    solvers = c("fn", "qfnb"),
    reference_solver = "fn",
    control = scaling_control,
    sample_sizes = nrow(legacy_data),
    repetitions = 2L,
    warmup = 1L,
    randomize_order = TRUE,
    checkpoint_path = partial_checkpoint,
    seed = 2L
  )
  order_columns <- c("sample_n", "solver", "repetition", "run_order")
  canonical_order <- function(raw) {
    raw <- raw[order(raw$repetition, raw$run_order), order_columns]
    rownames(raw) <- NULL
    raw
  }
  stopifnot(identical(
    canonical_order(full_scaling$raw), canonical_order(resumed_scaling$raw)
  ))

  # Plausibly shaped but semantically corrupt scaling state is rejected.
  valid_scaling <- readRDS(partial_checkpoint)
  scaling_validator <- internal("valid_scaling_checkpoint")
  scaling_validator_args <- list(
    signature = valid_scaling$signature,
    sample_sizes = nrow(legacy_data),
    solvers = c("fn", "qfnb"),
    repetitions = 2L,
    reported_quantiles = scaling_control$reported_quantiles,
    conditional_quantiles = scaling_control$conditional_quantiles,
    design_columns = 2L,
    point_workers_requested = 1L,
    point_workers = 1L
  )
  scaling_mutations <- list(
    function(x) { x$rows[[1L]]$run_order <- -99; x },
    function(x) { x$rows[[1L]]$elapsed_seconds <- NA_real_; x },
    function(x) { x$rows[[1L]]$peak_r_heap_mb <- -Inf; x },
    function(x) { x$rows[[1L]]$point_workers_requested <- -5L; x },
    function(x) { x$rows[[1L]]$point_workers <- "bogus"; x },
    function(x) {
      key <- names(x$fits)[[1L]]
      x$fits[[key]]$effects[] <- NA_real_
      x
    },
    function(x) {
      key <- names(x$fits)[[1L]]
      x$fits[[key]]$group0_coefficients[] <- NA_real_
      x
    }
  )
  for (mutate_checkpoint in scaling_mutations) {
    candidate <- mutate_checkpoint(unserialize(serialize(valid_scaling, NULL)))
    stopifnot(!do.call(
      scaling_validator, c(list(checkpoint = candidate), scaling_validator_args)
    ))
  }

  # Matching-signature malformed scaling state is discarded and rebuilt.
  malformed <- unserialize(serialize(valid_scaling, NULL))
  malformed$rows[[1L]]$run_order <- -99
  malformed$rows[[1L]]$elapsed_seconds <- NA_real_
  malformed$rows[[1L]]$peak_r_heap_mb <- -Inf
  malformed$rows[[1L]]$point_workers_requested <- -5L
  key <- names(malformed$fits)[[1L]]
  malformed$fits[[key]]$effects[] <- NA_real_
  saveRDS(malformed, partial_checkpoint, compress = FALSE)
  rebuilt_scaling <- suppressWarnings(benchmark_qr_scaling(
    y ~ x,
    data = legacy_data,
    group = "group",
    weights = "weight",
    solvers = c("fn", "qfnb"),
    reference_solver = "fn",
    control = scaling_control,
    sample_sizes = nrow(legacy_data),
    repetitions = 2L,
    warmup = 1L,
    randomize_order = TRUE,
    checkpoint_path = partial_checkpoint,
    seed = 2L
  ))
  stopifnot(nrow(rebuilt_scaling$raw) == 4L)
  unlink(c(full_checkpoint, partial_checkpoint))
}

# Missing RSS is a supported processx platform outcome, not negative infinity.
rss_raw <- data.frame(
  sample_n = c(100L, 100L), solver = c("fn", "fn"), status = c("ok", "ok"),
  elapsed_seconds = c(1, 2), process_elapsed_seconds = c(1.1, 2.1),
  peak_process_rss_mb = c(NA_real_, NA_real_),
  max_abs_effect_difference = c(0, 0)
)
rss_summary <- internal("summarize_scaling_benchmark")(rss_raw, "fn")
stopifnot(
  is.na(rss_summary$median_peak_process_rss_mb),
  is.na(rss_summary$max_peak_process_rss_mb),
  is.na(rss_summary$max_observed_process_rss_mb)
)

cat("review inference/bootstrap regression tests passed\n")
