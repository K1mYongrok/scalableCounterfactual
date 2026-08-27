resolve_data_vector <- function(value, data, name) {
  if (is.character(value) && length(value) == 1L) {
    if (!value %in% names(data)) {
      stop(name, " column not found: ", value, call. = FALSE)
    }
    return(data[[value]])
  }
  if (length(value) != nrow(data)) {
    stop(name, " must be a column name or have nrow(data) elements", call. = FALSE)
  }
  value
}

resolve_optional_data_vector <- function(
    value, data, name, allow_scalar = FALSE, default = NULL) {
  if (is.null(value)) return(default)
  if (allow_scalar && length(value) == 1L && !is.character(value)) {
    return(rep(value, nrow(data)))
  }
  resolve_data_vector(value, data, name)
}

coerce_binary_group <- function(group) {
  if (is.factor(group)) group <- as.character(group)
  if (is.logical(group)) return(as.integer(group))
  if (is.character(group)) {
    if (!all(group %in% c("0", "1"))) {
      stop("group must contain only exact 0 and 1 values", call. = FALSE)
    }
    return(as.integer(group))
  }
  if (!is.numeric(group) || any(!is.finite(group)) ||
      !all(group %in% c(0, 1))) {
    stop("group must contain only exact 0 and 1 values", call. = FALSE)
  }
  as.integer(group)
}

reduce_common_design <- function(X0, X1, tolerance = 1e-10) {
  columns <- colnames(X0)
  intercept <- match("(Intercept)", columns)
  keep <- rep(TRUE, length(columns))

  has_variation <- function(x) {
    values <- range(x, finite = TRUE)
    length(values) == 2L && is.finite(values[[1L]]) &&
      is.finite(values[[2L]]) && (values[[2L]] - values[[1L]]) > tolerance
  }

  for (j in seq_along(columns)) {
    if (j != intercept &&
        (!has_variation(X0[, j]) || !has_variation(X1[, j]))) {
      keep[[j]] <- FALSE
    }
  }

  repeat {
    active <- which(keep)
    dependent <- integer()
    for (matrix_group in list(X0[, active, drop = FALSE],
                              X1[, active, drop = FALSE])) {
      qr_fit <- qr(matrix_group, tol = tolerance, LAPACK = FALSE)
      if (qr_fit$rank < ncol(matrix_group)) {
        aliased_local <- qr_fit$pivot[seq.int(qr_fit$rank + 1L,
                                             ncol(matrix_group))]
        dependent <- c(dependent, active[aliased_local])
      }
    }
    dependent <- setdiff(unique(dependent), intercept)
    if (!length(dependent)) break
    keep[dependent] <- FALSE
  }

  if (!keep[[intercept]]) {
    stop("the intercept cannot be removed from the common design", call. = FALSE)
  }
  if (sum(keep) < 2L) {
    stop("no identified covariates remain after rank reduction", call. = FALSE)
  }
  if (nrow(X0) <= sum(keep) || nrow(X1) <= sum(keep)) {
    stop("each group must have more observations than retained columns",
         call. = FALSE)
  }

  list(
    X0 = X0[, keep, drop = FALSE],
    X1 = X1[, keep, drop = FALSE],
    retained = columns[keep],
    dropped = columns[!keep]
  )
}

reduce_prepared_design <- function(prepared) {
  reduced <- reduce_common_design(prepared$X0, prepared$X1)
  prepared$X0 <- reduced$X0
  prepared$X1 <- reduced$X1
  prepared$design_columns <- reduced$retained
  prepared$dropped_design_columns <- unique(c(
    prepared$dropped_design_columns,
    reduced$dropped
  ))
  prepared
}

prepare_cf_data <- function(
    formula, data, group, weights = NULL, model = "qr",
    censoring = NULL, event = NULL) {
  if (!inherits(formula, "formula")) stop("formula must be a formula", call. = FALSE)
  data <- as.data.frame(data)
  model_frame <- stats::model.frame(
    formula,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = TRUE
  )
  terms_object <- stats::terms(model_frame)
  response <- stats::model.response(model_frame)
  if (!is.numeric(response)) {
    stop("the formula response must be numeric", call. = FALSE)
  }
  y <- as.numeric(response)
  X <- stats::model.matrix(terms_object, model_frame)
  group_vector <- resolve_data_vector(group, data, "group")
  weight_vector <- if (is.null(weights)) rep(1, nrow(data)) else {
    resolve_data_vector(weights, data, "weights")
  }
  censoring_vector <- resolve_optional_data_vector(
    censoring, data, "censoring", allow_scalar = TRUE
  )
  event_vector <- resolve_optional_data_vector(event, data, "event")
  if (model == "cqr" && is.null(censoring_vector)) {
    stop("censoring is required when model = 'cqr'", call. = FALSE)
  }
  if (model != "cqr" && !is.null(censoring_vector)) {
    stop("censoring is only valid when model = 'cqr'", call. = FALSE)
  }
  if (model != "cox" && !is.null(event_vector)) {
    stop("event is only valid when model = 'cox'", call. = FALSE)
  }
  if (model == "cox" && is.null(event_vector)) {
    event_vector <- rep(1L, nrow(data))
  }
  if (length(group_vector) != length(y) || length(weight_vector) != length(y)) {
    stop("formula, group, and weights must refer to the same rows", call. = FALSE)
  }
  keep <- stats::complete.cases(X, y, group_vector, weight_vector) &
    is.finite(y) & is.finite(weight_vector)
  if (!is.null(censoring_vector)) {
    if (!is.numeric(censoring_vector)) {
      stop("censoring must be numeric", call. = FALSE)
    }
    keep <- keep & is.finite(censoring_vector)
  }
  if (!is.null(event_vector)) {
    if (is.logical(event_vector)) event_vector <- as.integer(event_vector)
    if (!is.numeric(event_vector) || any(
      !is.na(event_vector) & !event_vector %in% c(0, 1)
    )) {
      stop("event must contain only 0 and 1", call. = FALSE)
    }
    keep <- keep & is.finite(event_vector)
  }
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  group_vector <- coerce_binary_group(group_vector[keep])
  weight_vector <- normalize_weights(weight_vector[keep])
  if (!is.null(censoring_vector)) censoring_vector <- censoring_vector[keep]
  if (!is.null(event_vector)) event_vector <- as.integer(event_vector[keep])
  if (!identical(sort(unique(group_vector)), 0:1)) {
    stop("group must contain both 0 and 1", call. = FALSE)
  }
  if (!"(Intercept)" %in% colnames(X)) {
    stop("the current implementation requires an intercept", call. = FALSE)
  }
  index0 <- which(group_vector == 0L)
  index1 <- which(group_vector == 1L)
  prepared <- structure(list(
    formula = formula,
    terms = terms_object,
    omitted_rows = sum(!keep),
    X0 = X[index0, , drop = FALSE],
    y0 = y[index0],
    w0 = weight_vector[index0],
    censoring0 = if (is.null(censoring_vector)) NULL else censoring_vector[index0],
    event0 = if (is.null(event_vector)) NULL else event_vector[index0],
    X1 = X[index1, , drop = FALSE],
    y1 = y[index1],
    w1 = weight_vector[index1],
    censoring1 = if (is.null(censoring_vector)) NULL else censoring_vector[index1],
    event1 = if (is.null(event_vector)) NULL else event_vector[index1],
    design_columns = colnames(X),
    dropped_design_columns = character(),
    n = nrow(X),
    n0 = length(index0),
    n1 = length(index1)
  ), class = "cf_prepared_data")
  reduce_prepared_design(prepared)
}

subsample_prepared_data <- function(prepared, sample_n, seed = 1L) {
  sample_n <- assert_scalar_integer(sample_n, "sample_n", 2L)
  if (sample_n >= prepared$n) return(prepared)
  set.seed(seed)
  minimum_per_group <- ncol(prepared$X0) + 1L
  if (sample_n < 2L * minimum_per_group) {
    stop(
      "sample_n must allocate at least ", minimum_per_group,
      " observations to each group for the retained design",
      call. = FALSE
    )
  }
  proportional_target1 <- as.integer(round(
    as.numeric(sample_n) * as.numeric(prepared$n1) / as.numeric(prepared$n)
  ))
  lower_target1 <- max(minimum_per_group, sample_n - prepared$n0)
  upper_target1 <- min(prepared$n1, sample_n - minimum_per_group)
  if (lower_target1 > upper_target1) {
    stop(
      "sample_n is incompatible with group sizes and the retained design",
      call. = FALSE
    )
  }
  target1 <- max(lower_target1, min(upper_target1, proportional_target1))
  target0 <- sample_n - target1
  i0 <- sample.int(prepared$n0, target0, replace = FALSE)
  i1 <- sample.int(prepared$n1, target1, replace = FALSE)
  prepared$X0 <- prepared$X0[i0, , drop = FALSE]
  prepared$y0 <- prepared$y0[i0]
  prepared$w0 <- prepared$w0[i0]
  if (!is.null(prepared$censoring0)) {
    prepared$censoring0 <- prepared$censoring0[i0]
  }
  if (!is.null(prepared$event0)) prepared$event0 <- prepared$event0[i0]
  prepared$X1 <- prepared$X1[i1, , drop = FALSE]
  prepared$y1 <- prepared$y1[i1]
  prepared$w1 <- prepared$w1[i1]
  if (!is.null(prepared$censoring1)) {
    prepared$censoring1 <- prepared$censoring1[i1]
  }
  if (!is.null(prepared$event1)) prepared$event1 <- prepared$event1[i1]
  prepared$n0 <- length(i0)
  prepared$n1 <- length(i1)
  prepared$n <- prepared$n0 + prepared$n1
  reduce_prepared_design(prepared)
}
