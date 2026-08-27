library(scalableCounterfactual)

# Frozen from code/vendor/Counterfactual_1.2/R/counterfactual.r using the
# audit harness code/validation/validate_counterfactual_1_2_parity.R.
set.seed(20260809L)
n <- 1200L
group <- rep(0:1, each = n / 2L)
x1 <- rnorm(n, mean = 0.2 * group)
x2 <- rbinom(n, 1, plogis(-0.4 + 0.5 * group))
weight <- runif(n, 0.5, 2)
error <- (0.7 + 0.15 * x1) * stats::rt(n, df = 7)
y <- 1 + 0.35 * x1 - 0.2 * x1^2 + 0.25 * x2 +
  group * (0.12 + 0.08 * x1) + error
data <- data.frame(y, x1, x2, group, weight)

fit <- counterfactual_decompose(
  y ~ x1 + I(x1^2) + x2,
  data,
  group = "group",
  weights = "weight",
  model = "qr",
  solver = "br",
  control = cf_control(
    nreg = 19L,
    trimming = 0.05,
    reported_quantiles = seq(0.1, 0.9, 0.1),
    legacy_qr_shift = TRUE,
    legacy_weighted_quantile = TRUE,
    qr_precondition = FALSE,
    marginal_method = "matrix",
    crossing_diagnostics = FALSE
  ),
  seed = 20260809L
)

reference <- rbind(
  structure = c(
    0.0811126201456794, 0.0903377086576258, 0.114027924237932,
    0.142878952027587, 0.149135955827721, 0.149413814920535,
    0.150823651432217, 0.156635294959477, 0.153343165568289
  ),
  composition = c(
    0.0281020541207119, 0.0184171151157444, 0.019575962818602,
    0.0309928677136879, 0.0339721935966129, 0.0419156741036746,
    0.0531751665322378, 0.063831365579148, 0.049639487971985
  ),
  total = c(
    0.109214674266391, 0.10875482377337, 0.133603887056534,
    0.173871819741275, 0.183108149424334, 0.191329489024209,
    0.203998817964455, 0.220466660538625, 0.202982653540274
  )
)
stopifnot(max(abs(fit$point$effects - reference)) < 1e-7)
