## Empirical 99% quantile of the null Gaussian FOCUS "both" statistic, for
## series lengths n = 10, 20, ..., 1000. This gives gamma(n) directly from
## simulation at every length, without assuming any parametric scaling form
## (e.g. gamma = c * log(n)) -- read off (or fit) the right threshold for
## any n from the resulting table.
##
## "Both" here is the full-series statistic used by SVP's subtests = "both"
## pruning rule for test = "gaussian_mean": the max over every start s and
## every contained split of the single-change FOCUS statistic -- i.e.
## max_s max_{s < tau < t} Z_{s,tau,t}^2, computed by the same
## gaussian_critical_values_cpp() helper used to calibrate the main power
## study (simulations/power_gaussian/calibrate_modes.R).
##
## Run from the package root:
##   SVP_RUN_EXPERIMENT=true Rscript simulations/experiment/quantile_both_all_lengths.R
## Optionally set SVP_EXPERIMENT_WORKERS to parallelize (mclapply, not
## available on Windows).

library(svpChange)
library(Rcpp)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

EXPERIMENT_ROOT <- file.path("simulations", "experiment")
BOTH_QUANTILE_LENGTHS <- seq(10L, 1000L, by = 10L)
BOTH_QUANTILE_REPS <- 5000L
BOTH_QUANTILE_TARGET <- 0.99
BOTH_QUANTILE_SEED <- 314159L

compile_both_focus <- function() {
  include <- normalizePath("src", mustWork = TRUE)
  old_flags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  on.exit({
    if (is.na(old_flags)) Sys.unsetenv("PKG_CPPFLAGS") else
      Sys.setenv(PKG_CPPFLAGS = old_flags)
  }, add = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste("-I", shQuote(include)))
  Rcpp::sourceCpp(file.path(
    "simulations", "power_gaussian", "gaussian_critical_values.cpp"
  ))
}

experiment_workers <- function() {
  requested <- suppressWarnings(as.integer(Sys.getenv(
    "SVP_EXPERIMENT_WORKERS", "1"
  )))
  if (length(requested) != 1L || is.na(requested) || requested < 1L) 1L else requested
}

run_both_quantile_experiment <- function(
    lengths = BOTH_QUANTILE_LENGTHS,
    reps = BOTH_QUANTILE_REPS,
    target = BOTH_QUANTILE_TARGET,
    seed = BOTH_QUANTILE_SEED,
    workers = experiment_workers()) {
  dir.create(EXPERIMENT_ROOT, recursive = TRUE, showWarnings = FALSE)
  compile_both_focus()
  rows <- lapply(seq_along(lengths), function(j) {
    n <- lengths[[j]]
    values <- unlist(power_lapply(seq_len(reps), function(i) {
      set_power_seed(seed + j * 1000000L + i)
      gaussian_critical_values_cpp(stats::rnorm(n))[["both"]]
    }, workers = workers))
    q <- as.numeric(stats::quantile(values, probs = target, type = 1, names = FALSE))
    data.frame(
      n = n,
      gamma_both = q,
      c_both = q / log(n),
      reps = reps,
      target = target,
      mean = mean(values),
      sd = stats::sd(values),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$seed <- seed
  utils::write.csv(
    out, file.path(EXPERIMENT_ROOT, "both_focus_quantile_by_length.csv"),
    row.names = FALSE
  )
  saveRDS(out, file.path(EXPERIMENT_ROOT, "both_focus_quantile_by_length.rds"))
  out
}

## Fit gamma(n) = A*log(n) + B*log(log(n)) + C by ordinary least squares.
## The log-log term is motivated by Darling-Erdos-type boundary-crossing
## theory for the running maximum of a CUSUM process (see the "both_focus_lil"
## subtest's own iterated-logarithm penalty for the single-interval analogue);
## it fits this data far better than a plain gamma = c*log(n) form (lower
## AIC, higher R2).
fit_both_quantile_curve <- function(
    out = readRDS(file.path(
      EXPERIMENT_ROOT, "both_focus_quantile_by_length.rds"
    ))) {
  out$logn <- log(out$n)
  out$loglogn <- log(log(out$n))
  stats::lm(gamma_both ~ logn + loglogn, data = out)
}

plot_both_quantile_experiment <- function(
    out = readRDS(file.path(
      EXPERIMENT_ROOT, "both_focus_quantile_by_length.rds"
    )),
    fit = fit_both_quantile_curve(out)) {
  curve <- data.frame(n = seq(min(out$n), max(out$n), length.out = 400L))
  curve$logn <- log(curve$n)
  curve$loglogn <- log(log(curve$n))
  curve$gamma_fit <- stats::predict(fit, newdata = curve)
  curve$c_fit <- curve$gamma_fit / curve$logn

  coef_fit <- stats::coef(fit)
  r2 <- summary(fit)$r.squared
  fit_label <- sprintf(
    "gamma(n) = %.3f log(n) + %.3f log(log(n)) + %.3f  (R2 = %.4f)",
    coef_fit[["logn"]], coef_fit[["loglogn"]], coef_fit[["(Intercept)"]], r2
  )

  p1 <- ggplot2::ggplot(out, ggplot2::aes(n, gamma_both)) +
    ggplot2::geom_point(size = 0.9, colour = "#0072B2") +
    ggplot2::geom_line(
      data = curve, ggplot2::aes(n, gamma_fit),
      colour = "black", linetype = 2, linewidth = 0.7
    ) +
    ggplot2::labs(
      x = "Series length n", y = "Empirical 99% null quantile of 'both'",
      title = "gamma(n) for SVP both, from direct simulation",
      subtitle = fit_label
    ) +
    ggplot2::theme_minimal(base_size = 13)
  p2 <- ggplot2::ggplot(out, ggplot2::aes(n, c_both)) +
    ggplot2::geom_point(size = 0.9, colour = "#D55E00") +
    ggplot2::geom_line(
      data = curve, ggplot2::aes(n, c_fit),
      colour = "black", linetype = 2, linewidth = 0.7
    ) +
    ggplot2::labs(
      x = "Series length n", y = "Implied c = gamma(n) / log(n)",
      title = "How far a single c * log(n) formula would be from exact",
      subtitle = fit_label
    ) +
    ggplot2::theme_minimal(base_size = 13)
  ggplot2::ggsave(
    file.path(EXPERIMENT_ROOT, "both_focus_gamma_vs_n.pdf"), p1,
    width = 8, height = 5
  )
  ggplot2::ggsave(
    file.path(EXPERIMENT_ROOT, "both_focus_c_vs_n.pdf"), p2,
    width = 8, height = 5
  )
  invisible(list(gamma_plot = p1, c_plot = p2, fit = fit))
}

if (identical(tolower(Sys.getenv("SVP_RUN_EXPERIMENT")), "true")) {
  out <- run_both_quantile_experiment()
  print(out)
  plot_both_quantile_experiment(out)
}
