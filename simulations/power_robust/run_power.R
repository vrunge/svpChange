## Heavy-tail power study: corrected Student-t(2) Section 5.2 study.

library(svpChange)
library(changepoint)
library(robseg)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "power_robust", "calibrate_robust.R"))

ROBUST_ROOT <- file.path("simulations", "power_robust")

format_robust_parameter <- function(x) {
  trimws(formatC(x, format = "fg", digits = 6))
}

robust_method_labels <- function(calibration) {
  selected <- calibration$selected
  c(
    "PELT", "RFPOP",
    paste0("SVP MedianMood / right / alpha = ",
           format_robust_parameter(selected$mood_alpha)),
    paste0("SVP Wilcoxon / right / c = ",
           format_robust_parameter(selected$wilcoxon_constant))
  )
}

fit_robust_methods <- function(y, oracle_segments, calibration) {
  n <- length(y)
  selected <- calibration$selected
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
  )
  rfpop <- fit_rfpop(y, selected$rfpop_constant)
  mood_label <- robust_method_labels(calibration)[[3L]]
  wilcoxon_label <- robust_method_labels(calibration)[[4L]]
  boundaries <- list(
    PELT = normalise_boundaries(changepoint::cpts(pelt), n),
    RFPOP = rfpop$boundaries
  )
  boundaries[[mood_label]] <- normalise_boundaries(
    SVP(
      y, mood_threshold(n, oracle_segments, selected$mood_alpha),
      "MedianMoodCost", subtests = "right"
    )$changepoints, n
  )
  boundaries[[wilcoxon_label]] <- normalise_boundaries(
    SVP(
      y, wilcoxon_threshold(n, oracle_segments,
                            selected$wilcoxon_constant),
      "WilcoxonCost", subtests = "right"
    )$changepoints, n
  )
  fitted <- lapply(boundaries, function(x) fitted_piecewise_mean(y, x))
  fitted[["RFPOP"]] <- rfpop$fitted
  list(boundaries = boundaries, fitted = fitted)
}

run_robust_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 4, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L,
    calibration = NULL) {
  if (is.null(calibration)) {
    calibration <- if (file.exists(file.path(
      ROBUST_ROOT, "robust_calibration.rds"
    ))) {
      readRDS(file.path(ROBUST_ROOT, "robust_calibration.rds"))
    } else {
      run_robust_calibration(n = n, workers = workers)
    }
  }
  run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = function(size) ts_generator(
      chpts = size, parameters = 0, type = "student", df = 2, scale = 1
    ),
    fit_methods = function(y, oracle_segments) {
      fit_robust_methods(y, oracle_segments, calibration)
    },
    tolerance = round(n * 0.0025),
    workers = workers,
    seed = seed
  )
}

run_and_save_robust <- function(workers = power_default_workers()) {
  calibration_path <- file.path(ROBUST_ROOT, "robust_calibration.rds")
  calibration <- if (file.exists(calibration_path)) {
    readRDS(calibration_path)
  } else {
    run_robust_calibration(workers = workers)
  }
  results <- run_robust_power(workers = workers, calibration = calibration)
  save_power_outputs(
    results, ROBUST_ROOT,
    power_scenario_plot(
      1000L, 0.6, function(size) ts_generator(
        chpts = size, parameters = 0, type = "student", df = 2, scale = 1
      ),
      y_limits = c(-10, 10)
    ),
    plot_results = set_algorithm_order(
      results, robust_method_labels(calibration)
    )
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_robust()
}
