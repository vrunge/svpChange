## Corrected Section 5.3 AR(1) power study.

library(svpChange)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_ar1", "calibrate_ar1.R"))

AR1_ROOT <- file.path("simulations", "power_ar1")
AR1_JUMP_SIZES <- seq(1, 5, 0.2)
AR1_REPS <- 100L

fit_ar1_methods <- function(y, rho = AR1_RHO, calibration) {
  if (missing(calibration) || is.null(calibration)) {
    calibration_path <- file.path(AR1_ROOT, "ar1_calibration.rds")
    if (!file.exists(calibration_path)) {
      stop("AR(1) calibration is required before fitting the power study")
    }
    calibration <- readRDS(calibration_path)
  }
  n <- length(y)
  selected <- calibration$selected
  fits <- list(
    pelt_ar1_approximate_boundaries(
      y, rho, penalty = selected$approximate_constant * log(n)
    ),
    pelt_inflated_boundaries(y, rho),
    decafs_ar1_boundaries(
      y, rho, penalty = selected$decafs_constant * log(n)
    ),
    svp_ar1focus_boundaries(
      y, rho, constant = selected$svp_constant, cost = "ar1"
    )
  )
  stats::setNames(fits, ar1_method_labels(calibration))
}

run_ar1_power <- function(
    n = AR1_N,
    rho = AR1_RHO,
    jump_sizes = AR1_JUMP_SIZES,
    reps = AR1_REPS,
    workers = power_default_workers(),
    seed = 123L,
    calibration = NULL) {
  require_decafs()
  if (is.null(calibration)) {
    calibration_path <- file.path(AR1_ROOT, "ar1_calibration.rds")
    calibration <- if (file.exists(calibration_path)) {
      readRDS(calibration_path)
    } else {
      run_ar1_calibration(n = n, rho = rho, workers = workers)
    }
  }
  results <- run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = function(size) simulate_ar1_noise(size, rho),
    fit_methods = function(y, true_segments) {
      fit_ar1_methods(y, rho = rho, calibration = calibration)
    },
    tolerance = AR1_TOLERANCE,
    workers = workers,
    seed = seed
  )
  expected_rows <- length(POWER_PATTERNS) * length(jump_sizes) * reps * 4L
  if (nrow(results) != expected_rows) {
    stop("AR(1) power study returned ", nrow(results),
         " rows; expected ", expected_rows)
  }
  results
}

run_and_save_ar1 <- function(
    workers = power_default_workers(), calibration = NULL) {
  require_decafs()
  if (is.null(calibration)) {
    calibration_path <- file.path(AR1_ROOT, "ar1_calibration.rds")
    calibration <- if (file.exists(calibration_path)) {
      readRDS(calibration_path)
    } else {
      run_ar1_calibration(workers = workers)
    }
  }
  results <- run_ar1_power(workers = workers, calibration = calibration)
  save_power_outputs(
    results, AR1_ROOT,
    power_scenario_plot(
      AR1_N, 1.4, function(size) simulate_ar1_noise(size, AR1_RHO)
    ),
    selected_jump = 1.4,
    plot_results = set_algorithm_order(results, ar1_method_labels(calibration))
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_ar1()
}
