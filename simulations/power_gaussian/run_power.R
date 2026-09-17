## Gaussian power study: SVP paper Figures 1, 2, 6, and 8, plus two metrics.

library(svpChange)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_ROOT <- file.path("simulations", "power_gaussian")

fit_gaussian_methods <- function(
    y, constants = read_gaussian_calibration()) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
  )
  methods <- list(
    normalise_boundaries(changepoint::cpts(pelt), n),
    normalise_boundaries(
      SVP(y, 2 * log(n), "gaussian_mean", subtests = "right")$changepoints,
      n
    ),
    normalise_boundaries(
      SVP(y, constants[["right"]] * log(n), "gaussian_mean",
          subtests = "right")$changepoints,
      n
    ),
    normalise_boundaries(
      SVP(y, constants[["both"]] * log(n), "gaussian_mean",
          subtests = "both")$changepoints,
      n
    )
  )
  names(methods) <- c(
    "PELT",
    gaussian_method_label("right", 2),
    gaussian_method_label("right", constants[["right"]]),
    gaussian_method_label("both", constants[["both"]])
  )
  methods
}

run_gaussian_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L,
    constants = read_gaussian_calibration()) {
  run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) {
      fit_gaussian_methods(y, constants)
    },
    tolerance = round(n * 0.0025),
    workers = workers,
    seed = seed
  )
}

run_and_save_gaussian <- function(
    workers = power_default_workers(), root = GAUSSIAN_ROOT) {
  constants <- read_gaussian_calibration()
  results <- run_gaussian_power(workers = workers, constants = constants)
  save_power_outputs(
    results, root,
    power_scenario_plot(1000L, 0.6, stats::rnorm),
    plot_results = set_algorithm_order(
      results,
      c(
        "PELT",
        gaussian_method_label("right", 2),
        gaussian_method_label("right", constants[["right"]]),
        gaussian_method_label("both", constants[["both"]])
      )
    )
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_gaussian()
}
