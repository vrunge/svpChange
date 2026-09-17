## Extended Gaussian power study for the current inclusive-left implementation.
## The canonical Section 5.1 study remains in run_power.R; this file evaluates
## a separately calibrated subtests = "both" threshold.

source_path <- file.path("simulations", "power_gaussian", "run_power.R")
old_flag <- Sys.getenv("SVP_RUN_SIMULATIONS", unset = NA_character_)
Sys.setenv(SVP_RUN_SIMULATIONS = "false")
source(source_path)
if (is.na(old_flag)) Sys.unsetenv("SVP_RUN_SIMULATIONS") else
  Sys.setenv(SVP_RUN_SIMULATIONS = old_flag)

GAUSSIAN_BOTH_CONSTANT <- 1.9
GAUSSIAN_BOTH_ROOT <- file.path(
  "simulations", "power_gaussian", "reproduced_v020_both_calibrations"
)

fit_both_calibrated <- function(y, constant = GAUSSIAN_BOTH_CONSTANT) {
  n <- length(y)
  normalise_boundaries(
    SVP(y, constant * log(n), "gaussian_mean", subtests = "both")$changepoints,
    n
  )
}

run_both_calibrated_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L,
    constant = GAUSSIAN_BOTH_CONSTANT) {
  results <- run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) {
      setNames(
        list(fit_both_calibrated(y, constant)),
        paste0("SVP both c=", format(constant, trim = TRUE))
      )
    },
    tolerance = round(n * 0.0025),
    workers = workers,
    seed = seed
  )
  attr(results, "both_constant") <- constant
  results
}

run_both_comparison <- function(
    constants = c(1.9, 2.0), n = 1000L,
    jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L) {
  dplyr::bind_rows(lapply(constants, function(constant) {
    run_both_calibrated_power(
      n = n, jump_sizes = jump_sizes, reps = reps,
      workers = workers, seed = seed, constant = constant
    )
  }))
}

calibrate_both_null <- function(
    constants = seq(1.4, 2.4, 0.1), n = 1000L, reps = 100L,
    workers = power_default_workers(), seed = 700000L) {
  jobs <- expand.grid(
    multiplier = constants, rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  rows <- power_lapply(seq_len(nrow(jobs)), function(i) {
    task <- jobs[i, ]
    set_power_seed(seed + i)
    y <- stats::rnorm(n)
    boundaries <- fit_both_calibrated(y, task$multiplier)
    metric <- paper_metrics(n, boundaries, round(0.0025 * n))
    data.frame(
      multiplier = task$multiplier, rep = task$rep,
      F1 = metric[["F1"]],
      CorrectNumCP = metric[["CorrectNumCP"]],
      NumSegments = length(boundaries)
    )
  }, workers = workers)
  dplyr::bind_rows(rows)
}

save_both_null_calibration <- function(
    calibration = calibrate_both_null(),
    root = GAUSSIAN_BOTH_ROOT) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    calibration, file.path(root, "both_null_calibration.csv"),
    row.names = FALSE
  )
  summary <- calibration |>
    dplyr::group_by(multiplier) |>
    dplyr::summarise(
      F1 = mean(F1), CorrectNumCP = mean(CorrectNumCP),
      mean_segments = mean(NumSegments), .groups = "drop"
    )
  utils::write.csv(
    summary, file.path(root, "both_null_calibration_summary.csv"),
    row.names = FALSE
  )
  invisible(summary)
}

run_and_save_both_calibrated <- function(
    workers = power_default_workers(),
    constant = GAUSSIAN_BOTH_CONSTANT) {
  results <- run_both_calibrated_power(
    workers = workers, constant = constant
  )
  save_power_outputs(
    results, GAUSSIAN_BOTH_ROOT,
    power_scenario_plot(1000L, 0.6, stats::rnorm),
    plot_results = results
  )
  writeLines(
    c(
      "Current svpChange 0.2.0 source.",
      paste0("subtests = both, gamma = ", constant, " log(n)."),
      "Calibration target: independent no-change F1 at least 0.99.",
      "Design: n=1000, four patterns, jumps 0.1:2.0, 100 replicates, seed=123."
    ),
    file.path(GAUSSIAN_BOTH_ROOT, "README.txt")
  )
  invisible(results)
}

run_and_save_both_comparison <- function(
    workers = power_default_workers(),
    constants = c(1.9, 2.0), root = GAUSSIAN_BOTH_ROOT) {
  results <- run_both_comparison(workers = workers, constants = constants)
  save_power_outputs(
    results, root, power_scenario_plot(1000L, 0.6, stats::rnorm),
    plot_results = results
  )
  summary <- results |>
    dplyr::group_by(algorithm, pattern) |>
    dplyr::summarise(
      F1 = mean(F1), CorrectNumCP = mean(CorrectNumCP),
      MSE = mean(MSE), NumSegments = mean(NumSegments), .groups = "drop"
    )
  utils::write.csv(summary, file.path(root, "summary.csv"), row.names = FALSE)
  invisible(results)
}
