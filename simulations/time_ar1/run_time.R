## Runtime scaling under stationary Gaussian AR(1) noise.

library(svpChange)
library(changepoint)
library(DeCAFS)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "time_common.R"))
source(file.path("simulations", "power_ar1", "calibrate_ar1.R"))

AR1_TIME_ROOT <- file.path("simulations", "time_ar1")
AR1_TIME_METHODS <- function(calibration) ar1_method_labels(calibration)

simulate_ar1_time_data <- function(
    n, changes = 0L, rho = AR1_RHO, jump = 10) {
  sizes <- balanced_segment_sizes(n, changes)
  ts_generator(
    chpts = cumsum(sizes),
    parameters = rep(c(0, jump), length.out = changes + 1L),
    sd_noise = sqrt(1 - rho^2),
    rho = rho,
    type = "gaussAR1"
  )
}

fit_ar1_time_method <- function(y, method, calibration, rho = AR1_RHO) {
  n <- length(y)
  selected <- calibration$selected
  if (identical(method, "PELT AR1 approximate")) {
    return(pelt_ar1_approximate_boundaries(
      y, rho, penalty = selected$approximate_constant * log(n)
    ))
  }
  if (identical(method, "PELT inflated")) {
    return(pelt_inflated_boundaries(y, rho))
  }
  if (identical(method, "DeCAFS AR1")) {
    return(decafs_ar1_boundaries(
      y, rho, penalty = selected$decafs_constant * log(n)
    ))
  }
  svp_label <- AR1_TIME_METHODS(calibration)[[4L]]
  if (identical(method, svp_label)) {
    return(svp_ar1focus_boundaries(
      y, rho, constant = selected$svp_constant, cost = "ar1"
    ))
  }
  stop("Unknown AR(1) timing method: ", method)
}

run_time_ar1 <- function(
    n_values = round(2^(seq(8, 15, length.out = 24))),
    k_values = 0:99,
    fixed_n = 2000L,
    reps_n = 64L,
    reps_k = 64L,
    seed = 810L,
    min_batch_time = 0.05,
    rho = AR1_RHO,
    jump = 10,
    calibration = readRDS(file.path(
      "simulations", "power_ar1", "ar1_calibration.rds"
    ))) {
  require_decafs()
  methods <- AR1_TIME_METHODS(calibration)
  set.seed(seed - 1L)
  warm_data <- simulate_ar1_time_data(256L, 0L, rho = rho)
  invisible(lapply(methods, function(method) {
    fit_ar1_time_method(warm_data, method, calibration, rho)
  }))
  benchmark_time_study(
    n_values, k_values, fixed_n, reps_n, reps_k,
    methods = methods,
    simulate_data = function(n, changes) {
      simulate_ar1_time_data(n, changes, rho = rho, jump = jump)
    },
    fit_method = function(y, method, changes) {
      fit_ar1_time_method(y, method, calibration, rho)
    },
    seed = seed,
    min_batch_time = min_batch_time,
    randomize_method_order = TRUE
  )
}

run_and_save_time_ar1_detected <- function(
    existing_results = file.path(
      AR1_TIME_ROOT, "ar1_time_results.rds"
    ),
    k_values = 0:99,
    fixed_n = 2000L,
    reps_k = 64L,
    seed = 810L,
    min_batch_time = 0.05,
    rho = AR1_RHO,
    jump = 10,
    calibration = readRDS(file.path(
      "simulations", "power_ar1", "ar1_calibration.rds"
    ))) {
  require_decafs()
  started <- proc.time()[["elapsed"]]
  methods <- AR1_TIME_METHODS(calibration)
  set.seed(seed - 1L)
  warm_data <- simulate_ar1_time_data(256L, 0L, rho = rho, jump = jump)
  invisible(lapply(methods, function(method) {
    fit_ar1_time_method(warm_data, method, calibration, rho)
  }))
  detected_results <- benchmark_time_study(
    n_values = integer(), k_values = k_values, fixed_n = fixed_n,
    reps_n = 0L, reps_k = reps_k, methods = methods,
    simulate_data = function(n, changes) {
      simulate_ar1_time_data(n, changes, rho = rho, jump = jump)
    },
    fit_method = function(y, method, changes) {
      fit_ar1_time_method(y, method, calibration, rho)
    },
    seed = seed, min_batch_time = min_batch_time,
    randomize_method_order = TRUE
  )
  if (file.exists(existing_results)) {
    previous <- readRDS(existing_results)
    previous <- dplyr::filter(previous, experiment != "vary_k")
    results <- dplyr::bind_rows(previous, detected_results)
  } else {
    results <- detected_results
  }
  run_elapsed <- proc.time()[["elapsed"]] - started
  dir.create(AR1_TIME_ROOT, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, existing_results)
  utils::write.csv(
    results, file.path(AR1_TIME_ROOT, "ar1_time_results.csv"),
    row.names = FALSE
  )
  summaries <- plot_paper_time_study(
    results,
    root = AR1_TIME_ROOT,
    paper_root = normalizePath(
      file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
    ),
    prefix = "ar1_time",
    paper_prefix = "2_ar1",
    slope_min_n = 1000L,
    detected_max = 100L,
    slope_label_factors = c(12.0, 6.0, 3.0, 1.5)
  )
  utils::write.csv(
    summaries$summary_n,
    file.path(AR1_TIME_ROOT, "ar1_time_summary_by_n.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$summary_detected,
    file.path(AR1_TIME_ROOT, "ar1_time_summary_by_detected.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$slopes,
    file.path(AR1_TIME_ROOT, "ar1_time_slopes.csv"),
    row.names = FALSE
  )
  manifest <- c(
    paste0("timestamp = ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("svpChange = ", as.character(utils::packageVersion("svpChange"))),
    paste0("DeCAFS = ", as.character(utils::packageVersion("DeCAFS"))),
    paste0("changepoint = ", as.character(utils::packageVersion("changepoint"))),
    paste0("R = ", R.version.string),
    paste0("platform = ", R.version$platform),
    paste0("git_sha = ", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste0("noise = stationary Gaussian AR(1), rho = ", rho,
           ", marginal variance = 1"),
    "vary_n results reused from the previous AR(1) timing study",
    "vary_k = 0:99; fixed_n = 2000; jump = 10",
    paste0("replicates = ", reps_k, " for each method and k"),
    "timer = proc.time elapsed; adaptive batching target 0.05 s",
    "displayed detected-segment range = 1:100",
    paste0("simulation_elapsed_seconds = ", sprintf("%.3f", run_elapsed)),
    paste0("selected_approximate_c = ", calibration$selected$approximate_constant),
    paste0("selected_decafs_c = ", calibration$selected$decafs_constant),
    paste0("selected_svp_c = ", calibration$selected$svp_constant),
    "PELT inflated = 3 * log(n) * (1 + rho) / (1 - rho)"
  )
  writeLines(
    manifest,
    file.path(AR1_TIME_ROOT, "ar1_time_manifest.txt")
  )
  invisible(results)
}

run_and_save_time_ar1_pelt_detected <- function(
    existing_results = file.path(
      AR1_TIME_ROOT, "ar1_time_results.rds"
    ),
    k_values = 0:99,
    fixed_n = 2000L,
    reps_k = 16L,
    seed = 1810L,
    min_batch_time = 0.2,
    rho = AR1_RHO,
    jump = 10,
    calibration = readRDS(file.path(
      "simulations", "power_ar1", "ar1_calibration.rds"
    ))) {
  require_decafs()
  started <- proc.time()[["elapsed"]]
  methods <- AR1_TIME_METHODS(calibration)[1:2]
  set.seed(seed - 1L)
  warm_data <- simulate_ar1_time_data(256L, 0L, rho = rho, jump = jump)
  invisible(lapply(methods, function(method) {
    fit_ar1_time_method(warm_data, method, calibration, rho)
  }))
  pelt_results <- benchmark_time_study(
    n_values = integer(), k_values = k_values, fixed_n = fixed_n,
    reps_n = 0L, reps_k = reps_k, methods = methods,
    simulate_data = function(n, changes) {
      simulate_ar1_time_data(n, changes, rho = rho, jump = jump)
    },
    fit_method = function(y, method, changes) {
      fit_ar1_time_method(y, method, calibration, rho)
    },
    seed = seed, min_batch_time = min_batch_time,
    randomize_method_order = TRUE
  )
  if (!file.exists(existing_results)) {
    stop("Existing AR(1) timing results are required for a PELT-only rerun")
  }
  previous <- readRDS(existing_results)
  previous <- dplyr::filter(
    previous,
    !(experiment == "vary_k" & method %in% methods)
  )
  results <- dplyr::bind_rows(previous, pelt_results)
  run_elapsed <- proc.time()[["elapsed"]] - started
  dir.create(AR1_TIME_ROOT, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, existing_results)
  utils::write.csv(
    results, file.path(AR1_TIME_ROOT, "ar1_time_results.csv"),
    row.names = FALSE
  )
  summaries <- plot_paper_time_study(
    results,
    root = AR1_TIME_ROOT,
    paper_root = normalizePath(
      file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
    ),
    prefix = "ar1_time",
    paper_prefix = "2_ar1",
    slope_min_n = 1000L,
    detected_max = 100L,
    slope_label_factors = c(12.0, 6.0, 3.0, 1.5)
  )
  utils::write.csv(
    summaries$summary_n,
    file.path(AR1_TIME_ROOT, "ar1_time_summary_by_n.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$summary_detected,
    file.path(AR1_TIME_ROOT, "ar1_time_summary_by_detected.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$slopes,
    file.path(AR1_TIME_ROOT, "ar1_time_slopes.csv"),
    row.names = FALSE
  )
  manifest <- c(
    paste0("timestamp = ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("svpChange = ", as.character(utils::packageVersion("svpChange"))),
    paste0("DeCAFS = ", as.character(utils::packageVersion("DeCAFS"))),
    paste0("changepoint = ", as.character(utils::packageVersion("changepoint"))),
    paste0("R = ", R.version.string),
    paste0("platform = ", R.version$platform),
    paste0("git_sha = ", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste0("noise = stationary Gaussian AR(1), rho = ", rho,
           ", marginal variance = 1"),
    "vary_n results and DeCAFS/SVP vary_k results reused",
    "PELT vary_k = 0:99; fixed_n = 2000; jump = 10",
    paste0("PELT replicates = ", reps_k, " for each method and k"),
    paste0("PELT timer batch target = ", min_batch_time, " s"),
    "displayed detected-segment range = 1:100",
    paste0("pelt_simulation_elapsed_seconds = ", sprintf("%.3f", run_elapsed)),
    paste0("selected_approximate_c = ", calibration$selected$approximate_constant),
    paste0("selected_decafs_c = ", calibration$selected$decafs_constant),
    paste0("selected_svp_c = ", calibration$selected$svp_constant),
    "PELT inflated = 3 * log(n) * (1 + rho) / (1 - rho)"
  )
  writeLines(
    manifest,
    file.path(AR1_TIME_ROOT, "ar1_time_manifest.txt")
  )
  invisible(results)
}

run_and_save_time_ar1 <- function(
    calibration = readRDS(file.path(
      "simulations", "power_ar1", "ar1_calibration.rds"
    ))) {
  require_decafs()
  started <- proc.time()[["elapsed"]]
  results <- run_time_ar1(calibration = calibration)
  run_elapsed <- proc.time()[["elapsed"]] - started
  dir.create(AR1_TIME_ROOT, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(AR1_TIME_ROOT, "ar1_time_results.rds"))
  utils::write.csv(
    results, file.path(AR1_TIME_ROOT, "ar1_time_results.csv"),
    row.names = FALSE
  )
  summaries <- plot_paper_time_study(
    results,
    root = AR1_TIME_ROOT,
    paper_root = normalizePath(
      file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
    ),
    prefix = "ar1_time",
    paper_prefix = "2_ar1",
    slope_min_n = 1000L,
    detected_max = 100L,
    slope_label_factors = c(12.0, 6.0, 3.0, 1.5)
  )
  utils::write.csv(
    summaries$summary_n,
    file.path(AR1_TIME_ROOT, "ar1_time_summary_by_n.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$summary_detected,
    file.path(AR1_TIME_ROOT, "ar1_time_summary_by_detected.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$slopes,
    file.path(AR1_TIME_ROOT, "ar1_time_slopes.csv"),
    row.names = FALSE
  )
  manifest <- c(
    paste0("timestamp = ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("svpChange = ", as.character(utils::packageVersion("svpChange"))),
    paste0("DeCAFS = ", as.character(utils::packageVersion("DeCAFS"))),
    paste0("changepoint = ", as.character(utils::packageVersion("changepoint"))),
    paste0("R = ", R.version.string),
    paste0("platform = ", R.version$platform),
    paste0("git_sha = ", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste0("noise = stationary Gaussian AR(1), rho = ",
           calibration$selected$rho,
           ", marginal variance = 1"),
    "n_grid = 256 to 32768 (24 values); fixed_n = 2000; k_grid = 0:99; jump = 10",
    "replicates = 64 for each method and experiment",
    "timer = proc.time elapsed; adaptive batching target 0.05 s",
    paste0("simulation_elapsed_seconds = ", sprintf("%.3f", run_elapsed)),
    paste0("selected_approximate_c = ", calibration$selected$approximate_constant),
    paste0("selected_decafs_c = ", calibration$selected$decafs_constant),
    paste0("selected_svp_c = ", calibration$selected$svp_constant),
    "PELT inflated = 3 * log(n) * (1 + rho) / (1 - rho)"
  )
  writeLines(
    manifest,
    file.path(AR1_TIME_ROOT, "ar1_time_manifest.txt")
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_time_ar1()
}
