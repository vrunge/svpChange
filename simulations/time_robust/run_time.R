## Runtime scaling under Student-t(2) noise.

library(svpChange)
library(changepoint)
library(robseg)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "time_common.R"))
source(file.path("simulations", "power_robust", "calibrate_robust.R"))

ROBUST_TIME_ROOT <- file.path("simulations", "time_robust")
ROBUST_TIME_MAX_N <- 4000L
format_robust_time_parameter <- function(x) {
  trimws(formatC(x, digits = 6L, format = "fg"))
}

ROBUST_TIME_METHODS <- function(calibration) {
  selected <- calibration$selected
  c(
    "PELT", "RFPOP",
    paste0(
      "SVP MedianMood / right / alpha = ",
      format_robust_time_parameter(selected$mood_alpha)
    ),
    paste0(
      "SVP Wilcoxon / right / c = ",
      format_robust_time_parameter(selected$wilcoxon_constant)
    )
  )
}

simulate_robust_time_data <- function(n, changes = 0L, jump = 4) {
  sizes <- balanced_segment_sizes(n, changes)
  ts_generator(
    chpts = cumsum(sizes),
    parameters = rep(c(0, jump), length.out = changes + 1L),
    type = "student", df = 2, scale = 1
  )
}

fit_robust_time_method <- function(y, method, changes, calibration) {
  n <- length(y)
  selected <- calibration$selected
  methods <- ROBUST_TIME_METHODS(calibration)
  if (method == "PELT") {
    fit <- changepoint::cpt.mean(
      y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
    )
    return(normalise_boundaries(changepoint::cpts(fit), n))
  }
  if (method == "RFPOP") {
    return(fit_rfpop(y, selected$rfpop_constant)$boundaries)
  }
  oracle_segments <- changes + 1L
  if (identical(method, methods[[3L]])) {
    fit <- SVP(
      y, mood_threshold(n, oracle_segments, selected$mood_alpha),
      "MedianMoodCost", subtests = "right"
    )
    return(normalise_boundaries(fit$changepoints, n))
  }
  if (identical(method, methods[[4L]])) {
    fit <- SVP(
      y, wilcoxon_threshold(n, oracle_segments,
                            selected$wilcoxon_constant),
      "WilcoxonCost", subtests = "right"
    )
    return(normalise_boundaries(fit$changepoints, n))
  }
  stop("Unknown robust timing method: ", method)
}

run_time_robust <- function(
    n_values = round(2^(seq(8, log2(ROBUST_TIME_MAX_N), length.out = 16))),
    k_values = seq(0L, 40L, by = 2L),
    fixed_n = 2000L,
    reps_n = 64L,
    reps_k = 64L,
    seed = 910L,
    min_batch_time = 0.05,
    randomize_method_order = FALSE,
    calibration = readRDS(file.path(
      "simulations", "power_robust", "robust_calibration.rds"
    ))) {
  if (any(n_values > ROBUST_TIME_MAX_N) || fixed_n > ROBUST_TIME_MAX_N) {
    stop("the robust time study is limited to n <= ", ROBUST_TIME_MAX_N)
  }
  methods <- ROBUST_TIME_METHODS(calibration)
  set.seed(seed - 1L)
  warm_data <- simulate_robust_time_data(256L, 0L)
  invisible(lapply(methods, function(method) {
    fit_robust_time_method(warm_data, method, 0L, calibration)
  }))
  benchmark_time_study(
    n_values, k_values, fixed_n, reps_n, reps_k,
    methods = methods,
    simulate_data = simulate_robust_time_data,
    fit_method = function(y, method, changes) {
      fit_robust_time_method(y, method, changes, calibration)
    },
    seed = seed,
    min_batch_time = min_batch_time,
    randomize_method_order = randomize_method_order
  )
}

run_and_save_time_robust <- function(
    workers = 1L, calibration = readRDS(file.path(
      "simulations", "power_robust", "robust_calibration.rds"
    ))) {
  started <- proc.time()[["elapsed"]]
  results <- run_time_robust(calibration = calibration)
  run_elapsed <- proc.time()[["elapsed"]] - started
  dir.create(ROBUST_TIME_ROOT, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(ROBUST_TIME_ROOT, "robust_time_results.rds"))
  utils::write.csv(
    results, file.path(ROBUST_TIME_ROOT, "robust_time_results.csv"),
    row.names = FALSE
  )
  summaries <- plot_paper_time_study(
    results,
    root = ROBUST_TIME_ROOT,
    paper_root = normalizePath(
      file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
    ),
    prefix = "robust_time",
    paper_prefix = "2_robust",
    slope_min_n = 500L,
    detected_methods = setdiff(ROBUST_TIME_METHODS(calibration), "PELT")
  )
  utils::write.csv(
    summaries$summary_n,
    file.path(ROBUST_TIME_ROOT, "robust_time_summary_by_n.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$summary_detected,
    file.path(ROBUST_TIME_ROOT, "robust_time_summary_by_detected.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$slopes,
    file.path(ROBUST_TIME_ROOT, "robust_time_slopes.csv"),
    row.names = FALSE
  )
  manifest <- c(
    paste0("timestamp = ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("svpChange = ", as.character(utils::packageVersion("svpChange"))),
    paste0("robseg = ", as.character(utils::packageVersion("robseg"))),
    paste0("changepoint = ", as.character(utils::packageVersion("changepoint"))),
    paste0("R = ", R.version.string),
    paste0("platform = ", R.version$platform),
    paste0("git_sha = ", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    "noise = Student-t(df = 2, scale = 1)",
    "n_grid = 256 to 4000 (16 values); fixed_n = 2000; k_grid = 0:2:40",
    "replicates = 64 for each method and experiment",
    "timer = proc.time elapsed; adaptive batching target 0.05 s; fixed method order",
    paste0("simulation_elapsed_seconds = ", sprintf("%.3f", run_elapsed)),
    paste0("selected_rfpop_c = ", calibration$selected$rfpop_constant),
    paste0("selected_mood_alpha = ", calibration$selected$mood_alpha),
    paste0("selected_wilcoxon_c = ", calibration$selected$wilcoxon_constant),
    "SVP oracle_segments = simulated changes + 1"
  )
  writeLines(
    manifest,
    file.path(ROBUST_TIME_ROOT, "robust_time_manifest.txt")
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_time_robust()
}
