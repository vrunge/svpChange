## Re-time only the PELT rows in the fixed-n, varying-k timing studies.
## The data seeds and fitting functions are unchanged; the timing pass is
## isolated from the other algorithms and uses a longer measurement batch.

library(svpChange)
library(changepoint)

source(file.path("simulations", "time_common.R"))

time_one_fit <- function(data, fit_method, target_seconds = 0.05) {
  batch <- 1L
  repeat {
    started <- proc.time()[["elapsed"]]
    for (i in seq_len(batch)) {
      boundaries <- fit_method(data)
    }
    block_elapsed <- proc.time()[["elapsed"]] - started
    if (block_elapsed >= target_seconds || batch >= 100000L) break
    growth <- if (block_elapsed <= 0) {
      10L
    } else {
      max(2L, ceiling(1.1 * target_seconds / block_elapsed))
    }
    batch <- min(100000L, as.integer(batch * growth))
  }
  list(
    time = block_elapsed / batch,
    batch = batch,
    block_elapsed = block_elapsed,
    boundaries = boundaries
  )
}

rerun_gaussian_pelt_fixed_k <- function(
    results, k_values = 0:100, reps = 12L, seed = 456L,
    target_seconds = 0.05) {
  source(file.path("simulations", "gaussian_common.R"))
  source(file.path("simulations", "time_gaussian", "run_time.R"))
  rows <- vector("list", length(k_values) * reps)
  index <- 0L
  for (k in k_values) {
    for (rep in seq_len(reps)) {
      set.seed(seed + 100000L + k * 100L + rep)
      data <- simulate_gaussian_time_data(10000L, k, jump = 10)
      gc(FALSE)
      timed <- time_one_fit(
        data,
        fit_method = function(y) fit_gaussian_time_method(y, "PELT"),
        target_seconds = target_seconds
      )
      index <- index + 1L
      rows[[index]] <- data.frame(
        experiment = "vary_k", n = 10000L, k = k, rep = rep,
        method = "PELT", time = timed$time, batch = timed$batch,
        block_elapsed = timed$block_elapsed,
        detected = length(timed$boundaries)
      )
    }
  }
  replacement <- dplyr::bind_rows(rows)
  old <- results[!(results$experiment == "vary_k" &
                     results$method == "PELT" &
                     results$k %in% k_values & results$rep <= reps), ,
                 drop = FALSE]
  dplyr::bind_rows(old, replacement)
}

rerun_ar1_pelt_inflated_fixed_k <- function(
    results, k_values = 0:40, reps = 64L, seed = 810L,
    target_seconds = 0.05) {
  source(file.path("simulations", "power_ar1", "calibrate_ar1.R"))
  source(file.path("simulations", "time_ar1", "run_time.R"))
  tasks <- expand.grid(
    k = k_values, rep = seq_len(reps), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  set.seed(seed + 2000000L)
  tasks <- tasks[sample.int(nrow(tasks)), , drop = FALSE]
  rows <- vector("list", nrow(tasks))
  index <- 0L
  for (task in seq_len(nrow(tasks))) {
    k <- tasks$k[task]
    rep <- tasks$rep[task]
    set.seed(seed + 100000L + k * 100L + rep)
    data <- simulate_ar1_time_data(2000L, k, rho = AR1_RHO)
    gc(FALSE)
    timed <- time_one_fit(
      data,
      fit_method = function(y) {
        pelt_inflated_boundaries(y, AR1_RHO)
      },
      target_seconds = target_seconds
    )
    index <- index + 1L
    rows[[index]] <- data.frame(
      experiment = "vary_k", n = 2000L, k = k, rep = rep,
      method = "PELT inflated", time = timed$time, batch = timed$batch,
      block_elapsed = timed$block_elapsed,
      detected = length(timed$boundaries)
    )
  }
  replacement <- dplyr::bind_rows(rows)
  old <- results[!(results$experiment == "vary_k" &
                     results$method == "PELT inflated" &
                     results$k %in% k_values & results$rep <= reps), ,
                 drop = FALSE]
  dplyr::bind_rows(old, replacement)
}

save_rerun_ar1_outputs <- function(target_seconds = 0.05) {
  ar1_path <- file.path("simulations", "time_ar1", "ar1_time_results.rds")
  ar1 <- readRDS(ar1_path)
  ar1 <- rerun_ar1_pelt_inflated_fixed_k(
    ar1, target_seconds = target_seconds
  )
  saveRDS(ar1, ar1_path)
  utils::write.csv(
    ar1, file.path("simulations", "time_ar1", "ar1_time_results.csv"),
    row.names = FALSE
  )
  source(file.path("simulations", "time_ar1", "run_time.R"))
  summaries <- plot_paper_time_study(
    ar1,
    root = "simulations/time_ar1",
    paper_root = normalizePath(file.path("..", "SVP_NEW_Figures2"),
                               mustWork = FALSE),
    prefix = "ar1_time",
    paper_prefix = "2_ar1",
    slope_min_n = 1000L,
    slope_label_factors = c(12.0, 6.0, 3.0, 1.5)
  )
  utils::write.csv(
    summaries$summary_n,
    file.path("simulations", "time_ar1", "ar1_time_summary_by_n.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$summary_detected,
    file.path("simulations", "time_ar1", "ar1_time_summary_by_detected.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$slopes,
    file.path("simulations", "time_ar1", "ar1_time_slopes.csv"),
    row.names = FALSE
  )
  invisible(ar1)
}

save_rerun_outputs <- function() {
  gaussian_path <- file.path(
    "simulations", "time_gaussian", "gaussian_time_results.rds"
  )
  gaussian <- readRDS(gaussian_path)
  gaussian <- rerun_gaussian_pelt_fixed_k(gaussian)
  saveRDS(gaussian, gaussian_path)
  utils::write.csv(
    gaussian,
    file.path("simulations", "time_gaussian", "gaussian_time_results.csv"),
    row.names = FALSE
  )
  source(file.path("simulations", "time_gaussian", "run_time.R"))
  plot_gaussian_paper_time(
    gaussian,
    normalizePath(file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE),
    "gaussian_time"
  )

  ar1_path <- file.path("simulations", "time_ar1", "ar1_time_results.rds")
  ar1 <- readRDS(ar1_path)
  ar1 <- rerun_ar1_pelt_inflated_fixed_k(ar1)
  saveRDS(ar1, ar1_path)
  utils::write.csv(
    ar1, file.path("simulations", "time_ar1", "ar1_time_results.csv"),
    row.names = FALSE
  )
  source(file.path("simulations", "time_ar1", "run_time.R"))
  calibration <- readRDS(file.path(
    "simulations", "power_ar1", "ar1_calibration.rds"
  ))
  summaries <- plot_paper_time_study(
    ar1,
    root = "simulations/time_ar1",
    paper_root = normalizePath(file.path("..", "SVP_NEW_Figures2"),
                               mustWork = FALSE),
    prefix = "ar1_time",
    paper_prefix = "2_ar1",
    slope_min_n = 1000L,
    slope_label_factors = c(12.0, 6.0, 3.0, 1.5)
  )
  utils::write.csv(
    summaries$summary_n,
    file.path("simulations", "time_ar1", "ar1_time_summary_by_n.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$summary_detected,
    file.path("simulations", "time_ar1", "ar1_time_summary_by_detected.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$slopes,
    file.path("simulations", "time_ar1", "ar1_time_slopes.csv"),
    row.names = FALSE
  )
  invisible(list(gaussian = gaussian, ar1 = ar1))
}

if (identical(tolower(Sys.getenv("SVP_RERUN_PELT_TIMING")), "true")) {
  save_rerun_outputs()
}
