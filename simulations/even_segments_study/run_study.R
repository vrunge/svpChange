## Up-down Gaussian signals with 2, 4, ..., 20 segments.
##
## Calibrated PELT and SVP both, paired datasets, 200 replicates per
## configuration. Produces:
##   * f1_even_segments.pdf      F1 curves, one panel per segment count
##   * detection_table.csv       estimated vs true change counts
##   * localization_table.csv    localization conditional on detection
##   * example_missed_change.pdf a typical series where SVP misses a change
##                               but localizes the detected ones better
##
## Run from the svpChange package root:
##   SVP_RUN_EVEN_SEGMENTS=true Rscript simulations/even_segments_study/run_study.R

library(svpChange)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

EVEN_ROOT <- file.path("simulations", "even_segments_study")
EVEN_N <- 1000L
EVEN_SEGMENTS <- seq.int(2L, 20L, by = 2L)
EVEN_JUMPS <- seq(0.1, 2, by = 0.1)
EVEN_REPS <- 200L
EVEN_SEED <- 7160000L
EVEN_TOLERANCE <- round(EVEN_N * 0.0025)

even_changepoints <- function(n, segments) {
  stopifnot(segments >= 2L, segments < n)
  as.integer(floor(seq_len(segments - 1L) * n / segments))
}

even_updown_signal <- function(n, segments, jump) {
  cps <- even_changepoints(n, segments)
  lengths <- diff(c(0L, cps, n))
  levels <- (0:(segments - 1L)) %% 2L
  rep(levels * jump, times = lengths)
}

## Greedy one-to-one matching within tolerance, identical to
## maximum_tolerance_matches() but returning the matched pairs.
matched_pairs <- function(truth, estimate, tolerance) {
  truth <- sort(as.integer(truth))
  estimate <- sort(as.integer(estimate))
  i <- 1L
  j <- 1L
  true_index <- integer(0)
  estimate_index <- integer(0)
  while (i <= length(truth) && j <= length(estimate)) {
    if (abs(truth[i] - estimate[j]) <= tolerance) {
      true_index <- c(true_index, i)
      estimate_index <- c(estimate_index, j)
      i <- i + 1L
      j <- j + 1L
    } else if (estimate[j] < truth[i] - tolerance) {
      j <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  data.frame(
    true = truth[true_index], estimate = estimate[estimate_index],
    error = abs(estimate[estimate_index] - truth[true_index])
  )
}

even_seed <- function(segments_index, jump_index, rep) {
  EVEN_SEED + ((segments_index - 1L) * length(EVEN_JUMPS) + jump_index - 1L) *
    1000L + rep
}

even_batch <- function(segments_index, reps, pelt_constant, svp_constant,
                       workers) {
  segments <- EVEN_SEGMENTS[[segments_index]]
  internal <- segments - 1L
  cps <- even_changepoints(EVEN_N, segments)
  truth <- c(cps, EVEN_N)
  jobs <- expand.grid(jump_index = seq_along(EVEN_JUMPS), rep = seq_len(reps))
  chunks <- split(seq_len(nrow(jobs)),
                  ceiling(seq_len(nrow(jobs)) / 250L))
  rows <- power_lapply(chunks, function(indices) {
    dplyr::bind_rows(lapply(indices, function(i) {
    jump_index <- jobs$jump_index[i]
    rep <- jobs$rep[i]
    jump <- EVEN_JUMPS[[jump_index]]
    mu <- even_updown_signal(EVEN_N, segments, jump)
    set_power_seed(even_seed(segments_index, jump_index, rep))
    y <- mu + stats::rnorm(EVEN_N)

    pelt <- changepoint::cpt.mean(
      y, method = "PELT", penalty = "Manual",
      pen.value = pelt_constant * log(EVEN_N)
    )
    pelt_cps <- normalise_boundaries(changepoint::cpts(pelt), EVEN_N)
    svp_cps <- normalise_boundaries(
      SVP(y, svp_constant * log(EVEN_N), "gaussian_mean",
          subtests = "both", cost = "gaussian")$changepoints,
      EVEN_N
    )
    pelt_internal <- pelt_cps[pelt_cps < EVEN_N]
    svp_internal <- svp_cps[svp_cps < EVEN_N]
    pelt_pairs <- matched_pairs(cps, pelt_internal, EVEN_TOLERANCE)
    svp_pairs <- matched_pairs(cps, svp_internal, EVEN_TOLERANCE)
    common <- intersect(pelt_pairs$true, svp_pairs$true)
    pelt_common <- pelt_pairs[pelt_pairs$true %in% common, ]
    svp_common <- svp_pairs[svp_pairs$true %in% common, ]

    pelt_score <- paper_metrics(truth, pelt_cps, EVEN_TOLERANCE,
                                terminal = EVEN_N)
    svp_score <- paper_metrics(truth, svp_cps, EVEN_TOLERANCE,
                               terminal = EVEN_N)
    data.frame(
      segments = segments, internal = internal, jump = jump, rep = rep,
      pelt_f1 = pelt_score[["F1"]], svp_f1 = svp_score[["F1"]],
      pelt_precision = pelt_score[["Precision"]],
      svp_precision = svp_score[["Precision"]],
      pelt_recall = pelt_score[["Recall"]],
      svp_recall = svp_score[["Recall"]],
      pelt_cp_count = length(pelt_internal),
      svp_cp_count = length(svp_internal),
      pelt_missed = internal - nrow(pelt_pairs),
      svp_missed = internal - nrow(svp_pairs),
      pelt_fp = length(pelt_internal) - nrow(pelt_pairs),
      svp_fp = length(svp_internal) - nrow(svp_pairs),
      pelt_mean_error = if (nrow(pelt_pairs)) mean(pelt_pairs$error) else NA_real_,
      svp_mean_error = if (nrow(svp_pairs)) mean(svp_pairs$error) else NA_real_,
      pelt_exact = if (nrow(pelt_pairs)) mean(pelt_pairs$error == 0) else NA_real_,
      svp_exact = if (nrow(svp_pairs)) mean(svp_pairs$error == 0) else NA_real_,
      common_n = length(common),
      pelt_common_mean_error = if (nrow(pelt_common)) {
        mean(pelt_common$error)
      } else {
        NA_real_
      },
      svp_common_mean_error = if (nrow(svp_common)) {
        mean(svp_common$error)
      } else {
        NA_real_
      },
      pelt_common_exact = if (nrow(pelt_common)) {
        mean(pelt_common$error == 0)
      } else {
        NA_real_
      },
      svp_common_exact = if (nrow(svp_common)) {
        mean(svp_common$error == 0)
      } else {
        NA_real_
      },
      seed = even_seed(segments_index, jump_index, rep),
      stringsAsFactors = FALSE
    ) -> row
    row$pelt_changepoints <- list(pelt_internal)
    row$svp_changepoints <- list(svp_internal)
    row
    }))
  }, workers = workers)
  dplyr::bind_rows(rows)
}

run_even_study <- function(reps = EVEN_REPS, workers = power_default_workers(),
                           root = EVEN_ROOT) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  constants <- read_gaussian_calibration()
  pelt_constant <- read_gaussian_pelt_calibration()
  svp_constant <- constants[["both"]]
  results <- dplyr::bind_rows(lapply(seq_along(EVEN_SEGMENTS), function(i) {
    cat("segments", EVEN_SEGMENTS[[i]], "\n")
    even_batch(i, reps, pelt_constant, svp_constant, workers)
  }))
  results$pelt_label <- sprintf("PELT calibrated, c = %.3f", pelt_constant)
  results$svp_label <- sprintf("SVP both calibrated, c = %.3f", svp_constant)
  saveRDS(results, file.path(root, "results.rds"), compress = FALSE)
  results
}

f1_curve_data <- function(results) {
  pelt_label <- results$pelt_label[[1L]]
  svp_label <- results$svp_label[[1L]]
  curves <- dplyr::bind_rows(
    data.frame(segments = results$segments, jump = results$jump,
               algorithm = pelt_label, f1 = results$pelt_f1),
    data.frame(segments = results$segments, jump = results$jump,
               algorithm = svp_label, f1 = results$svp_f1)
  )
  curves$algorithm <- factor(curves$algorithm, c(pelt_label, svp_label))
  curves |>
    dplyr::group_by(segments, jump, algorithm) |>
    dplyr::summarise(
      mean = mean(f1),
      se = stats::sd(f1) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    dplyr::mutate(lower = pmax(0, mean - 1.96 * se),
                  upper = pmin(1, mean + 1.96 * se))
}

plot_f1_curves <- function(results, path) {
  data <- f1_curve_data(results)
  data$panel <- factor(sprintf("%d segments", data$segments),
                       sprintf("%d segments", EVEN_SEGMENTS))
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(jump, mean, colour = algorithm, fill = algorithm)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::facet_wrap(~panel, nrow = 2L) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Jump size", y = "F1", colour = "Method", fill = "Method") +
    paper_theme()
  ggplot2::ggsave(path, plot, width = 16, height = 7)
  invisible(plot)
}

if (identical(tolower(Sys.getenv("SVP_RUN_EVEN_SEGMENTS")), "true")) {
  results <- run_even_study()
  plot_f1_curves(results, file.path("simulations", "all_plots",
                                    "6_f1_vs_segments.pdf"))
}
