## Standalone review of robust SVP calibration on the well-log data.

WELL_LOG_ROOT <- file.path("simulations", "z_well_log")
WELL_LOG_PLOT_ROOT <- file.path(WELL_LOG_ROOT, "plots")
WELL_LOG_N_CALIBRATION <- 5000L
WELL_LOG_N_VALIDATION <- 5000L
WELL_LOG_N_POWER <- 500L
WELL_LOG_N_BOTH_NULL_STRESS <- 50L
WELL_LOG_NULL_TARGET <- 0.99
WELL_LOG_MATCH_TOLERANCE <- 10L
WELL_LOG_CALIBRATION_SEED <- 981000L
WELL_LOG_VALIDATION_SEED <- 982000L
WELL_LOG_POWER_SEED <- 983000L
WELL_LOG_BOTH_NULL_SEED <- 984000L

require_well_log_packages <- function() {
  packages <- c("svpChange", "robseg", "changepoint", "ggplot2", "Rcpp")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "))
  }
  if (!("cost" %in% names(formals(svpChange::SVP)))) {
    stop("Install the current svpChange source: SVP() must expose `cost`.")
  }
  cache <- file.path(tempdir(), "svp-well-log-sourcecpp")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(
    file.path(WELL_LOG_ROOT, "prefix_scan.cpp"),
    cacheDir = cache, rebuild = FALSE, showOutput = FALSE
  )
  invisible(TRUE)
}

parallel_map <- function(x, fun, workers) {
  if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(x, fun, mc.cores = workers, mc.preschedule = TRUE,
                       mc.set.seed = FALSE)
  } else {
    lapply(x, fun)
  }
}

internal_boundaries <- function(boundaries, n) {
  boundaries <- sort(unique(as.integer(boundaries)))
  boundaries[boundaries > 0L & boundaries < n]
}

match_changepoints <- function(reference, detected, tolerance) {
  reference <- sort(as.integer(reference))
  detected <- sort(as.integer(detected))
  if (!length(reference) && !length(detected)) {
    return(c(tp = 0, precision = 1, recall = 1, F1 = 1))
  }
  if (!length(reference) || !length(detected)) {
    return(c(tp = 0, precision = 0, recall = 0, F1 = 0))
  }
  i <- 1L
  j <- 1L
  tp <- 0L
  while (i <= length(reference) && j <= length(detected)) {
    if (abs(reference[i] - detected[j]) <= tolerance) {
      tp <- tp + 1L
      i <- i + 1L
      j <- j + 1L
    } else if (detected[j] < reference[i] - tolerance) {
      j <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  precision <- tp / length(detected)
  recall <- tp / length(reference)
  c(tp = tp, precision = precision, recall = recall,
    F1 = if (precision + recall > 0) {
      2 * precision * recall / (precision + recall)
    } else 0)
}

wilson_interval <- function(successes, trials, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  centre <- (p + z^2 / (2 * trials)) / denominator
  half <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) /
    denominator
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

fit_rfpop <- function(y, scale) {
  fit <- robseg::Rob_seg.std(
    y / scale, loss = "Outlier", lambda = 70, lthreshold = 2
  )
  list(
    fit = fit,
    changepoints = internal_boundaries(fit$t.est, length(y))
  )
}

fit_pelt <- function(y, scale) {
  fit <- changepoint::cpt.mean(
    y / scale, method = "PELT", penalty = "Manual", pen.value = 70
  )
  internal_boundaries(changepoint::cpts(fit), length(y))
}

fit_svp <- function(y, gamma, test, mode) {
  fit <- svpChange::SVP(
    y, gamma = gamma, test = test, subtests = mode, cost = "gaussian"
  )
  internal_boundaries(fit$changepoints, length(y))
}

segment_median_fit <- function(y, changepoints) {
  ends <- c(changepoints, length(y))
  starts <- c(1L, changepoints + 1L)
  medians <- vapply(seq_along(ends), function(i) {
    stats::median(y[starts[i]:ends[i]])
  }, numeric(1))
  list(
    fitted = rep(medians, ends - starts + 1L),
    medians = medians,
    starts = starts,
    ends = ends
  )
}

select_block_length <- function(residuals, run = 5L, max_lag = 200L) {
  values <- as.vector(stats::acf(
    residuals, lag.max = max_lag, plot = FALSE
  )$acf)[-1L]
  limit <- 1.96 / sqrt(length(residuals))
  candidates <- seq_len(length(values) - run + 1L)
  selected <- candidates[vapply(candidates, function(i) {
    all(abs(values[i:(i + run - 1L)]) < limit)
  }, logical(1))]
  if (!length(selected)) max_lag else selected[1L]
}

valid_block_starts <- function(n, block_length, changepoints) {
  starts <- seq_len(n - block_length + 1L)
  starts[!vapply(starts, function(s) {
    any(changepoints >= s & changepoints < s + block_length - 1L)
  }, logical(1))]
}

sample_blocks <- function(residuals, block_length, starts, seed) {
  set.seed(seed)
  n <- length(residuals)
  selected <- sample(starts, ceiling(n / block_length), replace = TRUE)
  values <- unlist(lapply(selected, function(s) {
    residuals[s:(s + block_length - 1L)]
  }), use.names = FALSE)
  values[seq_len(n)]
}

median_mood_statistic <- function(x) {
  n <- length(x)
  if (n < 2L) return(0)
  median_upper <- sort(x, partial = n %/% 2L + 1L)[n %/% 2L + 1L]
  below <- x < median_upper
  above <- x > median_upper
  total_below <- sum(below)
  total_above <- sum(above)
  effective_n <- as.numeric(total_below + total_above)
  if (effective_n == 0 || total_below == 0 || total_above == 0) return(0)
  prefix_below <- as.numeric(cumsum(below)[-n])
  prefix_above <- as.numeric(cumsum(above)[-n])
  n_a <- prefix_below + prefix_above
  n_b <- effective_n - n_a
  determinant <- prefix_below * (total_above - prefix_above) -
    prefix_above * (total_below - prefix_below)
  denominator <- n_a * n_b * total_below * total_above
  statistic <- effective_n * determinant^2 / denominator
  statistic[!is.finite(statistic)] <- 0
  max(statistic)
}

wilcoxon_statistic <- function(x) {
  n <- length(x)
  if (n < 2L) return(0)
  ranks <- rank(x, ties.method = "average")
  split <- seq_len(n - 1L)
  max(abs(cumsum(ranks)[split] - split * (n + 1) / 2))
}

simulate_null_statistics <- function(
    residuals, changepoints, block_length, reps, seed, workers) {
  starts <- valid_block_starts(length(residuals), block_length, changepoints)
  rows <- parallel_map(seq_len(reps), function(i) {
    sample <- sample_blocks(residuals, block_length, starts, seed + i)
    c(
      MedianMood = well_log_mood_right_critical_cpp(sample),
      Wilcoxon = well_log_wilcoxon_right_critical_cpp(sample)
    )
  }, workers)
  do.call(rbind, rows)
}

calibrate_null <- function(
    residuals, changepoints, block_lengths, calibration_reps,
    validation_reps, calibration_seed, validation_seed, workers) {
  tables <- list()
  statistics <- list()
  row_id <- 0L
  for (block_length in block_lengths) {
    calibration <- simulate_null_statistics(
      residuals, changepoints, block_length, calibration_reps,
      calibration_seed + 1000L * block_length, workers
    )
    validation <- simulate_null_statistics(
      residuals, changepoints, block_length, validation_reps,
      validation_seed + 1000L * block_length, workers
    )
    statistics[[as.character(block_length)]] <- list(
      calibration = calibration, validation = validation
    )
    for (test in colnames(calibration)) {
      raw_quantile <- unname(stats::quantile(
        calibration[, test], WELL_LOG_NULL_TARGET, type = 1
      ))
      gamma <- raw_quantile + max(1e-10, abs(raw_quantile) * 1e-12)
      calibration_success <- sum(calibration[, test] < gamma)
      validation_success <- sum(validation[, test] < gamma)
      interval <- wilson_interval(validation_success, validation_reps)
      row_id <- row_id + 1L
      tables[[row_id]] <- data.frame(
        test = test,
        subtests = "right",
        criterion = "maximum over all growing prefixes",
        block_length = block_length,
        gamma = gamma,
        raw_quantile = raw_quantile,
        target = WELL_LOG_NULL_TARGET,
        calibration_recovery = calibration_success / calibration_reps,
        validation_recovery = validation_success / validation_reps,
        validation_lower = interval[["lower"]],
        validation_upper = interval[["upper"]],
        calibration_reps = calibration_reps,
        validation_reps = validation_reps,
        calibration_seed = calibration_seed + 1000L * block_length,
        validation_seed = validation_seed + 1000L * block_length,
        stringsAsFactors = FALSE
      )
    }
  }
  list(summary = do.call(rbind, tables), statistics = statistics)
}

legacy_mood_threshold <- function(n) {
  segment_length <- n / 10 - 1
  alpha_split <- 1 - (1 - 0.01)^(1 / segment_length)
  stats::qchisq(1 - alpha_split, df = 10)
}

legacy_wilcoxon_threshold <- function(n) {
  1.5 * sqrt((n / 30)^3 / 12)
}

reverse_boundaries <- function(y, gamma, mode) {
  n <- length(y)
  sort(n - fit_svp(rev(y), gamma, "MedianMoodCost", mode))
}

observed_comparison <- function(
    y, reference, mood_gamma, wilcoxon_gamma, scale) {
  n <- length(y)
  fits <- list(
    "Gaussian PELT" = fit_pelt(y, scale),
    "RFPOP reference" = reference,
    "Legacy SVP Wilcoxon / both" = fit_svp(
      y, legacy_wilcoxon_threshold(n), "WilcoxonCost", "both"
    ),
    "Legacy SVP MedianMood / both" = fit_svp(
      y, legacy_mood_threshold(n), "MedianMoodCost", "both"
    ),
    "Calibrated SVP Wilcoxon / right" = fit_svp(
      y, wilcoxon_gamma, "WilcoxonCost", "right"
    ),
    "Calibrated SVP MedianMood / right" = fit_svp(
      y, mood_gamma, "MedianMoodCost", "right"
    ),
    "SVP MedianMood / both at right gamma" = fit_svp(
      y, mood_gamma, "MedianMoodCost", "both"
    )
  )
  rows <- lapply(names(fits), function(method) {
    cp <- fits[[method]]
    scores <- match_changepoints(reference, cp, WELL_LOG_MATCH_TOLERANCE)
    data.frame(
      method = method,
      changepoints = paste(cp, collapse = ","),
      n_changepoints = length(cp),
      matched = unname(scores[["tp"]]),
      precision = unname(scores[["precision"]]),
      recall = unname(scores[["recall"]]),
      F1 = unname(scores[["F1"]]),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  table$reverse_changepoints <- NA_character_
  table$reversal_exact <- NA
  table$reversal_F1 <- NA_real_
  method_names <- c(
    right = "Calibrated SVP MedianMood / right",
    both = "SVP MedianMood / both at right gamma"
  )
  for (mode in names(method_names)) {
    method <- unname(method_names[[mode]])
    forward <- fits[[method]]
    reversed <- reverse_boundaries(y, mood_gamma, mode)
    reverse_scores <- match_changepoints(
      forward, reversed, WELL_LOG_MATCH_TOLERANCE
    )
    index <- table$method == method
    table$reverse_changepoints[index] <- paste(reversed, collapse = ",")
    table$reversal_exact[index] <- identical(forward, reversed)
    table$reversal_F1[index] <- unname(reverse_scores[["F1"]])
  }
  list(table = table, fits = fits)
}

run_mode_bootstrap <- function(
    signal, residuals, changepoints, block_length, gamma, reps, seed, workers) {
  starts <- valid_block_starts(length(residuals), block_length, changepoints)
  results <- parallel_map(seq_len(reps), function(i) {
    simulated <- signal + sample_blocks(
      residuals, block_length, starts, seed + i
    )
    lapply(c("right", "both"), function(mode) {
      cp <- fit_svp(simulated, gamma, "MedianMoodCost", mode)
      scores <- match_changepoints(
        changepoints, cp, WELL_LOG_MATCH_TOLERANCE
      )
      list(
        row = data.frame(
          rep = i, mode = mode, n_changepoints = length(cp),
          precision = unname(scores[["precision"]]),
          recall = unname(scores[["recall"]]),
          F1 = unname(scores[["F1"]]),
          changepoints = paste(cp, collapse = ","),
          stringsAsFactors = FALSE
        ),
        changepoints = cp
      )
    })
  }, workers)
  flattened <- do.call(c, results)
  rows <- do.call(rbind, lapply(flattened, `[[`, "row"))
  list(
    rows = rows,
    changepoints = lapply(flattened, `[[`, "changepoints")
  )
}

run_both_null_stress <- function(
    residuals, changepoints, block_length, gamma, reps, seed, workers) {
  starts <- valid_block_starts(length(residuals), block_length, changepoints)
  values <- parallel_map(seq_len(reps), function(i) {
    sample <- sample_blocks(residuals, block_length, starts, seed + i)
    length(fit_svp(sample, gamma, "MedianMoodCost", "both")) == 0L
  }, workers)
  successes <- sum(unlist(values, use.names = FALSE))
  interval <- wilson_interval(successes, reps)
  data.frame(
    subtests = "both",
    gamma_source = "right-calibrated MedianMood threshold",
    gamma = gamma,
    one_segment_recovery = successes / reps,
    lower = interval[["lower"]],
    upper = interval[["upper"]],
    reps = reps,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

selection_frequencies <- function(bootstrap, reference, observed_both) {
  extras <- observed_both[!vapply(observed_both, function(cp) {
    any(abs(reference - cp) <= WELL_LOG_MATCH_TOLERANCE)
  }, logical(1))]
  targets <- sort(c(reference, extras))
  target_type <- ifelse(targets %in% reference, "RFPOP reference", "SVP only")
  rows <- list()
  k <- 0L
  for (mode in c("right", "both")) {
    indices <- which(bootstrap$rows$mode == mode)
    for (j in seq_along(targets)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        mode = mode,
        target = targets[j],
        target_type = target_type[j],
        selection_frequency = mean(vapply(indices, function(i) {
          any(abs(bootstrap$changepoints[[i]] - targets[j]) <=
                WELL_LOG_MATCH_TOLERANCE)
        }, logical(1))),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

threshold_sensitivity <- function(y, reference, gamma, workers) {
  grid <- seq(max(1, floor(gamma) - 20), ceiling(gamma) + 20, by = 1)
  results <- parallel_map(grid, function(value) {
    do.call(rbind, lapply(c("right", "both"), function(mode) {
      cp <- fit_svp(y, value, "MedianMoodCost", mode)
      scores <- match_changepoints(
        reference, cp, WELL_LOG_MATCH_TOLERANCE
      )
      data.frame(
        gamma = value, mode = mode, n_changepoints = length(cp),
        precision = unname(scores[["precision"]]),
        recall = unname(scores[["recall"]]),
        F1 = unname(scores[["F1"]]),
        stringsAsFactors = FALSE
      )
    }))
  }, workers)
  do.call(rbind, results)
}

piecewise_median <- function(y, changepoints) {
  segment_median_fit(y, changepoints)$fitted
}

write_review_plots <- function(
    y, observed, calibration, selected_gamma, bootstrap,
    frequencies, sensitivity) {
  dir.create(WELL_LOG_PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(WELL_LOG_PLOT_ROOT, "01_segmentations.pdf"),
      width = 10, height = 8)
  old_par <- par(mfrow = c(4, 1), mar = c(2.6, 4.2, 2.2, 0.8),
                 oma = c(2, 0, 0, 0))
  panels <- c(
    "Gaussian PELT", "RFPOP reference",
    "Calibrated SVP MedianMood / right",
    "SVP MedianMood / both at right gamma"
  )
  colors <- c("#4477AA", "#228833", "#CC6677", "#AA3377")
  for (i in seq_along(panels)) {
    cp <- observed$fits[[panels[i]]]
    plot(y, type = "p", pch = ".", cex = 1.8, xlab = "", ylab = "",
         main = panels[i])
    lines(piecewise_median(y, cp), col = colors[i], lwd = 2)
    abline(v = cp, col = "grey25", lty = 3)
  }
  mtext("Sequence position", side = 1, outer = TRUE)
  par(old_par)
  dev.off()

  curves <- list()
  k <- 0L
  probabilities <- seq(0.90, 0.999, length.out = 80)
  for (block in names(calibration$statistics)) {
    values <- calibration$statistics[[block]]$calibration
    for (test in colnames(values)) {
      k <- k + 1L
      curves[[k]] <- data.frame(
        gamma = as.numeric(stats::quantile(values[, test], probabilities,
                                           type = 1)),
        recovery = probabilities,
        test = test,
        block_length = factor(block)
      )
    }
  }
  curve_data <- do.call(rbind, curves)
  p <- ggplot2::ggplot(
    curve_data,
    ggplot2::aes(gamma, recovery, color = block_length)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = WELL_LOG_NULL_TARGET, linetype = 2) +
    ggplot2::facet_wrap(~test, scales = "free_x") +
    ggplot2::labs(x = "Validity threshold", y = "Null one-segment recovery",
                  color = "Block length") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(WELL_LOG_PLOT_ROOT, "02_null_calibration.pdf"),
                  p, width = 8, height = 4.5)

  stacked <- rbind(
    data.frame(mode = bootstrap$rows$mode, metric = "Precision",
               value = bootstrap$rows$precision),
    data.frame(mode = bootstrap$rows$mode, metric = "Recall",
               value = bootstrap$rows$recall),
    data.frame(mode = bootstrap$rows$mode, metric = "F1",
               value = bootstrap$rows$F1)
  )
  p <- ggplot2::ggplot(stacked, ggplot2::aes(mode, value, fill = mode)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.size = 0.5) +
    ggplot2::facet_wrap(~metric) +
    ggplot2::scale_fill_manual(values = c(right = "#CC6677", both = "#AA3377")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "SVP subtests mode", y = "Reference-bootstrap score") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(WELL_LOG_PLOT_ROOT, "03_mode_bootstrap.pdf"),
                  p, width = 8, height = 4.2)

  p <- ggplot2::ggplot(
    frequencies,
    ggplot2::aes(factor(target), selection_frequency, fill = mode)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(. ~ target_type, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = c(right = "#CC6677", both = "#AA3377")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Candidate changepoint", y = "Selection frequency",
                  fill = "SVP mode") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")
  ggplot2::ggsave(file.path(WELL_LOG_PLOT_ROOT, "04_selection_frequency.pdf"),
                  p, width = 9, height = 4.5)

  threshold_long <- rbind(
    data.frame(gamma = sensitivity$gamma, mode = sensitivity$mode,
               metric = "F1", value = sensitivity$F1),
    data.frame(gamma = sensitivity$gamma, mode = sensitivity$mode,
               metric = "Number of changepoints",
               value = sensitivity$n_changepoints)
  )
  p <- ggplot2::ggplot(
    threshold_long, ggplot2::aes(gamma, value, color = mode)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = selected_gamma, linetype = 2) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(right = "#CC6677", both = "#AA3377")) +
    ggplot2::labs(x = "Median-Mood threshold", y = NULL, color = "SVP mode") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(WELL_LOG_PLOT_ROOT, "05_threshold_sensitivity.pdf"),
                  p, width = 8, height = 4.2)
}

validate_scan_implementations <- function(
    residuals, changepoints, block_length, seed) {
  starts <- valid_block_starts(length(residuals), block_length, changepoints)
  sample <- sample_blocks(residuals, block_length, starts, seed)
  short <- sample[seq_len(80L)]
  mood_reference <- max(vapply(seq.int(2L, length(short)), function(end) {
    median_mood_statistic(short[seq_len(end)])
  }, numeric(1)))
  wilcoxon_reference <- max(vapply(seq.int(2L, length(short)), function(end) {
    wilcoxon_statistic(short[seq_len(end)])
  }, numeric(1)))
  stopifnot(
    isTRUE(all.equal(
      mood_reference, well_log_mood_right_critical_cpp(short),
      tolerance = 1e-12
    )),
    isTRUE(all.equal(
      wilcoxon_reference, well_log_wilcoxon_right_critical_cpp(short),
      tolerance = 1e-12
    ))
  )
  for (test in c("MedianMoodCost", "WilcoxonCost")) {
    statistic <- if (test == "MedianMoodCost") {
      well_log_mood_right_critical_cpp(sample)
    } else {
      well_log_wilcoxon_right_critical_cpp(sample)
    }
    high <- fit_svp(sample, statistic + max(1e-8, abs(statistic) * 1e-10),
                    test, "right")
    low <- fit_svp(sample, max(1e-12,
                              statistic - max(1e-8, abs(statistic) * 1e-10)),
                   test, "right")
    stopifnot(length(high) == 0L, length(low) > 0L)
  }
  invisible(TRUE)
}

run_well_log_review <- function(
    workers = 8L,
    calibration_reps = WELL_LOG_N_CALIBRATION,
    validation_reps = WELL_LOG_N_VALIDATION,
    power_reps = WELL_LOG_N_POWER,
    both_null_reps = WELL_LOG_N_BOTH_NULL_STRESS) {
  require_well_log_packages()
  dir.create(WELL_LOG_ROOT, recursive = TRUE, showWarnings = FALSE)
  dir.create(WELL_LOG_PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)

  y <- as.numeric(DeCAFS::oilWell)
  n <- length(y)
  primary_scale <- stats::mad(y, constant = 1.4826)
  difference_scale <- stats::mad(diff(y), constant = 1.4826) / sqrt(2)
  reference_fit <- fit_rfpop(y, primary_scale)
  reference <- reference_fit$changepoints
  signal <- segment_median_fit(y / primary_scale, reference)
  residuals <- y / primary_scale - signal$fitted
  selected_block <- select_block_length(residuals)
  block_lengths <- sort(unique(c(20L, selected_block, 64L)))

  calibration <- calibrate_null(
    residuals, reference, block_lengths,
    calibration_reps, validation_reps,
    WELL_LOG_CALIBRATION_SEED, WELL_LOG_VALIDATION_SEED, workers
  )
  selected_row <- calibration$summary[
    calibration$summary$test == "MedianMood" &
      calibration$summary$block_length == selected_block,
    , drop = FALSE
  ]
  stopifnot(nrow(selected_row) == 1L)
  selected_gamma <- selected_row$gamma
  selected_wilcoxon_row <- calibration$summary[
    calibration$summary$test == "Wilcoxon" &
      calibration$summary$block_length == selected_block,
    , drop = FALSE
  ]
  stopifnot(nrow(selected_wilcoxon_row) == 1L)
  selected_wilcoxon_gamma <- selected_wilcoxon_row$gamma

  validate_scan_implementations(
    residuals, reference, selected_block, WELL_LOG_VALIDATION_SEED
  )
  observed <- observed_comparison(
    y, reference, selected_gamma, selected_wilcoxon_gamma, primary_scale
  )
  bootstrap <- run_mode_bootstrap(
    signal$fitted, residuals, reference, selected_block, selected_gamma,
    power_reps, WELL_LOG_POWER_SEED, workers
  )
  frequencies <- selection_frequencies(
    bootstrap, reference,
    observed$fits[["SVP MedianMood / both at right gamma"]]
  )
  both_null_stress <- run_both_null_stress(
    residuals, reference, selected_block, selected_gamma,
    both_null_reps, WELL_LOG_BOTH_NULL_SEED, workers
  )
  sensitivity <- threshold_sensitivity(
    y, reference, selected_gamma, workers
  )

  rfpop_sensitivity <- do.call(rbind, lapply(list(
    package_global_mad = primary_scale,
    first_difference_mad = difference_scale
  ), function(scale) {
    fit <- fit_rfpop(y, scale)
    data.frame(
      scale = scale,
      n_changepoints = length(fit$changepoints),
      changepoints = paste(fit$changepoints, collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  rfpop_sensitivity$scale_method <- rownames(rfpop_sensitivity)
  rownames(rfpop_sensitivity) <- NULL

  bootstrap_summary <- do.call(rbind, lapply(c("right", "both"), function(mode) {
    data <- bootstrap$rows[bootstrap$rows$mode == mode, ]
    data.frame(
      mode = mode,
      mean_precision = mean(data$precision),
      mean_recall = mean(data$recall),
      mean_F1 = mean(data$F1),
      median_F1 = stats::median(data$F1),
      F1_lower_10 = unname(stats::quantile(data$F1, 0.1)),
      F1_upper_90 = unname(stats::quantile(data$F1, 0.9)),
      mean_changepoints = mean(data$n_changepoints),
      reps = nrow(data),
      stringsAsFactors = FALSE
    )
  }))

  write_review_plots(
    y, observed, calibration, selected_gamma, bootstrap,
    frequencies, sensitivity
  )

  metadata <- list(
    generated = as.character(Sys.time()),
    n = n,
    data_sha256 = if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(y, algo = "sha256", serialize = TRUE)
    } else NA_character_,
    svpChange_version = as.character(utils::packageVersion("svpChange")),
    svpChange_library = find.package("svpChange"),
    robseg_version = as.character(utils::packageVersion("robseg")),
    primary_scale = primary_scale,
    difference_scale = difference_scale,
    selected_block_length = selected_block,
    selected_gamma = selected_gamma,
    selected_wilcoxon_gamma = selected_wilcoxon_gamma,
    null_target = WELL_LOG_NULL_TARGET,
    match_tolerance = WELL_LOG_MATCH_TOLERANCE,
    calibration_reps = calibration_reps,
    validation_reps = validation_reps,
    power_reps = power_reps,
    both_null_reps = both_null_reps,
    calibration_seed = WELL_LOG_CALIBRATION_SEED,
    validation_seed = WELL_LOG_VALIDATION_SEED,
    power_seed = WELL_LOG_POWER_SEED,
    both_null_seed = WELL_LOG_BOTH_NULL_SEED
  )
  result <- list(
    metadata = metadata,
    reference = reference,
    calibration = calibration,
    observed = observed,
    bootstrap = bootstrap,
    bootstrap_summary = bootstrap_summary,
    both_null_stress = both_null_stress,
    frequencies = frequencies,
    threshold_sensitivity = sensitivity,
    rfpop_scale_sensitivity = rfpop_sensitivity
  )
  saveRDS(result, file.path(WELL_LOG_ROOT, "well_log_review.rds"))
  utils::write.csv(calibration$summary,
                   file.path(WELL_LOG_ROOT, "calibration_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(observed$table,
                   file.path(WELL_LOG_ROOT, "observed_comparison.csv"),
                   row.names = FALSE)
  utils::write.csv(bootstrap$rows,
                   file.path(WELL_LOG_ROOT, "bootstrap_mode_comparison.csv"),
                   row.names = FALSE)
  utils::write.csv(bootstrap_summary,
                   file.path(WELL_LOG_ROOT, "bootstrap_mode_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(both_null_stress,
                   file.path(WELL_LOG_ROOT, "both_null_stress.csv"),
                   row.names = FALSE)
  utils::write.csv(frequencies,
                   file.path(WELL_LOG_ROOT, "changepoint_frequencies.csv"),
                   row.names = FALSE)
  utils::write.csv(sensitivity,
                   file.path(WELL_LOG_ROOT, "threshold_sensitivity.csv"),
                   row.names = FALSE)
  utils::write.csv(rfpop_sensitivity,
                   file.path(WELL_LOG_ROOT, "rfpop_scale_sensitivity.csv"),
                   row.names = FALSE)

  stopifnot(
    n == 4050L,
    selected_block == 35L,
    nrow(bootstrap$rows) == 2L * power_reps
  )
  if (calibration_reps == WELL_LOG_N_CALIBRATION &&
      validation_reps == WELL_LOG_N_VALIDATION &&
      power_reps == WELL_LOG_N_POWER) {
    stopifnot(
      abs(selected_row$validation_recovery - WELL_LOG_NULL_TARGET) <= 0.01
    )
  }
  invisible(result)
}

if (identical(tolower(Sys.getenv("SVP_RUN_WELL_LOG_REVIEW")), "true")) {
  run_well_log_review()
}
