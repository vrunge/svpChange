## Shared infrastructure for the three paper power studies.

POWER_PATTERNS <- c("none", "up", "updown", "rand1")

power_default_workers <- function() {
  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected)) return(1L)
  max(1L, as.integer(detected) - 1L)
}

power_lapply <- function(X, FUN, workers = power_default_workers(), ...) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be one positive integer")
  }
  if (workers == 1L || .Platform$OS.type == "windows") {
    return(lapply(X, FUN, ...))
  }
  parallel::mclapply(
    X, FUN, ..., mc.cores = workers, mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

set_power_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
}

generate_power_signal <- function(n, pattern, jump, segments = 8L) {
  pattern <- match.arg(pattern, POWER_PATTERNS)
  if (pattern == "none") return(rep(0, n))
  if (pattern == "rand1") {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(42)
    raw_sizes <- stats::rpois(segments, 10)
    sizes <- pmax(round(raw_sizes * n / sum(raw_sizes)), 1)
    sizes[segments] <- sizes[segments] + n - sum(sizes)
    set.seed(43)
    levels <- stats::runif(segments, 0.5, 1) *
      sample(c(-1, 1), segments, replace = TRUE)
    return(rep(levels * jump, sizes))
  }
  if (n %% segments != 0) stop("n must be divisible by segments")
  levels <- if (pattern == "up") {
    seq(0, segments - 1) * jump
  } else {
    (seq(0, segments - 1) %% 2) * jump
  }
  rep(levels, each = n / segments)
}

localization_error <- function(cp_true, cp_est) {
  cp_true <- sort(as.integer(cp_true))
  cp_est <- sort(as.integer(cp_est))
  if (!length(cp_true) || length(cp_true) != length(cp_est)) return(NA_real_)
  max(abs(cp_est - cp_true))
}

normalise_boundaries <- function(boundaries, n) {
  sort(unique(c(as.integer(unlist(boundaries)), n)))
}

run_power_grid <- function(n, jump_sizes, reps, simulate_noise, fit_methods,
                           tolerance, workers, seed) {
  design <- expand.grid(
    pattern = POWER_PATTERNS,
    jump = jump_sizes,
    rep = seq_len(reps),
    stringsAsFactors = FALSE
  )
  rows <- power_lapply(seq_len(nrow(design)), function(i) {
    task <- design[i, ]
    set_power_seed(seed + i)
    mu <- generate_power_signal(n, task$pattern, task$jump)
    y <- mu + simulate_noise(n)
    truth <- normalise_boundaries(which(diff(mu) != 0), n)
    fits <- fit_methods(y, length(truth))
    if (is.null(fits$boundaries)) {
      fits <- list(
        boundaries = fits,
        fitted = lapply(fits, function(x) fitted_piecewise_mean(y, x))
      )
    }
    output <- dplyr::bind_rows(lapply(names(fits$boundaries), function(method) {
      result_row(
        task$pattern, task$jump, task$rep, method,
        fits$boundaries[[method]], truth, mu, fits$fitted[[method]], tolerance
      )
    }))
    output$changepoints <- unname(fits$boundaries)
    output
  }, workers = workers)
  dplyr::bind_rows(rows)
}

power_scenario_plot <- function(n, jump, simulate_noise, seed = 999L,
                                y_limits = NULL) {
  set_power_seed(seed)
  data <- dplyr::bind_rows(lapply(POWER_PATTERNS, function(pattern) {
    mu <- generate_power_signal(n, pattern, jump)
    data.frame(
      time = seq_len(n),
      value = mu + simulate_noise(n),
      mean = mu,
      pattern = pattern
    )
  }))
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  plot <- ggplot2::ggplot(data, ggplot2::aes(time, value)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.35) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean), colour = "red", linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~pattern, nrow = 2L) +
    ggplot2::labs(x = "Time (Index)", y = "Value") +
    paper_theme() +
    ggplot2::theme(legend.position = "none")
  if (!is.null(y_limits)) {
    plot <- plot + ggplot2::coord_cartesian(ylim = y_limits)
  }
  plot
}

internal_changepoints <- function(boundaries, terminal) {
  sort(unique(as.integer(boundaries[as.integer(boundaries) != terminal])))
}

maximum_tolerance_matches <- function(truth, estimate, tolerance) {
  truth <- sort(as.integer(truth))
  estimate <- sort(as.integer(estimate))
  i <- 1L
  j <- 1L
  matched <- 0L
  while (i <= length(truth) && j <= length(estimate)) {
    if (abs(truth[i] - estimate[j]) <= tolerance) {
      matched <- matched + 1L
      i <- i + 1L
      j <- j + 1L
    } else if (estimate[j] < truth[i] - tolerance) {
      j <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  matched
}

paper_metrics <- function(true_boundaries, estimated_boundaries, tolerance,
                          terminal = max(c(true_boundaries,
                                           estimated_boundaries))) {
  truth <- internal_changepoints(true_boundaries, terminal)
  estimate <- internal_changepoints(estimated_boundaries, terminal)
  matched <- maximum_tolerance_matches(truth, estimate, tolerance)

  if (!length(truth) && !length(estimate)) {
    precision <- recall <- f1 <- 1
  } else if (!length(truth) || !length(estimate)) {
    precision <- recall <- f1 <- 0
  } else {
    precision <- matched / length(estimate)
    recall <- matched / length(truth)
    f1 <- if (precision + recall > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      0
    }
  }
  c(Precision = precision, Recall = recall, F1 = f1,
    CorrectNumCP = as.numeric(length(estimate) == length(truth)),
    LocalizationError = localization_error(truth, estimate))
}

fitted_piecewise_mean <- function(y, boundaries) {
  starts <- c(1L, head(boundaries, -1L) + 1L)
  unlist(Map(function(first, last) {
    rep(mean(y[first:last]), last - first + 1L)
  }, starts, boundaries), use.names = FALSE)
}

result_row <- function(pattern, jump, rep, algorithm, boundaries,
                       truth, mu, fitted, tolerance) {
  metric <- paper_metrics(truth, boundaries, tolerance, terminal = length(mu))
  data.frame(
    pattern = pattern, jump = jump, rep = rep, algorithm = algorithm,
    Precision = metric[["Precision"]], Recall = metric[["Recall"]],
    F1 = metric[["F1"]], CorrectNumCP = metric[["CorrectNumCP"]],
    LocalizationError = metric[["LocalizationError"]],
    MSE = mean((mu - fitted)^2), NumSegments = length(boundaries),
    stringsAsFactors = FALSE
  )
}

paper_theme <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

empty_power_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message) +
    ggplot2::theme_void()
}

summarise_metric <- function(results, metric) {
  summary <- results |>
    dplyr::group_by(algorithm, pattern, jump) |>
    dplyr::summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),
      se = stats::sd(.data[[metric]], na.rm = TRUE) /
        sqrt(sum(is.finite(.data[[metric]]))),
      lower = mean - 1.96 * se, upper = mean + 1.96 * se,
      .groups = "drop"
    )
  if (metric %in% c("Precision", "Recall", "F1", "CorrectNumCP")) {
    summary$lower <- pmax(0, summary$lower)
    summary$upper <- pmin(1, summary$upper)
  } else if (metric == "MSE") {
    summary$lower <- pmax(0, summary$lower)
  }
  summary
}

metric_plot <- function(results, metric, y_label = metric,
                        probability_scale = FALSE) {
  data <- summarise_metric(results, metric)
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(jump, mean, colour = algorithm, fill = algorithm)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::facet_wrap(~pattern, nrow = 2L) +
    ggplot2::labs(x = "Jump size", y = y_label,
                  colour = "algorithm", fill = "algorithm") +
    paper_theme()
  if (probability_scale) plot <- plot + ggplot2::coord_cartesian(ylim = c(0, 1))
  plot
}

precision_recall_plot <- function(results) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required for the combined precision/recall figure")
  }
  precision <- metric_plot(results, "Precision", "Precision", TRUE) +
    ggplot2::theme(legend.position = "none")
  recall <- metric_plot(results, "Recall", "Recall", TRUE)
  patchwork::wrap_plots(precision, recall, ncol = 1L,
                        guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

localization_error_plot <- function(results) {
  results <- complete_localization_error(results)
  eligible <- results[
    results$CorrectNumCP == 1 & results$pattern != "none" &
      is.finite(results$LocalizationError),
  ]
  if (!nrow(eligible)) {
    return(empty_power_plot(
      "No fits with the correct number of changepoints"
    ))
  }
  eligible$pattern <- factor(
    eligible$pattern, c("up", "updown", "rand1")
  )
  summary <- eligible |>
    dplyr::group_by(algorithm, pattern, jump) |>
    dplyr::summarise(
      mean = mean(LocalizationError),
      se = if (dplyr::n() > 1L) {
        stats::sd(LocalizationError) / sqrt(dplyr::n())
      } else {
        0
      },
      lower = pmax(0, mean - 1.96 * se),
      upper = mean + 1.96 * se,
      .groups = "drop"
    )
  ggplot2::ggplot(
    summary, ggplot2::aes(jump, mean, colour = algorithm, fill = algorithm)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::facet_wrap(~pattern, nrow = 1L) +
    ggplot2::labs(
      x = "Jump size",
      y = "Maximum localization error | correct number of changes",
      colour = "algorithm", fill = "algorithm"
    ) +
    paper_theme()
}

add_localization_error <- function(results, n) {
  results$LocalizationError <- vapply(seq_len(nrow(results)), function(i) {
    mu <- generate_power_signal(n, results$pattern[i], results$jump[i])
    truth <- which(diff(mu) != 0)
    estimate <- as.integer(results$changepoints[[i]])
    estimate <- estimate[estimate < n]
    localization_error(truth, estimate)
  }, numeric(1))
  results
}

complete_localization_error <- function(results) {
  if (!is.list(results$changepoints)) {
    stop("A changepoints list column is required for localization errors")
  }
  terminal_boundaries <- vapply(
    results$changepoints,
    function(x) if (length(x)) max(as.integer(x)) else NA_integer_,
    integer(1)
  )
  n <- max(terminal_boundaries, na.rm = TRUE)
  add_localization_error(results, n)
}

set_algorithm_order <- function(results, order, drop = character()) {
  results <- results[!(as.character(results$algorithm) %in% drop), ,
                     drop = FALSE]
  algorithms <- as.character(results$algorithm)
  unexpected <- setdiff(unique(algorithms), order)
  if (length(unexpected)) {
    stop("Missing display order for: ", paste(unexpected, collapse = ", "))
  }
  results$algorithm <- factor(algorithms, levels = order)
  results
}

changepoint_distribution_plot <- function(results, n, selected_jump = 0.6,
                                          bins = 100L) {
  selected <- results[
    abs(results$jump - selected_jump) < 1e-10 & results$pattern != "none",
  ]
  rows <- lapply(seq_len(nrow(selected)), function(i) {
    cps <- selected$changepoints[[i]]
    cps <- cps[cps < n]
    if (!length(cps)) return(NULL)
    data.frame(pattern = selected$pattern[i],
               algorithm = selected$algorithm[i], changepoint = cps)
  })
  data <- dplyr::bind_rows(rows)
  if (!nrow(data)) {
    return(empty_power_plot(
      "No changepoints detected for the selected jump"
    ))
  }
  data$pattern <- factor(data$pattern, c("rand1", "up", "updown"))
  if (is.factor(results$algorithm)) {
    data$algorithm <- factor(as.character(data$algorithm),
                             levels = levels(results$algorithm))
  }
  ggplot2::ggplot(data, ggplot2::aes(changepoint)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue",
                            alpha = 0.75, colour = "black", linewidth = 0.15) +
    ggplot2::facet_grid(pattern ~ algorithm) +
    ggplot2::labs(x = "Sequence Position (Time)",
                  y = "Frequency of Detected Change") +
    paper_theme() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

save_power_outputs <- function(results, root, scenario_plot,
                               selected_jump = 0.6,
                               plot_results = results,
                               save_results = TRUE) {
  dir.create(file.path(root, "plots"), recursive = TRUE, showWarnings = FALSE)
  results <- complete_localization_error(results)
  plot_results <- complete_localization_error(plot_results)
  if (save_results) {
    saveRDS(results, file.path(root, "results.rds"))
    utils::write.csv(dplyr::select(results, -changepoints),
                     file.path(root, "results.csv"), row.names = FALSE)
  }
  mse_plot <- metric_plot(plot_results, "MSE", "Signal MSE")
  mse_summary <- summarise_metric(plot_results, "MSE")
  positive <- mse_summary$mean[is.finite(mse_summary$mean) & mse_summary$mean > 0]
  if (length(positive) && max(positive) > 20 * stats::median(positive)) {
    mse_plot <- mse_plot +
      ggplot2::scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 0.01, base = 10),
        breaks = c(0, 0.1, 1, 10, 100, 1000, 10000),
        labels = c("0", "0.1", "1", "10", "100", "1,000", "10,000")
      ) +
      ggplot2::labs(y = "Signal MSE (pseudo-log scale)")
  }
  plots <- list(
    "01_scenarios.pdf" = scenario_plot,
    "02_f1.pdf" = metric_plot(plot_results, "F1", "F1", TRUE),
    "03_precision_recall.pdf" = precision_recall_plot(plot_results),
    "04_changepoint_distributions.pdf" =
      changepoint_distribution_plot(
        plot_results, max(unlist(plot_results$changepoints)), selected_jump
      ),
    "05_correct_number_probability.pdf" =
      metric_plot(plot_results, "CorrectNumCP",
                  "Probability of the correct number of changepoints", TRUE),
    "06_signal_mse.pdf" = mse_plot,
    "07_loc_error.pdf" = localization_error_plot(plot_results)
  )
  sizes <- list(c(10, 5), c(10, 5), c(10, 10), c(15, 10), c(10, 5),
                c(10, 5), c(12, 4.5))
  for (i in seq_along(plots)) {
    ggplot2::ggsave(file.path(root, "plots", names(plots)[i]), plots[[i]],
                    width = sizes[[i]][1], height = sizes[[i]][2])
  }
  invisible(plots)
}
