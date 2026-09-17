## Runtime scaling under independent Gaussian noise.

library(svpChange)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "time_common.R"))
source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_TIME_ROOT <- file.path("simulations", "time_gaussian")

gaussian_time_methods <- function(constants = read_gaussian_calibration()) {
  c(
    "PELT",
    gaussian_method_label("none", constants[["none"]]),
    gaussian_method_label("right", constants[["right"]]),
    gaussian_method_label("both", constants[["both"]])
  )
}

simulate_gaussian_time_data <- function(n, changes = 0L, jump = 10) {
  sizes <- balanced_segment_sizes(n, changes)
  ts_generator(
    chpts = cumsum(sizes),
    parameters = rep(c(0, jump), length.out = changes + 1L),
    sd_noise = 1,
    type = "gauss"
  )
}

fit_gaussian_time_method <- function(
  y, method, constants = read_gaussian_calibration()
) {
  n <- length(y)
  if (method == "PELT") {
    fit <- changepoint::cpt.mean(
      y,
      method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
    )
    return(sort(unique(c(changepoint::cpts(fit), n))))
  }
  mode <- names(constants)[vapply(names(constants), function(candidate) {
    identical(method, gaussian_method_label(candidate, constants[[candidate]]))
  }, logical(1))]
  if (length(mode) != 1L) stop("Unknown Gaussian timing method: ", method)
  fit <- SVP(
    y, constants[[mode]] * log(n), "gaussian_mean",
    subtests = mode
  )
  extract_svp_boundaries(fit, n)
}

run_time_gaussian <- function(
  n_values = round(2^(seq(8, 14, length.out = 24))),
  k_values = 0:100,
  fixed_n = 10000L,
  reps_n = 12L,
  reps_k = 12L,
  reps_none_n = 48L,
  reps_none_k = 48L,
  seed = 456L, min_batch_time = 0.05,
  constants = read_gaussian_calibration()
) {
  if (reps_none_n < reps_n || reps_none_k < reps_k) {
    stop("SVP none must have at least as many replicates as the other methods")
  }
  methods <- gaussian_time_methods(constants)
  none_method <- gaussian_method_label("none", constants[["none"]])
  set.seed(seed - 1L)
  warm_data <- simulate_gaussian_time_data(256L, 0L)
  invisible(lapply(methods, function(method) {
    fit_gaussian_time_method(warm_data, method, constants)
  }))
  common <- benchmark_time_study(
    n_values, k_values, fixed_n, reps_n, reps_k,
    methods = methods,
    simulate_data = simulate_gaussian_time_data,
    fit_method = function(y, method, changes) {
      fit_gaussian_time_method(y, method, constants)
    },
    seed = seed, min_batch_time = min_batch_time,
    randomize_method_order = TRUE
  )
  extra_none <- benchmark_time_study(
    n_values, k_values, fixed_n,
    reps_none_n - reps_n, reps_none_k - reps_k,
    methods = none_method,
    simulate_data = simulate_gaussian_time_data,
    fit_method = function(y, method, changes) {
      fit_gaussian_time_method(y, method, constants)
    },
    seed = seed + 1000000L, min_batch_time = min_batch_time,
    randomize_method_order = FALSE
  )
  if (nrow(extra_none)) {
    extra_none$rep <- extra_none$rep + ifelse(
      extra_none$experiment == "vary_n", reps_n, reps_k
    )
  }
  dplyr::bind_rows(common, extra_none)
}

bootstrap_median_interval <- function(x, draws = 1000L) {
  values <- replicate(draws, stats::median(sample(x, replace = TRUE)))
  stats::quantile(values, c(0.025, 0.975), names = FALSE)
}

summarise_gaussian_time_n <- function(results, draws = 1000L, seed = 901L) {
  set.seed(seed)
  results |>
    dplyr::filter(experiment == "vary_n", is.finite(time), time > 0) |>
    dplyr::group_by(method, n) |>
    dplyr::group_modify(function(data, key) {
      interval <- bootstrap_median_interval(data$time, draws)
      data.frame(
        replicates = nrow(data),
        mean_time = mean(data$time),
        median_time = stats::median(data$time),
        lower = interval[1L], upper = interval[2L]
      )
    }) |>
    dplyr::ungroup()
}

summarise_gaussian_time_detected <- function(
    results, draws = 1000L, seed = 903L, min_replicates = 4L) {
  set.seed(seed)
  results |>
    dplyr::filter(
      experiment == "vary_k", is.finite(time), time > 0
    ) |>
    dplyr::group_by(method, detected) |>
    dplyr::filter(dplyr::n() >= min_replicates) |>
    dplyr::group_modify(function(data, key) {
      interval <- bootstrap_median_interval(data$time, draws)
      data.frame(
        replicates = nrow(data),
        mean_time = mean(data$time),
        median_time = stats::median(data$time),
        lower = interval[1L], upper = interval[2L]
      )
    }) |>
    dplyr::ungroup()
}

summarise_gaussian_time_true_changes <- function(
    results, draws = 1000L, seed = 904L) {
  set.seed(seed)
  results |>
    dplyr::filter(
      experiment == "vary_k", is.finite(time), time > 0
    ) |>
    dplyr::group_by(method, k) |>
    dplyr::group_modify(function(data, key) {
      interval <- bootstrap_median_interval(data$time, draws)
      data.frame(
        replicates = nrow(data),
        mean_time = mean(data$time),
        median_time = stats::median(data$time),
        lower = interval[1L], upper = interval[2L]
      )
    }) |>
    dplyr::ungroup()
}

plot_gaussian_time_true_changes <- function(
    results, root = GAUSSIAN_TIME_ROOT,
    paper_root = normalizePath(file.path("..", "SVP_NEW_Figures2"),
                               mustWork = FALSE),
    prefix = "gaussian_time") {
  summary_k <- summarise_gaussian_time_true_changes(results)
  method_colours <- time_method_palette(summary_k$method)
  panel <- ggplot2::ggplot(
    summary_k, ggplot2::aes(k, median_time, colour = method, fill = method)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "True number of changes",
      y = "Median elapsed time (s)",
      colour = "Method", fill = "Method"
    ) +
    ggplot2::scale_colour_manual(values = method_colours) +
    ggplot2::scale_fill_manual(values = method_colours) +
    ggplot2::theme_minimal(base_size = 14L) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )
  dir.create(file.path(root, "plots"), recursive = TRUE,
             showWarnings = FALSE)
  dir.create(paper_root, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(root, "plots", paste0(prefix, "_time_vs_true_changes.pdf")),
    panel, width = 9, height = 5.2
  )
  ggplot2::ggsave(
    file.path(paper_root, "2_gaussian_time_vs_true_changes.pdf"),
    panel, width = 9, height = 5.2
  )
  utils::write.csv(
    summary_k,
    file.path(root, "gaussian_time_summary_by_true_changes.csv"),
    row.names = FALSE
  )
  invisible(list(panel = panel, summary = summary_k))
}

bootstrap_log_slope <- function(data, tail_only = FALSE, min_n = NULL,
                                draws = 1000L) {
  n_values <- sort(unique(data$n))
  if (!is.null(min_n)) {
    n_values <- n_values[n_values >= min_n]
  } else if (tail_only) {
    n_values <- n_values[n_values >= stats::median(n_values)]
  }
  if (length(n_values) < 2L) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  estimate <- function(resample = FALSE) {
    medians <- vapply(n_values, function(n_value) {
      values <- data$time[data$n == n_value]
      if (resample) values <- sample(values, replace = TRUE)
      stats::median(values)
    }, numeric(1))
    unname(stats::coef(stats::lm(log(medians) ~ log(n_values)))[2L])
  }
  samples <- replicate(draws, estimate(TRUE))
  interval <- stats::quantile(
    samples, c(0.025, 0.975),
    names = FALSE, na.rm = TRUE
  )
  c(estimate = estimate(FALSE), lower = interval[1L], upper = interval[2L])
}

gaussian_time_slopes <- function(results, draws = 1000L, seed = 902L,
                                 min_n = NULL) {
  set.seed(seed)
  data <- dplyr::filter(
    results, experiment == "vary_n", is.finite(time), time > 0
  )
  dplyr::bind_rows(lapply(unique(data$method), function(method) {
    method_data <- data[data$method == method, , drop = FALSE]
    ranges <- if (is.null(min_n)) c(FALSE, TRUE) else FALSE
    dplyr::bind_rows(lapply(ranges, function(tail_only) {
      values <- bootstrap_log_slope(
        method_data, tail_only, min_n = min_n, draws = draws
      )
      data.frame(
        method = method,
        range = if (!is.null(min_n)) {
          paste0("n >= ", min_n)
        } else if (tail_only) {
          "upper half"
        } else {
          "full grid"
        },
        exponent = values[["estimate"]],
        lower = values[["lower"]], upper = values[["upper"]]
      )
    }))
  }))
}

plot_gaussian_paper_time <- function(results, root, prefix = "gaussian_time") {
  paper_base_size <- 14
  summary_n <- summarise_gaussian_time_n(results)
  summary_detected <- summarise_gaussian_time_detected(results)
  method_colours <- time_method_palette(summary_n$method)
  slopes <- gaussian_time_slopes(results)
  tail_min_n <- 1000L
  tail_slopes <- gaussian_time_slopes(results, min_n = tail_min_n)
  slope_labels <- summary_n[summary_n$n == max(summary_n$n), , drop = FALSE]
  slope_labels <- merge(
    slope_labels,
    tail_slopes[, c("method", "exponent")],
    by = "method",
    all.x = TRUE,
    sort = FALSE
  )
  slope_labels$label <- paste0(
    sub(" / c =.*$", "", sub(
      "SVP Gaussian / ", "", slope_labels$method, fixed = TRUE
    )),
    ": ", sprintf("%.3f", slope_labels$exponent)
  )
  slope_labels$vjust <- ifelse(
    slope_labels$method == "PELT", -0.8,
    ifelse(grepl(" / none / ", slope_labels$method), 1.8,
      ifelse(grepl(" / right / ", slope_labels$method), -1.0, 0.5)
    )
  )
  panel_a <- ggplot2::ggplot(
    summary_n, ggplot2::aes(n, median_time, colour = method, fill = method)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::geom_label(
      data = slope_labels,
      ggplot2::aes(label = label),
      hjust = 1.05,
      vjust = slope_labels$vjust,
      fill = "white", alpha = 0.85, linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Sequence length", y = "Median elapsed time (s)", colour = "Method",
      fill = "Method", title = NULL, subtitle = NULL, caption = NULL
    ) +
    ggplot2::scale_colour_manual(values = method_colours) +
    ggplot2::scale_fill_manual(values = method_colours) +
    ggplot2::theme_minimal(base_size = paper_base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )

  panel_b <- ggplot2::ggplot(
    summary_detected,
    ggplot2::aes(
      detected, median_time, colour = method, fill = method
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Detected number of segments", y = "Median elapsed time (s)",
      colour = "Method", fill = "Method"
    ) +
    ggplot2::scale_colour_manual(values = method_colours) +
    ggplot2::scale_fill_manual(values = method_colours) +
    ggplot2::theme_minimal(base_size = paper_base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )

  combined <- combine_time_panels(
    panel_a, panel_b, base_size = paper_base_size
  )

  plot_root <- file.path(GAUSSIAN_TIME_ROOT, "plots")
  dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(plot_root, paste0(prefix, "_time_vs_n.pdf")),
    panel_a,
    width = 9, height = 5.2
  )
  ggplot2::ggsave(file.path(plot_root, paste0(prefix, "_time_vs_detected.pdf")),
    panel_b,
    width = 9, height = 5.2
  )
  ggplot2::ggsave(file.path(root, "2_exp1_time_vs_n.pdf"),
    panel_a,
    width = 9, height = 5.2
  )
  ggplot2::ggsave(file.path(root, "2_time_vs_detected_changes.pdf"),
    panel_b,
    width = 9, height = 5.2
  )
  ggplot2::ggsave(file.path(root, "2_gaussian_time_vs_n.pdf"),
    panel_a,
    width = 9, height = 5.2
  )
  ggplot2::ggsave(file.path(root, "2_gaussian_time_vs_detected.pdf"),
    panel_b,
    width = 9, height = 5.2
  )
  ggplot2::ggsave(file.path(plot_root, paste0(prefix, "_time_pair.pdf")),
    combined,
    width = 13.2, height = 5.6
  )
  ggplot2::ggsave(file.path(root, "2_gaussian_time.pdf"),
    combined,
    width = 13.2, height = 5.6
  )
  utils::write.csv(summary_n, file.path(
    GAUSSIAN_TIME_ROOT,
    "gaussian_time_summary_by_n.csv"
  ),
  row.names = FALSE
  )
  utils::write.csv(slopes, file.path(
    GAUSSIAN_TIME_ROOT,
    "gaussian_time_slopes.csv"
  ),
  row.names = FALSE
  )
  utils::write.csv(summary_detected, file.path(
    GAUSSIAN_TIME_ROOT,
    "gaussian_time_summary_by_detected.csv"
  ),
  row.names = FALSE
  )
  invisible(list(
    panel_a = panel_a, panel_b = panel_b, combined = combined,
    slopes = slopes
  ))
}

run_and_save_time_gaussian <- function(prefix = "gaussian_time") {
  started <- proc.time()[["elapsed"]]
  results <- run_time_gaussian()
  run_elapsed <- proc.time()[["elapsed"]] - started
  dir.create(GAUSSIAN_TIME_ROOT, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(GAUSSIAN_TIME_ROOT, paste0(prefix, "_results.rds")))
  utils::write.csv(
    results, file.path(GAUSSIAN_TIME_ROOT, paste0(prefix, "_results.csv")),
    row.names = FALSE
  )
  plot_gaussian_paper_time(
    results,
    normalizePath(file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE),
    prefix
  )
  manifest <- c(
    paste0("timestamp = ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("svpChange = ", as.character(utils::packageVersion("svpChange"))),
    paste0("changepoint = ", as.character(utils::packageVersion("changepoint"))),
    paste0("R = ", R.version.string),
    paste0("platform = ", R.version$platform),
    paste0("machine = ", Sys.info()[["machine"]]),
    paste0("system = ", Sys.info()[["sysname"]], " ", Sys.info()[["release"]]),
    paste0("git_sha = ", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    "n_grid = 256 to 16384 (24 values); fixed_n = 10000; k_grid = 0:100",
    "replicates = PELT/right/both: 12; none: 48",
    "timer = proc.time elapsed; adaptive batching target 0.05 s",
    paste0("simulation_elapsed_seconds = ", sprintf("%.3f", run_elapsed))
  )
  writeLines(
    manifest,
    file.path(GAUSSIAN_TIME_ROOT, paste0(prefix, "_manifest.txt"))
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_time_gaussian()
}
