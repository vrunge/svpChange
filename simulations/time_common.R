## Shared infrastructure for the three runtime studies.

balanced_segment_sizes <- function(n, changes) {
  segments <- changes + 1L
  sizes <- rep(n %/% segments, segments)
  remainder <- n - sum(sizes)
  if (remainder > 0L) {
    sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
  }
  sizes
}

alternating_mean <- function(n, changes = 0L, jump = 4) {
  sizes <- balanced_segment_sizes(n, changes)
  rep(rep(c(0, jump), length.out = changes + 1L), sizes)
}

extract_svp_boundaries <- function(fit, n) {
  if (is.null(fit) || is.null(fit$changepoints)) return(n)
  sort(unique(c(as.integer(unlist(fit$changepoints)), n)))
}

benchmark_time_study <- function(
    n_values, k_values, fixed_n, reps_n, reps_k, methods,
    simulate_data, fit_method, seed, min_batch_time = 0,
    randomize_method_order = FALSE) {
  total <- length(methods) *
    (length(n_values) * reps_n + length(k_values) * reps_k)
  if (total == 0L) return(data.frame())

  progress <- utils::txtProgressBar(min = 0, max = total, style = 3)
  on.exit(close(progress), add = TRUE)
  rows <- vector("list", total)
  completed <- 0L

  add_tasks <- function(experiment, values, reps, n_for, k_for, seed_for) {
    for (value in values) {
      for (rep in seq_len(reps)) {
        n <- n_for(value)
        k <- k_for(value)
        set.seed(seed_for(value, rep))
        data <- simulate_data(n, k)
        method_order <- if (randomize_method_order) sample(methods) else methods
        for (method in method_order) {
          gc(FALSE)
          batch <- 1L
          repeat {
            started <- proc.time()[["elapsed"]]
            for (batch_index in seq_len(batch)) {
              boundaries <- fit_method(data, method, k)
            }
            block_elapsed <- proc.time()[["elapsed"]] - started
            if (min_batch_time <= 0 || block_elapsed >= min_batch_time ||
                batch >= 10000L) break
            growth <- if (block_elapsed <= 0) {
              10
            } else {
              max(2, ceiling(1.1 * min_batch_time / block_elapsed))
            }
            batch <- min(10000L, as.integer(batch * growth))
          }
          elapsed <- block_elapsed / batch
          completed <<- completed + 1L
          rows[[completed]] <<- data.frame(
            experiment = experiment,
            n = n,
            k = k,
            rep = rep,
            method = method,
            time = elapsed,
            batch = batch,
            block_elapsed = block_elapsed,
            detected = length(boundaries)
          )
          utils::setTxtProgressBar(progress, completed)
        }
      }
    }
  }

  add_tasks(
    "vary_n", n_values, reps_n,
    n_for = function(value) value,
    k_for = function(value) 0L,
    seed_for = function(value, rep) seed + value + rep
  )
  add_tasks(
    "vary_k", k_values, reps_k,
    n_for = function(value) fixed_n,
    k_for = function(value) value,
    seed_for = function(value, rep) seed + 100000L + value * 100L + rep
  )
  dplyr::bind_rows(rows)
}

bootstrap_median_interval_time <- function(x, draws = 1000L) {
  if (!length(x)) return(c(NA_real_, NA_real_))
  values <- replicate(draws, stats::median(sample(x, replace = TRUE)))
  stats::quantile(values, c(0.025, 0.975), names = FALSE)
}

empty_time_summary <- function(group) {
  groups <- if (identical(group, "n")) {
    list(method = character(), n = integer())
  } else {
    list(method = character(), detected = integer())
  }
  data.frame(
    groups,
    replicates = integer(),
    mean_time = numeric(),
    median_time = numeric(),
    lower = numeric(),
    upper = numeric()
  )
}

summarise_time_by_n <- function(results, draws = 1000L, seed = 901L) {
  set.seed(seed)
  data <- results |>
    dplyr::filter(experiment == "vary_n", is.finite(time), time > 0) |>
    dplyr::ungroup()
  if (!nrow(data)) return(empty_time_summary("n"))
  data |>
    dplyr::group_by(method, n) |>
    dplyr::group_modify(function(data, key) {
      interval <- bootstrap_median_interval_time(data$time, draws)
      data.frame(
        replicates = nrow(data),
        mean_time = mean(data$time),
        median_time = stats::median(data$time),
        lower = interval[1L],
        upper = interval[2L]
      )
    }) |>
    dplyr::ungroup()
}

summarise_time_by_detected <- function(
    results, draws = 1000L, seed = 903L, min_replicates = 4L) {
  set.seed(seed)
  data <- results |>
    dplyr::filter(
      experiment == "vary_k", is.finite(time), time > 0,
      is.finite(detected)
    ) |>
    dplyr::ungroup()
  if (!nrow(data)) return(empty_time_summary("detected"))
  data <- data |>
    dplyr::group_by(method, detected) |>
    dplyr::mutate(group_replicates = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(group_replicates >= min_replicates)
  if (!nrow(data)) return(empty_time_summary("detected"))
  data |>
    dplyr::group_by(method, detected) |>
    dplyr::group_modify(function(data, key) {
      interval <- bootstrap_median_interval_time(data$time, draws)
      data.frame(
        replicates = nrow(data),
        mean_time = mean(data$time),
        median_time = stats::median(data$time),
        lower = interval[1L],
        upper = interval[2L]
      )
    }) |>
    dplyr::ungroup()
}

bootstrap_time_log_slope <- function(
    data, min_n = NULL, draws = 1000L) {
  n_values <- sort(unique(data$n))
  if (is.null(min_n)) min_n <- stats::median(n_values)
  n_values <- n_values[n_values >= min_n]
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
    samples, c(0.025, 0.975), names = FALSE, na.rm = TRUE
  )
  c(
    estimate = estimate(FALSE),
    lower = interval[1L],
    upper = interval[2L]
  )
}

summarise_time_slopes <- function(
    results, draws = 1000L, seed = 902L, min_n = NULL) {
  set.seed(seed)
  data <- dplyr::filter(
    results, experiment == "vary_n", is.finite(time), time > 0
  )
  if (is.null(min_n)) {
    min_n <- stats::median(unique(data$n))
  }
  dplyr::bind_rows(lapply(unique(data$method), function(method) {
    method_data <- data[data$method == method, , drop = FALSE]
    values <- bootstrap_time_log_slope(method_data, min_n, draws)
    data.frame(
      method = method,
      range = paste0("n >= ", min_n),
      exponent = values[["estimate"]],
      lower = values[["lower"]],
      upper = values[["upper"]]
    )
  }))
}

time_method_palette <- function(methods) {
  methods <- unique(as.character(methods))
  colours <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#000000"
  )
  stats::setNames(colours[seq_along(methods)], methods)
}

combine_time_panels <- function(panel_a, panel_b, base_size = 14L) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required for the combined runtime figure")
  }
  panel_b <- panel_b + ggplot2::guides(
    colour = "none", fill = "none"
  )
  combined <- patchwork::wrap_plots(
    panel_a, panel_b, ncol = 2L, guides = "collect"
  ) + patchwork::plot_annotation(
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")"
  )
  combined &
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      text = ggplot2::element_text(size = base_size),
      plot.tag = ggplot2::element_text(size = base_size, face = "bold")
    )
}

plot_paper_time_study <- function(
    results, root, paper_root, prefix, paper_prefix,
    slope_min_n = NULL, base_size = 14L, detected_methods = NULL,
    detected_max = NULL, slope_label_factors = NULL) {
  summary_n <- summarise_time_by_n(results)
  summary_detected <- summarise_time_by_detected(results)
  plot_summary_detected <- summary_detected
  if (!is.null(detected_methods)) {
    plot_summary_detected <- dplyr::filter(
      plot_summary_detected, method %in% detected_methods
    )
  }
  if (!is.null(detected_max)) {
    if (length(detected_max) != 1L || !is.finite(detected_max) ||
        detected_max < 1) {
      stop("detected_max must be one finite positive value")
    }
    plot_summary_detected <- dplyr::filter(
      plot_summary_detected, detected <= detected_max
    )
  }
  slopes <- summarise_time_slopes(results, min_n = slope_min_n)
  method_colours <- time_method_palette(summary_n$method)

  slope_labels <- summary_n[
    summary_n$n == max(summary_n$n), , drop = FALSE
  ]
  slope_labels <- merge(
    slope_labels, slopes[, c("method", "exponent")],
    by = "method", all.x = TRUE, sort = FALSE
  )
  slope_labels <- slope_labels[order(slope_labels$method), , drop = FALSE]
  short_method_label <- function(method) {
    method <- sub("^SVP MedianMood.*$", "Mood", method)
    method <- sub("^SVP Wilcoxon.*$", "Wilcoxon", method)
    method <- sub("^SVP AR1Focus.*$", "SVP AR1Focus", method)
    method
  }
  slope_labels$label <- paste0(
    short_method_label(slope_labels$method),
    ": ", sprintf("%.3f", slope_labels$exponent)
  )
  if (is.null(slope_label_factors)) {
    slope_label_factors <- 1.08 /
      1.6^(seq_len(nrow(slope_labels)) - 1L)
  }
  if (length(slope_label_factors) != nrow(slope_labels) ||
      any(!is.finite(slope_label_factors)) ||
      any(slope_label_factors <= 0)) {
    stop("slope_label_factors must have one positive finite value per method")
  }
  slope_labels$annotation_y <- max(summary_n$median_time) *
    slope_label_factors

  panel_a <- ggplot2::ggplot(
    summary_n, ggplot2::aes(n, median_time, colour = method, fill = method)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::geom_label(
      data = slope_labels,
      ggplot2::aes(label = label, y = annotation_y),
      hjust = 1.05, vjust = 0.5, fill = "white", alpha = 0.8,
      linewidth = 0, show.legend = FALSE
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Sequence length", y = "Median elapsed time (s)",
      colour = "Method", fill = "Method"
    ) +
    ggplot2::scale_colour_manual(values = method_colours) +
    ggplot2::scale_fill_manual(values = method_colours) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )

  panel_b <- ggplot2::ggplot(
    plot_summary_detected,
    ggplot2::aes(detected, median_time, colour = method, fill = method)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.12, colour = NA
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
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )
  if (!is.null(detected_max)) {
    panel_b <- panel_b + ggplot2::scale_x_continuous(
      limits = c(1, detected_max),
      breaks = unique(c(1, seq(20, detected_max, by = 20)))
    )
  }

  combined <- combine_time_panels(panel_a, panel_b, base_size = base_size)

  plot_root <- file.path(root, "plots")
  dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(paper_root, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(plot_root, paste0(prefix, "_time_vs_n.pdf")),
    panel_a, width = 9, height = 5.2
  )
  ggplot2::ggsave(
    file.path(plot_root, paste0(prefix, "_time_vs_detected.pdf")),
    panel_b, width = 9, height = 5.2
  )
  ggplot2::ggsave(
    file.path(paper_root, paste0(paper_prefix, "_time_vs_n.pdf")),
    panel_a, width = 9, height = 5.2
  )
  ggplot2::ggsave(
    file.path(paper_root, paste0(paper_prefix, "_time_vs_detected.pdf")),
    panel_b, width = 9, height = 5.2
  )
  ggplot2::ggsave(
    file.path(plot_root, paste0(prefix, "_time_pair.pdf")),
    combined, width = 13.2, height = 5.6
  )
  ggplot2::ggsave(
    file.path(paper_root, paste0(paper_prefix, "_time.pdf")),
    combined, width = 13.2, height = 5.6
  )
  list(
    panel_a = panel_a,
    panel_b = panel_b,
    combined = combined,
    summary_n = summary_n,
    summary_detected = summary_detected,
    slopes = slopes
  )
}

plot_time_study <- function(results, root, prefix) {
  dir.create(file.path(root, "plots"), recursive = TRUE, showWarnings = FALSE)
  vary_n <- dplyr::filter(results, experiment == "vary_n")
  vary_k <- dplyr::filter(results, experiment == "vary_k")
  positive_vary_n <- dplyr::filter(vary_n, is.finite(time), time > 0)
  positive_results <- dplyr::filter(results, is.finite(time), time > 0)

  time_vs_n <- ggplot2::ggplot(
    positive_vary_n, ggplot2::aes(n, time, colour = method)
  ) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Sequence length (log scale)", y = "Time (s, log scale)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  time_vs_k <- ggplot2::ggplot(
    vary_k, ggplot2::aes(k, time, colour = method)
  ) +
    ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::labs(x = "True number of changes", y = "Time (s)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  time_vs_detected <- ggplot2::ggplot(
    positive_results, ggplot2::aes(detected, time, colour = method)
  ) +
    ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Detected number of segments", y = "Time (s, log scale)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  plots <- list(
    time_vs_n = time_vs_n,
    time_vs_k = time_vs_k,
    time_vs_detected = time_vs_detected
  )
  for (name in names(plots)) {
    ggplot2::ggsave(
      file.path(root, "plots", paste0(prefix, "_", name, ".pdf")),
      plots[[name]], width = 8, height = 5
    )
  }
  invisible(plots)
}

save_time_outputs <- function(results, root, prefix) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(root, paste0(prefix, "_results.rds")))
  utils::write.csv(
    results, file.path(root, paste0(prefix, "_results.csv")),
    row.names = FALSE
  )
  plot_time_study(results, root, prefix)
  invisible(results)
}
