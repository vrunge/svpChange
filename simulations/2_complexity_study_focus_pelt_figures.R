## Figure 4-style runtime comparison for PELT and the old SVP FOCuS implementation.
##
## This is a non-overwriting companion to 2_complexity_study.R. It uses the same
## two simulation settings but writes separate CSVs and figures.

library(future)
library(future.apply)
library(changepoint)
library(svpChange)
library(ggplot2)
library(dplyr)
library(tidyr)

################################################################################

RUN_WORKERS <- 4

EXP1_N_VALUES <- round(2^(seq(10, 14, length.out = 20)))
EXP1_REPS <- 20

EXP2_N <- 10000
EXP2_K_VALUES <- 0:100
EXP2_REPS <- 10

OUTPUT_EXP1 <- file.path("simulations", "complexity_focus_pelt_exp1_nochange.csv")
OUTPUT_EXP2 <- file.path("simulations", "complexity_focus_pelt_exp2_vary_k.csv")

PLOT_EXP1 <- file.path("simulations", "plots", "2_focus_pelt_time_vs_n.pdf")
PLOT_EXP2 <- file.path("simulations", "plots", "2_focus_pelt_time_vs_k.pdf")
PLOT_DETECTED <- file.path("simulations", "plots", "2_focus_pelt_time_vs_detected_changes.pdf")
PLOT_COMBINED <- file.path("simulations", "plots", "2_focus_pelt_complexity_figure.pdf")

################################################################################

time_and_count_pelt <- function(data, penalty)
{
  obj <- NULL
  elapsed <- system.time({
    obj <- tryCatch(
      cpt.mean(data, method = "PELT", penalty = "Manual", pen.value = penalty),
      error = function(e) NULL
    )
  })[["elapsed"]]

  nchanges <- NA_integer_
  if (!is.null(obj)) {
    ch <- tryCatch(cpts(obj), error = function(e) NULL)
    nchanges <- if (!is.null(ch)) length(as.integer(unlist(ch))) + 1L else NA_integer_
  }

  data.frame(
    method = "PELT",
    time = ifelse(is.finite(elapsed), elapsed, NA_real_),
    detected = nchanges
  )
}

time_and_count_svp_focus <- function(data, gamma)
{
  obj <- NULL
  elapsed <- system.time({
    obj <- tryCatch(
      SVP_old(data, gamma = gamma, test = "gaussian_mean"),
      error = function(e) NULL
    )
  })[["elapsed"]]

  nchanges <- NA_integer_
  if (!is.null(obj) && !is.null(obj$changepoints)) {
    nchanges <- length(as.integer(unlist(obj$changepoints)))
  }

  data.frame(
    method = "SVP FOCuS",
    time = ifelse(is.finite(elapsed), elapsed, NA_real_),
    detected = nchanges
  )
}

################################################################################

generate_piecewise_signal <- function(n, k, gap = 10)
{
  seg_count <- k + 1
  base_seg_len <- floor(n / seg_count)
  seg_sizes <- rep(base_seg_len, seg_count)
  remainder <- n - sum(seg_sizes)

  if (remainder > 0) {
    seg_sizes[seq_len(remainder)] <- seg_sizes[seq_len(remainder)] + 1
  }

  means <- rep(c(0, gap), length.out = seg_count)
  unlist(
    mapply(
      function(m, sz) rnorm(sz, mean = m, sd = 1),
      means,
      seg_sizes,
      SIMPLIFY = FALSE
    )
  )
}

################################################################################

run_exp1 <- function(n_values = EXP1_N_VALUES,
                     reps = EXP1_REPS,
                     penalty_large = 1000,
                     seed = 123)
{
  tasks <- expand.grid(
    n = n_values,
    rep = seq_len(reps),
    stringsAsFactors = FALSE
  )

  res_list <- future_lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, , drop = FALSE]
    set.seed(seed + i)

    x <- rnorm(task$n, mean = 0, sd = 1)

    bind_rows(
      time_and_count_pelt(x, penalty = penalty_large),
      time_and_count_svp_focus(x, gamma = penalty_large)
    ) %>%
      mutate(
        experiment = "no_change",
        n = task$n,
        k = 0L,
        rep = task$rep,
        .before = method
      )
  }, future.seed = TRUE)

  df <- bind_rows(res_list)
  write.csv(df, OUTPUT_EXP1, row.names = FALSE)
  message("[complexity_focus_pelt] Experiment 1 saved to: ", OUTPUT_EXP1)
  df
}

################################################################################

run_exp2 <- function(n = EXP2_N,
                     k_values = EXP2_K_VALUES,
                     reps = EXP2_REPS,
                     gap = 10,
                     penalty_pelt = NULL,
                     seed = 456)
{
  if (is.null(penalty_pelt)) {
    penalty_pelt <- 2 * log(n)
  }

  tasks <- expand.grid(
    k = k_values,
    rep = seq_len(reps),
    stringsAsFactors = FALSE
  )

  res_list <- future_lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, , drop = FALSE]
    set.seed(seed + i)

    x <- generate_piecewise_signal(n = n, k = task$k, gap = gap)
    penalty_svp <- 2 * penalty_pelt

    bind_rows(
      time_and_count_pelt(x, penalty = penalty_pelt),
      time_and_count_svp_focus(x, gamma = penalty_svp)
    ) %>%
      mutate(
        experiment = "vary_k",
        n = n,
        k = task$k,
        rep = task$rep,
        penalty_svp = penalty_svp,
        penalty_pelt = penalty_pelt,
        .before = method
      )
  }, future.seed = TRUE)

  df <- bind_rows(res_list)
  write.csv(df, OUTPUT_EXP2, row.names = FALSE)
  message("[complexity_focus_pelt] Experiment 2 saved to: ", OUTPUT_EXP2)
  df
}

################################################################################

summarise_runtime <- function(df, xvar)
{
  df %>%
    filter(is.finite(time), time > 0) %>%
    group_by(method, .data[[xvar]]) %>%
    summarise(
      mean_time = mean(time, na.rm = TRUE),
      sd_time = sd(time, na.rm = TRUE),
      n_rep = n(),
      se_time = ifelse(is.na(sd_time), 0, sd_time / sqrt(pmax(1, n_rep))),
      lower = pmax(mean_time - 1.96 * se_time, .Machine$double.eps),
      upper = mean_time + 1.96 * se_time,
      .groups = "drop"
    )
}

estimate_loglog_slopes <- function(df, xvar)
{
  df %>%
    filter(is.finite(time), time > 0, .data[[xvar]] > 0) %>%
    group_by(method) %>%
    summarise(
      slope = if (n() > 2) coef(lm(log(time) ~ log(.data[[xvar]])))[2] else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(label = paste0(method, ": ", round(slope, 3)))
}

plot_exp1 <- function(df1)
{
  stats <- summarise_runtime(df1, "n")
  slopes <- estimate_loglog_slopes(df1, "n")

  ggplot(stats, aes(x = n, y = mean_time, color = method, fill = method)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Signal length n (log scale)",
      y = "Mean runtime in seconds (log scale)",
      color = "Method",
      fill = "Method"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste(slopes$label, collapse = "\n"),
      hjust = 1.05,
      vjust = 1.15,
      size = 3.4
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_exp2 <- function(df2)
{
  stats <- summarise_runtime(df2, "k")

  ggplot(stats, aes(x = k, y = mean_time, color = method, fill = method)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_y_log10() +
    labs(
      x = "Number of true changepoints",
      y = "Mean runtime in seconds (log scale)",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_detected <- function(df1, df2)
{
  bind_rows(df1, df2) %>%
    filter(is.finite(time), time > 0, !is.na(detected)) %>%
    ggplot(aes(x = detected, y = time, color = method, fill = method)) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15) +
    scale_y_log10() +
    labs(
      x = "Detected number of segments",
      y = "Runtime in seconds (log scale)",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_combined <- function(df1, df2)
{
  stats1 <- summarise_runtime(df1, "n") %>%
    mutate(panel = "A. No-change signals: increasing n", x = n)
  stats2 <- summarise_runtime(df2, "k") %>%
    mutate(panel = "B. Fixed n: increasing number of changes", x = k)

  bind_rows(stats1, stats2) %>%
    ggplot(aes(x = x, y = mean_time, color = method, fill = method)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_y_log10() +
    facet_wrap(~ panel, scales = "free_x", nrow = 1) +
    labs(
      x = NULL,
      y = "Mean runtime in seconds (log scale)",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

################################################################################

run_focus_pelt_complexity_study <- function(workers = RUN_WORKERS)
{
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)

  message("[complexity_focus_pelt] Running Experiment 1 ...")
  df1 <- run_exp1()

  message("[complexity_focus_pelt] Running Experiment 2 ...")
  df2 <- run_exp2()

  invisible(list(exp1 = df1, exp2 = df2))
}

save_focus_pelt_plots <- function(df1 = NULL,
                                  df2 = NULL)
{
  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  if (is.null(df1) && file.exists(OUTPUT_EXP1)) {
    df1 <- read.csv(OUTPUT_EXP1)
  }
  if (is.null(df2) && file.exists(OUTPUT_EXP2)) {
    df2 <- read.csv(OUTPUT_EXP2)
  }

  if (!is.null(df1) && nrow(df1) > 0) {
    p1 <- plot_exp1(df1)
    ggsave(filename = PLOT_EXP1, plot = p1, width = 8, height = 5)
    message("[complexity_focus_pelt] Saved plot: ", PLOT_EXP1)
  }

  if (!is.null(df2) && nrow(df2) > 0) {
    p2 <- plot_exp2(df2)
    ggsave(filename = PLOT_EXP2, plot = p2, width = 8, height = 5)
    message("[complexity_focus_pelt] Saved plot: ", PLOT_EXP2)
  }

  if (!is.null(df1) && nrow(df1) > 0 && !is.null(df2) && nrow(df2) > 0) {
    p3 <- plot_detected(df1, df2)
    ggsave(filename = PLOT_DETECTED, plot = p3, width = 8, height = 5)
    message("[complexity_focus_pelt] Saved plot: ", PLOT_DETECTED)

    p4 <- plot_combined(df1, df2)
    ggsave(filename = PLOT_COMBINED, plot = p4, width = 10, height = 5)
    message("[complexity_focus_pelt] Saved plot: ", PLOT_COMBINED)
  }

  invisible(TRUE)
}

## End of file
