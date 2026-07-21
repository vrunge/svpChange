## Computational cost of SVP with Gaussian and robust cost tests.
##
## This follows the same simulation structure as the original Figure 4-style
## complexity benchmark:
## 1) no-change signals, increasing n;
## 2) fixed n, increasing number of equally-spaced changes.
##
## The comparison is restricted to SVP implementations so the extra cost of the
## robust Wilcoxon and Median-Mood validity functions is visible.

library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(tidyr)
library(svpChange)

################################################################################

RUN_WORKERS <- 4

DO_RUN <- FALSE
DO_PLOT <- TRUE

EXP1_N_VALUES <- round(2^(seq(8, 12, length.out = 12)))
EXP1_REPS <- 10

EXP2_N <- 2000
EXP2_K_VALUES <- 0:40
EXP2_REPS <- 5

OUTPUT_EXP1 <- file.path("simulations", "complexity_allCost_exp1_nochange.csv")
OUTPUT_EXP2 <- file.path("simulations", "complexity_allCost_exp2_vary_k.csv")

PLOT_EXP1 <- file.path("simulations", "plots", "2_allCost_time_vs_n.pdf")
PLOT_EXP2 <- file.path("simulations", "plots", "2_allCost_time_vs_k.pdf")
PLOT_COMBINED <- file.path("simulations", "plots", "2_allCost_complexity_figure.pdf")

SVP_METHODS <- c("SVP Gaussian", "SVP Wilcoxon", "SVP Median-Mood")

################################################################################

threshold_gaussian <- function(n)
{
  2 * log(n)
}

threshold_wilcoxon <- function(n, nbSeg)
{
  1.5 * sqrt((n / nbSeg)^3 / 12)
}

threshold_mood <- function(n, nbSeg, alpha = 0.01)
{
  m <- max(n / nbSeg - 1, 1)
  p_single <- 1 - (1 - alpha)^(1 / m)
  qchisq(1 - p_single, df = max(nbSeg, 1))
}

################################################################################

generate_piecewise_signal <- function(n,
                                      k,
                                      gap = 10,
                                      noise = c("gaussian", "student"))
{
  noise <- match.arg(noise)
  seg_count <- k + 1
  base_seg_len <- floor(n / seg_count)
  seg_sizes <- rep(base_seg_len, seg_count)
  remainder <- n - sum(seg_sizes)

  if (remainder > 0) {
    seg_sizes[seq_len(remainder)] <- seg_sizes[seq_len(remainder)] + 1
  }

  means <- rep(c(0, gap), length.out = seg_count)
  signal <- unlist(
    mapply(
      function(m, sz) rep(m, sz),
      means,
      seg_sizes,
      SIMPLIFY = FALSE
    )
  )

  eps <- switch(
    noise,
    gaussian = rnorm(n, mean = 0, sd = 1),
    student = rt(n, df = 2)
  )

  signal + eps
}

################################################################################

time_one_svp <- function(data,
                         method = SVP_METHODS,
                         nbSeg = 1)
{
  method <- match.arg(method)
  n <- length(data)

  obj <- NULL
  err <- NULL

  elapsed <- system.time({
    obj <- tryCatch(
      switch(
        method,
        "SVP Gaussian" = SVP(
          data,
          gamma = threshold_gaussian(n),
          test = "gaussian_mean"
        ),
        "SVP Wilcoxon" = SVP_costTests(
          data,
          gamma = threshold_wilcoxon(n, nbSeg),
          test = "WilcoxonCost"
        ),
        "SVP Median-Mood" = SVP_costTests(
          data,
          gamma = threshold_mood(n, nbSeg),
          test = "MedianMoodCost"
        )
      ),
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
  })[["elapsed"]]

  changepoints <- integer(0)
  mean_candidates <- NA_real_
  max_candidates <- NA_integer_

  if (!is.null(obj)) {
    if (!is.null(obj$changepoints)) {
      changepoints <- as.integer(unlist(obj$changepoints))
    }
    if (!is.null(obj$nb)) {
      mean_candidates <- mean(as.numeric(obj$nb), na.rm = TRUE)
      max_candidates <- max(as.integer(obj$nb), na.rm = TRUE)
    }
  }

  data.frame(
    method = method,
    time = ifelse(is.finite(elapsed), elapsed, NA_real_),
    nsegments = length(changepoints),
    mean_candidates = mean_candidates,
    max_candidates = max_candidates,
    error = ifelse(is.null(err), NA_character_, err)
  )
}

################################################################################

run_exp1 <- function(n_values = EXP1_N_VALUES,
                     reps = EXP1_REPS,
                     methods = SVP_METHODS,
                     noise = "student",
                     seed = 123)
{
  tasks <- expand.grid(
    n = n_values,
    rep = seq_len(reps),
    method = methods,
    stringsAsFactors = FALSE
  )

  res_list <- future_lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, , drop = FALSE]
    set.seed(seed + i)

    x <- generate_piecewise_signal(
      n = task$n,
      k = 0,
      noise = noise
    )
    out <- time_one_svp(x, method = task$method, nbSeg = 1)

    cbind(
      experiment = "no_change",
      n = task$n,
      k = 0L,
      rep = task$rep,
      noise = noise,
      out
    )
  }, future.seed = TRUE)

  df <- bind_rows(res_list)
  write.csv(df, OUTPUT_EXP1, row.names = FALSE)
  message("[complexity_allCost] Experiment 1 saved to: ", OUTPUT_EXP1)
  df
}

################################################################################

run_exp2 <- function(n = EXP2_N,
                     k_values = EXP2_K_VALUES,
                     reps = EXP2_REPS,
                     methods = SVP_METHODS,
                     noise = "student",
                     gap = 10,
                     seed = 456)
{
  tasks <- expand.grid(
    k = k_values,
    rep = seq_len(reps),
    method = methods,
    stringsAsFactors = FALSE
  )

  res_list <- future_lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, , drop = FALSE]
    set.seed(seed + i)

    x <- generate_piecewise_signal(
      n = n,
      k = task$k,
      gap = gap,
      noise = noise
    )
    nbSeg <- task$k + 1
    out <- time_one_svp(x, method = task$method, nbSeg = nbSeg)

    cbind(
      experiment = "vary_k",
      n = n,
      k = task$k,
      rep = task$rep,
      noise = noise,
      out
    )
  }, future.seed = TRUE)

  df <- bind_rows(res_list)
  write.csv(df, OUTPUT_EXP2, row.names = FALSE)
  message("[complexity_allCost] Experiment 2 saved to: ", OUTPUT_EXP2)
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
    mutate(label = paste0(method, ": ", round(slope, 2)))
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
      color = "Implementation",
      fill = "Implementation"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste(slopes$label, collapse = "\n"),
      hjust = 1.05,
      vjust = 1.15,
      size = 3.1
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
    geom_point(size = 1.4) +
    scale_y_log10() +
    labs(
      x = "Number of true changepoints",
      y = "Mean runtime in seconds (log scale)",
      color = "Implementation",
      fill = "Implementation"
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
      color = "Implementation",
      fill = "Implementation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

################################################################################

run_all_cost_complexity_study <- function(workers = RUN_WORKERS)
{
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)

  message("[complexity_allCost] Running Experiment 1 ...")
  df1 <- run_exp1()

  message("[complexity_allCost] Running Experiment 2 ...")
  df2 <- run_exp2()

  invisible(list(exp1 = df1, exp2 = df2))
}

save_all_cost_plots <- function(df1 = NULL,
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
    message("[complexity_allCost] Saved plot: ", PLOT_EXP1)
  }

  if (!is.null(df2) && nrow(df2) > 0) {
    p2 <- plot_exp2(df2)
    ggsave(filename = PLOT_EXP2, plot = p2, width = 8, height = 5)
    message("[complexity_allCost] Saved plot: ", PLOT_EXP2)
  }

  if (!is.null(df1) && nrow(df1) > 0 && !is.null(df2) && nrow(df2) > 0) {
    p3 <- plot_combined(df1, df2)
    ggsave(filename = PLOT_COMBINED, plot = p3, width = 10, height = 5)
    message("[complexity_allCost] Saved plot: ", PLOT_COMBINED)
  }

  invisible(TRUE)
}

################################################################################

if (DO_RUN)
{
  res <- run_all_cost_complexity_study()
  df1 <- res$exp1
  df2 <- res$exp2
}

if (DO_PLOT)
{
  save_all_cost_plots(
    df1 = if (exists("df1")) df1 else NULL,
    df2 = if (exists("df2")) df2 else NULL
  )
}

## End of file
