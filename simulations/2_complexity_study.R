#!/usr/bin/env Rscript
## simulations/complexity_study.R
## Benchmarking PELT vs SVP computational cost
## Two experiments:
## 1) No-change signals, huge penalty (effectively never detect change) for increasing n up to 10000
## 2) Fixed n = 10000, increase number of equally-spaced large changes from 1..100
##    SVP penalty = (1/avg_length_of_change) * log(n), PELT penalty = 2 * log(n)

library(future.apply)
library(future)
library(changepoint)
library(svpChange)
library(ggplot2)
library(dplyr)
library(tidyr)

## Utility: run algorithm and measure time + detected changepoints
time_and_count_pelt <- function(data, penalty) {
  res <- list(time = NA_real_, nchanges = NA_integer_)
  t <- system.time({
    obj <- tryCatch(cpt.mean(data, method = "PELT", penalty = "Manual", pen.value = penalty),
                    error = function(e) NULL)
  })[['elapsed']]
  res$time <- ifelse(is.finite(t), t, NA_real_)
  if (!is.null(obj)) {
    ch <- tryCatch(cpts(obj), error = function(e) NULL)
    res$nchanges <- if (!is.null(ch)) length(as.integer(unlist(ch))) else NA_integer_
    res$nchanges <- res$nchanges + 1 # including the last point at n 
  }
  res
}

time_and_count_svp <- function(data, penalty, cost = "gaussian_mean") {
  res <- list(time = NA_real_, nchanges = NA_integer_)
  t <- system.time({
    obj <- tryCatch(SVP(data, penalty, cost), error = function(e) NULL)
  })[['elapsed']]
  res$time <- ifelse(is.finite(t), t, NA_real_)
  if (!is.null(obj)) {
    ch <- NULL
    if (is.list(obj) && !is.null(obj$changepoints)) ch <- obj$changepoints
    if (is.null(ch)) ch <- tryCatch(obj$changepoints, error = function(e) NULL)
    res$nchanges <- if (!is.null(ch)) length(as.integer(unlist(ch))) else NA_integer_
  }
  res
}

## Experiment 1: no-change signals
run_exp1 <- function(n_values = round(seq(500, 10000, length.out = 20)),
                     reps = 20,
                     penalty_large = 1000,
                     seed = 123) {
  set.seed(seed)

  res_list <- future_lapply(n_values, function(n) {
    replicate(reps, {
      data <- rnorm(n, mean = 0, sd = 1)
      pelt <- time_and_count_pelt(data, penalty_large)
      svp  <- time_and_count_svp(data, penalty_large, cost = "gaussian_mean")
      list(n = n,
           PELT = pelt$time, PELT_nchanges = pelt$nchanges,
           SVP  = svp$time,  SVP_nchanges  = svp$nchanges)
    }, simplify = FALSE)
  }, future.seed = TRUE)

  results_flat <- do.call(c, res_list)
  df <- do.call(rbind, lapply(results_flat, as.data.frame))
  df$n <- as.integer(df$n)
  path <- file.path("simulations", "complexity_exp1_nochange.csv")
  write.csv(df, path, row.names = FALSE)
  message("Experiment 1 saved to: ", path)
  df
}

## Experiment 2: fixed n, vary number of equally-spaced large changes
run_exp2 <- function(n = 10000,
                     k_values = 1:100,
                     reps = 10,
                     gap = 10,
                     penalty_pelt = NULL,
                     seed = 456) {
  set.seed(seed)

  if (is.null(penalty_pelt)) penalty_pelt <- 2 * log(n)

  bench_one_k <- function(k) {
    # number of changepoints = k ; number of segments = k + 1
    seg_count <- k + 1
    base_seg_len <- floor(n / seg_count)
    seg_sizes <- rep(base_seg_len, seg_count)
    remainder <- n - sum(seg_sizes)
    if (remainder > 0) seg_sizes[1:remainder] <- seg_sizes[1:remainder] + 1

    #avg_len <- n / seg_count
    #penalty_svp <- 0.5 * log(avg_len)
    penalty_svp <- 2 * penalty_pelt

    replicate(reps, {
      # build piecewise-constant signal with large jumps
      means <- rep(c(0, gap), length.out = seg_count)
      x <- unlist(mapply(function(m, sz) rnorm(sz, mean = m, sd = 1), means, seg_sizes, SIMPLIFY = FALSE))

      svp <- time_and_count_svp(x, penalty_svp, cost = "gaussian_mean")
      pelt <- time_and_count_pelt(x, penalty_pelt)
      list(k = k, n = n,
           SVP = svp$time, SVP_nchanges = svp$nchanges,
           PELT = pelt$time, PELT_nchanges = pelt$nchanges,
           penalty_svp = penalty_svp, penalty_pelt = penalty_pelt)
    }, simplify = FALSE)
  }

  res_list <- future_lapply(k_values, bench_one_k, future.seed = TRUE)

  results_flat <- do.call(c, res_list)
  df <- do.call(rbind, lapply(results_flat, as.data.frame))
  df$k <- as.integer(df$k)
  path <- file.path("simulations", "complexity_exp2_vary_k.csv")
  write.csv(df, path, row.names = FALSE)
  message("Experiment 2 saved to: ", path)
  df
}



## Helper: summarise results (simple averages)
summarise_results <- function(df, by = c("n")) {
  library(dplyr)
  df %>% group_by(across(all_of(by))) %>%
    summarise(PELT_mean = mean(PELT, na.rm = TRUE),
              SVP_mean  = mean(SVP, na.rm = TRUE),
              .groups = "drop")
}

## Example: do not run both automatically; user can source and call

## Self-contained flags to run experiments and plot results
DO_RUN <- TRUE
DO_PLOT <- TRUE

if (DO_RUN) {
  ## Configure parallel plan so work is distributed across workers for all tasks
  RUN_WORKERS <- 4
  plan(multisession, workers = RUN_WORKERS)

  ## WARNING: these runs can be time-consuming. Adjust reps/k_values/workers as needed.
  message("[complexity_study] Running Experiment 1 (no-change scaling) ...")
  df1 <- run_exp1(n_values = round(2^(seq(10, 14, length.out = 20))), reps = 20)

  message("[complexity_study] Running Experiment 2 (varying k) ...")
  df2 <- run_exp2(n = 10000, k_values = 0:100, reps = 10, gap = 10)

  ## After runs, release plan (optional)
  plan(sequential)
}

if (DO_PLOT) {
  message("[complexity_study] Plotting results ...")
  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  # load data if not in memory
  if (!exists("df1")) {
    path1 <- file.path("simulations", "complexity_exp1_nochange.csv")
    if (file.exists(path1)) df1 <- read.csv(path1)
  }
  if (!exists("df2")) {
    path2 <- file.path("simulations", "complexity_exp2_vary_k.csv")
    if (file.exists(path2)) df2 <- read.csv(path2)
  }

  # Plot for Experiment 1: mean time vs n (log-log)
  if (exists("df1") && nrow(df1) > 0) {
    # Remove NA, NaN, or Inf values before fitting
    df1_clean <- df1[is.finite(df1$PELT) & is.finite(df1$SVP) & 
                      df1$PELT > 0 & df1$SVP > 0, ]
    
    # Calculate slopes for PELT and SVP on log-log scale
    if (nrow(df1_clean) > 2) {
      fit_pelt <- lm(log(PELT) ~ log(n), data = df1_clean)
      fit_svp <- lm(log(SVP) ~ log(n), data = df1_clean)
      slope_pelt <- round(coef(fit_pelt)[2], 3)
      slope_svp <- round(coef(fit_svp)[2], 3)
    } else {
      warning("Not enough valid data points to fit slopes")
      slope_pelt <- NA
      slope_svp <- NA
    }
    
    p1 <- ggplot(df1, aes(x = n, y = PELT, color = "PELT")) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2) +
      geom_smooth(data = df1, aes(x = n, y = SVP, color = "SVP"), 
                  method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2) +
      scale_x_log10() + scale_y_log10() +
      labs(x = "n (log)", y = "Time (s, log)") +
      theme_minimal()
    
    # Add slope annotations if valid
    if (!is.na(slope_pelt)) {
      p1 <- p1 + annotate("text", x = Inf, y = Inf, label = paste("PELT slope:", slope_pelt), 
                          hjust = 1.1, vjust = 20, color = "#F8766D", size = 4)
    }
    if (!is.na(slope_svp)) {
      p1 <- p1 + annotate("text", x = Inf, y = Inf, label = paste("SVP slope:", slope_svp), 
                          hjust = 1.1, vjust = 22, color = "#00BFC4", size = 4)
    }
    ggsave(filename = file.path("simulations", "plots", "exp1_time_vs_n.pdf"), plot = p1, width = 8, height = 5)
    message("[complexity_study] Saved plot: simulations/plots/exp1_time_vs_n.png")
  } else {
    message("[complexity_study] No data for Experiment 1 plotting (df1 missing).")
  }

  # Plot for Experiment 2: mean time vs k
  if (exists("df2") && nrow(df2) > 0) {
    p2 <- ggplot(df2, aes(x = k, y = PELT, color = "PELT")) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
      geom_smooth(data = df2, aes(x = k, y = SVP, color = "SVP"), 
                  method = "loess", se = TRUE, alpha = 0.2) +
      labs(x = "Number of changepoints (k)", y = "Time (s)") +
      theme_minimal()
    ggsave(filename = file.path("simulations", "plots", "exp2_time_vs_k.pdf"), plot = p2, width = 8, height = 5)
    message("[complexity_study] Saved plot: simulations/plots/exp2_time_vs_k.png")
  } else {
    message("[complexity_study] No data for Experiment 2 plotting (df2 missing).")
  }

  # New plot: runtime vs number of detected changes (combine both experiments)
  if ((exists("df1") && nrow(df1) > 0) || (exists("df2") && nrow(df2) > 0)) {
    build_long <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      rows <- lapply(seq_len(nrow(df)), function(i) {
        r <- df[i, , drop = FALSE]
        data.frame(true_k = if ("k" %in% names(r)) r$k else NA_integer_,
                   method = c("PELT", "SVP"),
                   time = c(r$PELT, r$SVP),
                   detected = c(r$PELT_nchanges, r$SVP_nchanges))
      })
      do.call(rbind, rows)
    }

    long1 <- build_long(df1)
    long2 <- build_long(df2)
    combined <- do.call(rbind, Filter(Negate(is.null), list(long1, long2)))

    # ensure detected is numeric
    combined$detected <- as.integer(combined$detected)
    combined$time <- as.numeric(combined$time)

    p3 <- ggplot(combined, aes(x = detected, y = time, color = method, fill = method)) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
      scale_y_log10() +
      labs(x = "Detected number of changes",
           y = "Time (s, log)") +
      theme_minimal()
    ggsave(filename = file.path("simulations", "plots", "time_vs_detected_changes.pdf"), plot = p3, width = 8, height = 5)
    message("[complexity_study] Saved plot: simulations/plots/time_vs_detected_changes.png")
  } else {
    message("[complexity_study] No data available to plot runtime vs detected changes.")
  }
}

## End of file
