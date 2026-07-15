
## simulations/power_study_robust_heavy_tails.R
## Power study comparing PELT vs SVP (Wilcoxon) vs SVP (MedianMood) vs rfpop on heavy-tailed noise.
## - n = 1000
## - patterns = c("none","up","updown","rand1")
## - jumpSize in seq(0.1, 2, by = 0.1)
## - reps per combination (default 100)
## Produces CSV and plots faceted by scenario.

library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(changepoint)
library(svpChange)
library(robseg) #devtools::install_github("guillemr/robust-fpop", force = TRUE)
library(progressr)


################################################################################

## generate_signal (copied/adapted from original study)
generate_signal <- function(n,
                            pattern = c("none", "up", "updown", "rand1"),
                            nbSeg = 8,
                            jumpSize = 1)
{
  type <- match.arg(pattern)

  if (type == "rand1") {
    set.seed(42)
    rand1CP <- rpois(nbSeg, lambda = 10)

    # Scale counts to total n
    r1 <- pmax(round(rand1CP * n / sum(rand1CP)), 1)

    # Adjust to sum exactly n
    diff <- n - sum(r1)
    r1[nbSeg] <- r1[nbSeg] + diff
    stopifnot(sum(r1) == n)

    set.seed(43)
    rand1Jump <- runif(nbSeg, min = 0.5, max = 1) * sample(c(-1, 1), nbSeg, replace = TRUE)
    set.seed(NULL)
  }

  switch(
    type,
    none = rep(0, n),
    up = rep(seq(0, nbSeg - 1) * jumpSize, each = n / nbSeg),
    updown = rep((seq(0, nbSeg - 1) %% 2) * jumpSize, each = n / nbSeg),
    rand1 = map2(rand1Jump, r1, ~rep(.x * jumpSize, .y)) %>% unlist()
  )
}


################################################################################

## cp_metrics: greedy 1-1 matching within tolerance
cp_metrics <- function(cp_true,
                       cp_est,
                       tol = 5)
{
  cp_true <- as.integer(unlist(cp_true))
  cp_est  <- as.integer(unlist(cp_est))

  if (length(cp_true) == 0 && length(cp_est) == 0) return(list(Precision = NA_real_, Recall = NA_real_, F1 = NA_real_, LocalizationError = NA_real_))
  if (length(cp_true) == 0) return(list(Precision = 0, Recall = NA_real_, F1 = NA_real_, LocalizationError = NA_real_))
  if (length(cp_est) == 0) return(list(Precision = NA_real_, Recall = 0, F1 = NA_real_, LocalizationError = NA_real_))

  D <- abs(outer(cp_true, cp_est, "-"))
  matched_est <- rep(FALSE, length(cp_est))
  matched_dist <- numeric(0)
  TP <- 0
  for (i in seq_along(cp_true))
  {
    j <- which.min(D[i, ])
    if (length(j) && D[i, j] <= tol && !matched_est[j])
    {
      TP <- TP + 1
      matched_est[j] <- TRUE
      matched_dist <- c(matched_dist, D[i, j])
    }
  }
  FP <- sum(!matched_est)
  FN <- length(cp_true) - TP

  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
  loc_err <- if (length(cp_est) == length(cp_true)) max(matched_dist) else NA_real_
  list(Precision = precision, Recall = recall, F1 = f1, LocalizationError = loc_err)
}


################################################################################

## Threshold function for Median-Mood test (as provided)
threshold_mood <- function(n, nbSeg, alpha = 0.01)
{
  m <- n / nbSeg - 1
  p_single <- 1 - (1 - alpha)^(1 / m)  # Dunn–Šidák
  qchisq(1 - p_single, df = nbSeg)
}


################################################################################

run_robust_power_study <- function(n = 1000,
                                   patterns = c("none", "up", "updown", "rand1"),
                                   jumpSizes = seq(0.1, 2, by = 0.1),
                                   reps = 100,
                                   workers = max(1, future::availableCores() - 1),
                                   seed = 123,
                                   nbSeg = 8) {
  set.seed(seed)
  plan(multisession, workers = workers)

  combos <- expand.grid(pattern = patterns, jump = jumpSizes, stringsAsFactors = FALSE)

  # for each combo, do `reps` replicates in parallel across combos
  handlers("txtprogressbar")
  with_progress({
    # one tick per replicate (combinations * reps)
    p <- progressor(steps = nrow(combos) * reps)
    res_list <- future_lapply(seq_len(nrow(combos)), function(idx) {
      pattern <- combos$pattern[idx]
      jump <- combos$jump[idx]
      replicate(reps, {
        p()
        mu <- generate_signal(n, pattern = pattern, nbSeg = nbSeg, jumpSize = jump)
        # Heavy-tailed noise: t-distribution with df=2
        y <- mu + rt(n, df = 2)
        cp_true <- c(which(diff(mu) != 0), length(mu))
        nseg_true <- length(cp_true)
        # infer the actual number of segments in this scenario
        nbSeg_current <- length(cp_true)

        # penalties
        penalty_pelt <- 2 * log(n)

        # run PELT
        obj_pelt <- tryCatch(cpt.mean(y, method = "PELT", penalty = "Manual", pen.value = penalty_pelt), error = function(e) NULL)
        cp_pelt <- if (!is.null(obj_pelt)) as.integer(cpts(obj_pelt)) else integer(0)
        # ensure we include the final index as a segment boundary (consistent with true cps)
        cp_pelt <- unique(c(cp_pelt, length(y)))
        metrics_pelt <- cp_metrics(cp_true, cp_pelt, tol = round(n * 0.0025))
        # compute number of segments and MSE for PELT
        nseg_pelt <- length(cp_pelt)
        mse_pelt <- NA_real_
        if (!is.null(obj_pelt)) {
          means_pelt <- tryCatch(param.est(obj_pelt)$mean, error = function(e) NULL)
          if (!is.null(means_pelt)) {
            mu_hat_pelt <- rep(means_pelt, times = diff(c(0, cp_pelt)))
            if (length(mu_hat_pelt) == length(mu)) mse_pelt <- mean((mu - mu_hat_pelt)^2)
          }
        }

        # run SVP with Wilcoxon cost test using SVP_costTests
        # gamma uses the inferred nbSeg_current through cp_true length already below
        resW <- tryCatch(SVP_costTests(y, gamma = 1.5 * sqrt((n / nbSeg_current)^3 / 12), test = "WilcoxonCost"), error = function(e) NULL)
        cp_svp_wilk <- integer(0)
        if (!is.null(resW)) {
          if (is.list(resW) && !is.null(resW$changepoints)) cp_svp_wilk <- as.integer(unlist(resW$changepoints))
        }
        cp_svp_wilk <- unique(c(cp_svp_wilk, length(y)))
        metrics_svp_wilk <- cp_metrics(cp_true, cp_svp_wilk, tol = round(n * 0.0025))
        nseg_svp_wilk <- length(cp_svp_wilk)
        mse_svp_wilk <- NA_real_
        if (length(cp_svp_wilk) > 0) {
          mu_hat_svp_wilk <- numeric(0)
          for (i in seq_along(cp_svp_wilk)) {
            if (i == 1) {
              mu_hat_svp_wilk <- rep(mean(y[1:cp_svp_wilk[i]]), cp_svp_wilk[i])
            } else {
              mu_hat_svp_wilk <- c(mu_hat_svp_wilk, rep(mean(y[(cp_svp_wilk[i - 1] + 1):cp_svp_wilk[i]]), cp_svp_wilk[i] - cp_svp_wilk[i - 1]))
            }
          }
          if (length(mu_hat_svp_wilk) == length(mu)) mse_svp_wilk <- mean((mu - mu_hat_svp_wilk)^2)
        }

        # run SVP with Median-Mood Cost test using the threshold function provided
        # Use nbSeg_current (inferred) for the threshold calculation
        gamma_mood <- tryCatch(threshold_mood(n = n, nbSeg = nbSeg_current, alpha = 0.01), error = function(e) NULL)
        resM <- NULL
        if (!is.null(gamma_mood)) {
          resM <- tryCatch(SVP_costTests(y, gamma = gamma_mood, test = "MedianMoodCost"), error = function(e) NULL)
        }
        cp_svp_mood <- integer(0)
        if (!is.null(resM)) {
          if (is.list(resM) && !is.null(resM$changepoints)) cp_svp_mood <- as.integer(unlist(resM$changepoints))
        }
        cp_svp_mood <- unique(c(cp_svp_mood, length(y)))
        metrics_svp_mood <- cp_metrics(cp_true, cp_svp_mood, tol = round(n * 0.0025))
        nseg_svp_mood <- length(cp_svp_mood)
        mse_svp_mood <- NA_real_
        if (length(cp_svp_mood) > 0) {
          mu_hat_svp_mood <- numeric(0)
          for (i in seq_along(cp_svp_mood)) {
            if (i == 1) {
              mu_hat_svp_mood <- rep(mean(y[1:cp_svp_mood[i]]), cp_svp_mood[i])
            } else {
              mu_hat_svp_mood <- c(mu_hat_svp_mood, rep(mean(y[(cp_svp_mood[i - 1] + 1):cp_svp_mood[i]]), cp_svp_mood[i] - cp_svp_mood[i - 1]))
            }
          }
          if (length(mu_hat_svp_mood) == length(mu)) mse_svp_mood <- mean((mu - mu_hat_svp_mood)^2)
        }

        # run rfpop from robseg
        # Estimate standard deviation robustly using MAD
        est_sd <- mad(y, constant = 1.4826)
        res_rfpop <- tryCatch({
          Rob_seg.std(x = y / est_sd,
                      loss = "Outlier",
                      lambda = 2 * log(length(y)),
                      lthreshold = 3)
        }, error = function(e) NULL)

        cp_rfpop <- integer(0)
        if (!is.null(res_rfpop) && !is.null(res_rfpop$t.est)) {
          # Remove the final boundary from t.est
          cp_rfpop <- as.integer(res_rfpop$t.est[-length(res_rfpop$t.est)])
        }
        # ensure final index included for consistency with cp_true
        cp_rfpop <- unique(c(cp_rfpop, length(y)))
        metrics_rfpop <- cp_metrics(cp_true, cp_rfpop, tol = round(n * 0.0025))
        nseg_rfpop <- length(cp_rfpop)
        mse_rfpop <- NA_real_
        if (length(cp_rfpop) > 0 && !is.null(res_rfpop) && !is.null(res_rfpop$smt)) {
          # Use the smoothed profile from rfpop, rescale back
          mu_hat_rfpop <- res_rfpop$smt * est_sd
          if (length(mu_hat_rfpop) == length(mu)) mse_rfpop <- mean((mu - mu_hat_rfpop)^2)
        }

        tib <- tibble::tibble(
          pattern = pattern,
          jump = jump,
          rep = as.integer(NA),
          algorithm = c("PELT", "SVP (Wilcoxon)", "SVP (MedianMood)", "rfpop"),
          Precision = c(metrics_pelt$Precision, metrics_svp_wilk$Precision, metrics_svp_mood$Precision, metrics_rfpop$Precision),
          Recall = c(metrics_pelt$Recall, metrics_svp_wilk$Recall, metrics_svp_mood$Recall, metrics_rfpop$Recall),
          F1 = c(metrics_pelt$F1, metrics_svp_wilk$F1, metrics_svp_mood$F1, metrics_rfpop$F1),
          NumSegments = c(nseg_pelt, nseg_svp_wilk, nseg_svp_mood, nseg_rfpop),
          MSE = c(mse_pelt, mse_svp_wilk, mse_svp_mood, mse_rfpop),
          LocalizationError = c(metrics_pelt$LocalizationError, metrics_svp_wilk$LocalizationError, metrics_svp_mood$LocalizationError, metrics_rfpop$LocalizationError),
          PropCorrectNumCP = c(nseg_pelt == nseg_true, nseg_svp_wilk == nseg_true, nseg_svp_mood == nseg_true, nseg_rfpop == nseg_true),
          changepoints = list(cp_pelt, cp_svp_wilk, cp_svp_mood, cp_rfpop)
        )
        tib
      }, simplify = FALSE)
    }, future.seed = TRUE)
  })

  plan(sequential)

  # flatten: produce one data.frame from the nested list of tibbles
  all_reps <- unlist(res_list, recursive = FALSE)
  results_df <- dplyr::bind_rows(all_reps)

  # Add rep id
  results_df <- results_df %>% group_by(pattern, jump, algorithm) %>% mutate(rep = row_number()) %>% ungroup()

  # Save a version without the list-column for inspection (drop changepoints for CSV export)
  results_csv <- results_df %>% select(-changepoints)
  out_path <- file.path("simulations", "power_study_robust_heavy_tails_results.csv")
  write.csv(results_csv, out_path, row.names = FALSE)
  # Also save full RDS (preserves changepoints list-column)
  out_rds <- file.path("simulations", "power_study_robust_heavy_tails_results.rds")
  saveRDS(results_df, out_rds)
  message("Saved results to CSV: ", out_path)
  message("Saved full results (RDS) with changepoints to: ", out_rds)

  # Return the full results_df with changepoints intact for use in plotting
  results_df
}

################################################################################

plot_robust_power_metrics <- function(results_df)
{

  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  summarise_ci <- function(d)
  {
    d %>% summarise(mean = mean(value, na.rm = TRUE),
                    sd = sd(value, na.rm = TRUE),
                    n = sum(!is.na(value)),
                    se = sd / sqrt(pmax(1, n)),
                    lower = mean - 1.96 * se,
                    upper = mean + 1.96 * se,
                    .groups = "drop")
  }
  long <- results_df %>%
    pivot_longer(cols = c(Precision, Recall, F1, NumSegments, MSE, LocalizationError, PropCorrectNumCP),
                 names_to = "metric", values_to = "value")
    # Ensure pattern is a factor with consistent levels
  long$pattern <- factor(long$pattern, levels = c("none", "up", "updown", "rand1"))
  stats_df <- long %>%
    dplyr::group_by(algorithm, pattern, jump, metric) %>%
    dplyr::summarise(
      mean  = mean(value, na.rm = TRUE),
      sd    = sd(value, na.rm = TRUE),
      n     = sum(!is.na(value)),
      se    = sd / sqrt(pmax(1, n)),
      lower = mean - 1.96 * se,
      upper = mean + 1.96 * se,
      .groups = "drop"
    )

  stopifnot("metric" %in% names(stats_df))

  metrics <- c("F1", "Precision", "Recall", "NumSegments", "MSE", "LocalizationError", "PropCorrectNumCP")
  for (m in metrics) {
    dat <- dplyr::filter(stats_df, metric == m)

    p <- ggplot(dat, aes(x = jump, y = mean, color = algorithm, fill = algorithm)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_wrap(~factor(pattern, levels = c("none", "up", "updown", "rand1")), scales = "fixed") +
      labs(x = "Jump size", y = m) +
      theme_minimal()

    out <- file.path("simulations", "plots",
                     paste0("3_", tolower(gsub("[^A-Za-z0-9]", "", m)), "_vs_jump.pdf"))
    ggsave(filename = out, plot = p, width = 10, height = 5)
  }
}


################################################################################


plot_robust_scenarios <- function(n = 1000,
                                  jumpSize = 0.75,
                                  nbSeg = 8)
{
  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  patterns <- c("none", "up", "updown", "rand1")
  set.seed(999)

  # Generate clean and noisy signals with heavy-tailed noise
  signal_data <- lapply(patterns, function(pat) {
    mu <- generate_signal(n, pattern = pat, nbSeg = nbSeg, jumpSize = jumpSize)
    y <- mu + rt(length(mu), df = 2)  # t-distribution with df=2
    data.frame(t = seq_along(mu), y = y, mu = mu, pattern = pat)
  })

  df_signals <- do.call(rbind, signal_data)
  df_signals$pattern <- factor(df_signals$pattern, levels = patterns)
  rownames(df_signals) <- NULL

  p <- ggplot(df_signals) +
    geom_point(aes(x = t, y = y), alpha = 0.2) +
    geom_line(aes(x = t, y = mu), col = "red", linewidth = 0.8) +
    facet_wrap(~factor(pattern)) +
    labs(x = "Time (Index)", y = "Value") +
    theme_minimal()

  out <- file.path("simulations", "plots", "3_signal_scenarios_robust.pdf")
  ggsave(filename = out, plot = p, width = 10, height = 5)
  message("Saved plot: ", out)
  p
}


################################################################################

plot_robust_changepoint_distributions <- function(results_df,
                                                  n_bins = 50,
                                                  selected_jumpsize = 0.6)
{
  # results_df should have columns: pattern, algorithm, jump, rep, changepoints (list-column)
  # Extract all changepoints across all replicates and create a histogram

  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  # Unnest changepoints and exclude the final boundary point (end of sequence)
  cp_data_list <- list()
  for (i in seq_len(nrow(results_df)))
  {
    row <- results_df[i, ]
    cps <- unlist(row$changepoints[[1]])
    pattern <- row$pattern
    algorithm <- row$algorithm

    # Exclude the final boundary point
    cps_interior <- cps[cps < 1000]  # assuming sequence length is 1000

    if (length(cps_interior) > 0) {
      cp_df <- data.frame(
        pattern = pattern,
        algorithm = algorithm,
        changepoint = cps_interior,
        jump = row$jump
      )
      cp_data_list[[length(cp_data_list) + 1]] <- cp_df
    }
  }

  if (length(cp_data_list) > 0)
  {
    cp_all <- do.call(rbind, cp_data_list)
    rownames(cp_all) <- NULL

    cp_plot <- cp_all |>
      dplyr::filter(round(jump, 10) %in% round(selected_jumpsize, 10),
                    pattern != "none")

    p <- ggplot(cp_plot, aes(x = changepoint)) +
      geom_histogram(bins = n_bins, fill = "steelblue", alpha = 0.7, color = "black") +
      facet_grid(pattern ~ algorithm) +
      geom_vline(xintercept = seq(125, 875, by = 125), linetype = "dashed", color = "grey") +
      labs(x = "Sequence Position (Time)", y = "Frequency of Detected Change") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


    out <- file.path("simulations", "plots", "3_changepoint_distributions_robust.pdf")
    ggsave(filename = out, plot = p, width = 14, height = 10)
    message("Saved plot: ", out)
    p
  }
  else
  {
    message("No changepoints detected in any simulation.")
    NULL
  }
}


################################################################################
################################################################################
################################################################################
################################################################################

## Self-contained run & plot
DO_RUN <- TRUE
DO_PLOT <- TRUE

if (DO_RUN)
{
  results_df <- run_robust_power_study(n = 1000,
                                       patterns = c("none", "up", "updown", "rand1"),
                                       jumpSizes = seq(0.1, 4, by = 0.1),
                                       reps = 100,
                                       nbSeg = 8)
}


################################################################################

if (DO_PLOT)
{
  if (!exists("results_df")) {
    # Try to load from an RDS file that preserves the list-column
    rds_path <- file.path("simulations", "power_study_robust_heavy_tails_results.rds")
    csv_path <- file.path("simulations", "power_study_robust_heavy_tails_results.csv")
    if (file.exists(rds_path)) {
      results_df <- readRDS(rds_path)
      message("Loaded results from RDS: ", rds_path)
    } else if (file.exists(csv_path)) {
      message("Warning: CSV file found but it doesn't contain changepoints.")
      message("To generate plots with changepoint distributions, please re-run with DO_RUN=TRUE")
      results_df <- read.csv(csv_path)
    } else {
      stop("No results file found. Please run with DO_RUN=TRUE first.")
    }
  }

  plot_robust_power_metrics(results_df)
  plot_robust_scenarios(n = 1000, jumpSize = 0.6, nbSeg = 8)
  plot_robust_changepoint_distributions(results_df, n_bins = 100, selected_jumpsize = 2)
}

## End

