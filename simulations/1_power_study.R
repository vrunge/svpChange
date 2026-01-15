
## simulations/1_power_study.R
##
## Power study comparing SVP vs PELT on four scenarios
## - n = 10000
## - patterns = c("none","up","updown","rand1")
## - jumpSize in 1:10
## - reps per combination (default 30)
## Produces CSV and three plots (F1, precision, recall) faceted by scenario.

library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(changepoint)
library(svpChange)
library(progressr)

################################################################################

generate_signal <- function(n,
                            pattern = c("none", "up", "updown", "rand1"),
                            nbSeg = 8,
                            jumpSize = 1)
{
  type <- match.arg(pattern)

  if (type == "rand1")
  {
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

  if (length(cp_true) == 0 && length(cp_est) == 0) return(list(Precision = NA_real_, Recall = NA_real_, F1 = NA_real_))
  if (length(cp_true) == 0) return(list(Precision = 0, Recall = NA_real_, F1 = NA_real_))
  if (length(cp_est) == 0) return(list(Precision = NA_real_, Recall = 0, F1 = NA_real_))

  D <- abs(outer(cp_true, cp_est, "-"))
  matched_est <- rep(FALSE, length(cp_est))
  TP <- 0
  for (i in seq_along(cp_true))
  {
    j <- which.min(D[i, ])
    if (length(j) && D[i, j] <= tol && !matched_est[j])
    {
      TP <- TP + 1
      matched_est[j] <- TRUE
    }
  }

  FP <- sum(!matched_est)
  FN <- length(cp_true) - TP

  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
  list(Precision = precision, Recall = recall, F1 = f1)
}

################################################################################

run_power_study <- function(n = 1000,
                            patterns = c("none", "up", "updown", "rand1"),
                            jumpSizes = 1:10,
                            reps = 30,
                            workers = max(1, future::availableCores() - 1),
                            seed = 123)
{
  set.seed(seed)
  plan(multisession, workers = workers)

  combos <- expand.grid(pattern = patterns, jump = jumpSizes, stringsAsFactors = FALSE)

  # for each combo, do `reps` replicates in parallel across combos
  handlers("txtprogressbar")
  with_progress({
    # one tick per replicate (combinations * reps)
    p <- progressor(steps = nrow(combos) * reps)

    res_list <- future_lapply(seq_len(nrow(combos)), function(idx)
    {
      pattern <- combos$pattern[idx]
      jump <- combos$jump[idx]
      replicate(reps, {
        p()
        mu <- generate_signal(n, pattern = pattern, nbSeg = 8, jumpSize = jump)
        # random seed inside replicate ensures variability across workers
        y <- mu + rnorm(length(mu), mean = 0, sd = 1)
        cp_true <- c(which(diff(mu) != 0), length(mu))

        # penalties
        penalty_pelt <- 2 * log(n)
        penalty_svp  <- 1.5 * log(n)

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

        # run SVP (default SVP penalty)
        obj_svp <- tryCatch(SVP(y, penalty_svp, "gaussian_mean"), error = function(e) NULL)
        cp_svp <- integer(0)
        if (!is.null(obj_svp)) {
          if (is.list(obj_svp) && !is.null(obj_svp$changepoints)) cp_svp <- as.integer(unlist(obj_svp$changepoints))
        }
        # ensure final index included for consistency with cp_true
        cp_svp <- unique(c(cp_svp, length(y)))
        metrics_svp <- cp_metrics(cp_true, cp_svp, tol = round(n * 0.0025))
        nseg_svp <- length(cp_svp)
        mse_svp <- NA_real_
        if (length(cp_svp) > 0) {
          # build mu_hat for SVP
          mu_hat_svp <- numeric(0)
          for (i in seq_along(cp_svp)) {
            if (i == 1) {
              mu_hat_svp <- rep(mean(y[1:cp_svp[i]]), cp_svp[i])
            } else {
              mu_hat_svp <- c(mu_hat_svp, rep(mean(y[(cp_svp[i - 1] + 1):cp_svp[i]]), cp_svp[i] - cp_svp[i - 1]))
            }
          }
          if (length(mu_hat_svp) == length(mu)) mse_svp <- mean((mu - mu_hat_svp)^2)
        }

        # run SVP with BIC-like penalty (use same penalty as PELT -> call this "SVP (BIC)")
        obj_svp_bic <- tryCatch(SVP(y, penalty_pelt, "gaussian_mean"), error = function(e) NULL)
        cp_svp_bic <- integer(0)
        if (!is.null(obj_svp_bic)) {
          if (is.list(obj_svp_bic) && !is.null(obj_svp_bic$changepoints)) cp_svp_bic <- as.integer(unlist(obj_svp_bic$changepoints))
        }
        cp_svp_bic <- unique(c(cp_svp_bic, length(y)))
        metrics_svp_bic <- cp_metrics(cp_true, cp_svp_bic, tol = round(n * 0.0025))
        nseg_svp_bic <- length(cp_svp_bic)
        mse_svp_bic <- NA_real_
        if (length(cp_svp_bic) > 0) {
          mu_hat_svp_bic <- numeric(0)
          for (i in seq_along(cp_svp_bic)) {
            if (i == 1) {
              mu_hat_svp_bic <- rep(mean(y[1:cp_svp_bic[i]]), cp_svp_bic[i])
            } else {
              mu_hat_svp_bic <- c(mu_hat_svp_bic, rep(mean(y[(cp_svp_bic[i - 1] + 1):cp_svp_bic[i]]), cp_svp_bic[i] - cp_svp_bic[i - 1]))
            }
          }
          if (length(mu_hat_svp_bic) == length(mu)) mse_svp_bic <- mean((mu - mu_hat_svp_bic)^2)
        }

        tib <- tibble::tibble(
          pattern = pattern,
          jump = jump,
          rep = as.integer(NA),
          algorithm = c("PELT", "SVP", "SVP (BIC)"),
          Precision = c(metrics_pelt$Precision, metrics_svp$Precision, metrics_svp_bic$Precision),
          Recall = c(metrics_pelt$Recall, metrics_svp$Recall, metrics_svp_bic$Recall),
          F1 = c(metrics_pelt$F1, metrics_svp$F1, metrics_svp_bic$F1),
          NumSegments = c(nseg_pelt, nseg_svp, nseg_svp_bic),
          MSE = c(mse_pelt, mse_svp, mse_svp_bic),
          changepoints = list(cp_pelt, cp_svp, cp_svp_bic)
        )
        tib
      }, simplify = FALSE)
    }, future.seed = TRUE)
  })

  plan(sequential)

  # flatten: produce one data.frame from the nested list of tibbles
  all_reps <- unlist(res_list, recursive = FALSE)
  df <- dplyr::bind_rows(all_reps)

  # Add rep id
  df <- df %>% group_by(pattern, jump, algorithm) %>% mutate(rep = row_number()) %>% ungroup()

  # Save a version without the list-column for inspection (drop changepoints for CSV export)
  df_csv <- df %>% select(-changepoints)
  out_path <- file.path("simulations", "power_study_results.csv")
  write.csv(df_csv, out_path, row.names = FALSE)
  message("Saved results to: ", out_path)

  # Return the full df with changepoints intact for use in plotting
  df
}


################################################################################


plot_power_metrics <- function(df)
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

  long <- df %>% pivot_longer(cols = c(Precision, Recall, F1, NumSegments, MSE), names_to = "metric", values_to = "value")
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

  metrics <- c("F1", "Precision", "Recall", "NumSegments", "MSE")
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
                     paste0("1_", tolower(gsub("[^A-Za-z0-9]", "", m)), "_vs_jump.pdf"))
    ggsave(filename = out, plot = p, width = 10, height = 5)
  }
}


################################################################################

plot_scenarios <- function(n = 1000,
                           jumpSize = 0.75,
                           nbSeg = 8)
{
  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  patterns <- c("none", "up", "updown", "rand1")
  set.seed(999)

  # Generate clean and noisy signals
  signal_data <- lapply(patterns, function(pat) {
    mu <- generate_signal(n, pattern = pat, nbSeg = nbSeg, jumpSize = jumpSize)
    y <- mu + rnorm(length(mu), mean = 0, sd = 1)
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

  out <- file.path("simulations", "plots", "1_signal_scenarios.pdf")
  ggsave(filename = out, plot = p, width = 10, height = 5)
  message("Saved plot: ", out)
  p
}

################################################################################

plot_changepoint_distributions <- function(df,
                                           n_bins = 50,
                                           selected_jumpsize = 2)
{
  # df should have columns: pattern, algorithm, jump, rep, changepoints (list-column)
  # Extract all changepoints across all replicates and create a histogram

  dir.create(file.path("simulations", "plots"), recursive = TRUE, showWarnings = FALSE)

  # Unnest changepoints and exclude the final boundary point (end of sequence)
  cp_data_list <- list()
  for (i in seq_len(nrow(df)))
  {
    row <- df[i, ]
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


    # Create histogram with ggplot2
    cp_plot <- cp_all |>
      dplyr::filter(round(jump, 10) %in% round(selected_jumpsize, 10),
                    pattern != "none")

    p <- ggplot(cp_plot, aes(x = changepoint)) +
      geom_histogram(bins = n_bins, fill = "steelblue", alpha = 0.7, color = "black") +
      facet_grid(pattern ~ algorithm) +
      geom_vline(xintercept = seq(125, 875, by = 125), linetype = "dashed", color = "grey") +
      labs(x = "Sequence Position (Time)",
           y = "Frequency of Detected Change") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    out <- file.path("simulations", "plots", "1_changepoint_distributions.pdf")
    ggsave(filename = out, plot = p, width = 14, height = 10)
    message("Saved plot: ", out)
    p
  } else {
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
  df <- run_power_study(n = 1000, patterns = c("none", "up", "updown", "rand1"), jumpSizes = c(seq(0.1, 2, by = 0.1)), reps = 100)
}

################################################################################


if (DO_PLOT)
{
  if (!exists("df"))
  {
    # Try to load from an RDS file that preserves the list-column
    rds_path <- file.path("simulations", "power_study_results.rds")
    csv_path <- file.path("simulations", "power_study_results.csv")
    if (file.exists(rds_path))
    {
      df <- readRDS(rds_path)
      message("Loaded results from RDS: ", rds_path)
    }
    else if (file.exists(csv_path))
    {
      message("Warning: CSV file found but it doesn't contain changepoints.")
      message("To generate plots with changepoint distributions, please re-run with DO_RUN=TRUE")
      df <- read.csv(csv_path)
    }
    else
    {
      stop("No results file found. Please run with DO_RUN=TRUE first.")
    }
  }
  else
  {
    # Save the full df with changepoints to RDS for future reloading
    rds_path <- file.path("simulations", "power_study_results.rds")
    saveRDS(df, rds_path)
    message("Saved full results (with changepoints) to: ", rds_path)
  }

  plot_power_metrics(df)
  plot_scenarios(n = 1000, jumpSize = 0.6, nbSeg = 8)

  plot_changepoint_distributions(df, n_bins = 100, selected_jumpsize = 0.6)

}

## End
