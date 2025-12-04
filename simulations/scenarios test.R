library(tidyverse)

generate_signal <- function(n, pattern = c("none", "up", "updown", "rand1"), nbSeg = 8, jumpSize = 1) {
  type <- match.arg(pattern)

  if (type == "rand1") {
    set.seed(42)
    rand1CP <- rpois(nbSeg, lambda = 10)

    # Scale counts to total n
    r1 <- pmax(round(rand1CP * n / sum(rand1CP)), 1)

    # Adjust to sum exactly n
    diff <- n - sum(r1)
    r1[nbSeg] <- r1[nbSeg] + diff  # put leftover into last segment

    stopifnot(sum(r1) == n)  # safety check

    set.seed(43)
    rand1Jump <- runif(nbSeg, min = 0.5, max = 1) *
      sample(c(-1, 1), nbSeg, replace = TRUE)
  }

  # Generate scenarios
  switch(
    type,
    none = rep(0, n),
    up = rep(seq(0, nbSeg - 1) * jumpSize, each = n / nbSeg),
    updown = rep((seq(0, nbSeg - 1) %% 2) * jumpSize, each = n / nbSeg),
    rand1 = map2(rand1Jump, r1, ~rep(.x * jumpSize, .y)) %>% unlist()
  )
}

sims <- expand_grid(pattern = c("none", "up", "updown", "rand1"), rep = 1:50)

full_seqs <- pmap(sims, \(pattern, rep, jumpSize = 1) {
  mu <- generate_signal(2e3, pattern)
  set.seed(rep)
  y <- mu + rnorm(length(mu))
  cps <- c(which(diff(mu) != 0), length(mu))  # true changepoints
  return(list(y = y, mu = mu, cps = cps, pattern = pattern))
})


library(changepoint)
library(tidyverse)
library(svpChange)

# --- Helper: Mean Squared Error ---
mse_loss <- function(mu_true, mu_hat) {
  mean((mu_true - mu_hat)^2)
}

# --- Robust CP metrics (greedy 1-1 matching within tolerance) ---
cp_metrics <- function(cp_true, cp_est, tol = 5) {
  # Force plain integer vectors (handles list-columns safely)
  cp_true <- as.integer(unlist(cp_true))
  cp_est  <- as.integer(unlist(cp_est))

  # Edge cases
  if (length(cp_true) == 0L && length(cp_est) == 0L)
    return(tibble(precision = NA_real_, recall = NA_real_, F1 = NA_real_))
  if (length(cp_true) == 0L) {          # only false alarms
    precision <- 0
    recall <- NA_real_
    return(tibble(precision = precision, recall = recall, F1 = NA_real_))
  }
  if (length(cp_est) == 0L) {           # all misses
    precision <- NA_real_
    recall <- 0
    return(tibble(precision = precision, recall = recall, F1 = NA_real_))
  }

  # Pairwise absolute distances
  D <- abs(outer(cp_true, cp_est, "-"))

  # Greedy matching: just pick the closest unmatched estimate within tolerance
  matched_est <- rep(FALSE, length(cp_est))
  TP <- 0L
  for (i in seq_along(cp_true)) {
    j <- which.min(D[i, ])
    if (length(j) && D[i, j] <= tol && !matched_est[j]) {
      TP <- TP + 1L
      matched_est[j] <- TRUE
    }
  }

  FP <- sum(!matched_est)
  FN <- length(cp_true) - TP

  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  F1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_

  tibble(precision = precision, recall = recall, F1 = F1)
}

evaluate_method <- function(y, mu_true, cps_true, method, penalty) {

  if (method == "SVP") {

    penalty <- switch(penalty,
                      AIC = 2,
                      BIC = 2 * log(length(y)),
                      MBIC = 3 * log(length(y)),
                      custom = 7 * log(length(y)/8),
                      custom_2 = 8 * log(length(y)/8)
    )

    fit <- SVP(y, penalty, "gaussian_mean", prune_if_unvalid = TRUE)
    cp_est <- fit$changepoints  # integer vector of estimated CPs

    # compute the mean estimates in between changepoints
    for (i in seq_along(cp_est)) {
      if (i == 1) {
        mu_hat <- rep(mean(y[1:cp_est[i]]), cp_est[i])
      } else {
        mu_hat <- c(mu_hat, rep(mean(y[(cp_est[i - 1] + 1):cp_est[i]]), cp_est[i] - cp_est[i - 1]))
      }
    }

    # get mu_hat as vector of same length as y

  } else if (method == "BinSeg" || method == "PELT") {
    fit <- cpt.mean(y, penalty = penalty, method = method, Q = 1000)
    cp_est <- cpts(fit)  # integer vector of estimated CPs
    cp_est <- unique(c(cp_est, length(y)))
    mu_hat <- rep(param.est(fit)$mean, times = diff(c(0, cp_est)))

  } else {
    stop("Unknown method")
  }


  tibble(
    algorithm = paste0(tolower(method), "_", tolower(penalty)),
    mse = mse_loss(mu_true, mu_hat),
    cps_est = list(cp_est)   # <-- keep as list-column
  ) %>%
    bind_cols(cp_metrics(cps_true, cp_est, tol = length(mu_true) * 0.005))
}

# --- Compare across all algorithms (coerce cps at the boundary as well) ---
compare_methods <- function(y, mu, cps, pattern) {
  cps <- as.integer(unlist(cps))         # <- important when coming from transpose(list)
  algos <- expand_grid(
    penalty = c("BIC", "MBIC"),
    method  = c("PELT", "BinSeg")
  ) %>% add_row(penalty = c("custom", "custom_2"), method = "SVP")

  results <- pmap(algos, \(penalty, method)
                  evaluate_method(y, mu, cps, method, penalty)
  ) %>% list_rbind()

  results$pattern <- pattern
  results
}

# --- Run full study ---
full_out <- pmap(transpose(full_seqs), compare_methods, .progress = TRUE) |> list_rbind()

# --- Summarise ---
final_summary <- full_out %>%
  group_by(algorithm, pattern) %>%
  summarise(across(c(mse, precision, recall, F1), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pattern, values_from = c(mse, precision, recall, F1))

final_summary
