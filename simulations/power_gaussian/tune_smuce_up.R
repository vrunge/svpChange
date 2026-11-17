## Tune SMUCE gamma for average F1 over the Gaussian up-scenario curve.
library(svpChange)
source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_SMUCE_UP_TUNING_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_smuce_up_tuning.csv"
)

make_gaussian_up_tuning_data <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L, seed = 123L) {
  pattern_index <- match("up", POWER_PATTERNS)
  lapply(seq_along(jump_sizes), function(j) {
    lapply(seq_len(reps), function(rep) {
      design_index <- pattern_index + (j - 1L) * length(POWER_PATTERNS) +
        (rep - 1L) * length(POWER_PATTERNS) * length(jump_sizes)
      set_power_seed(seed + design_index)
      mu <- generate_power_signal(n, "up", jump_sizes[j])
      list(y = mu + stats::rnorm(n), mu = mu,
           truth = normalise_boundaries(which(diff(mu) != 0), n),
           jump = jump_sizes[j], rep = rep)
    })
  })
}

tune_gaussian_smuce_up <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    seed = 123L, workers = power_default_workers(),
    upper_gamma = read_gaussian_smuce_calibration(n = n),
    tolerance = round(n * 0.0025)) {
  data <- unlist(make_gaussian_up_tuning_data(
    n, jump_sizes, reps, seed
  ), recursive = FALSE)
  message("Prepared ", length(data), " matched up-scenario series")

  evaluate_new_gammas <- function(gammas) {
    gammas <- sort(unique(round(gammas, 6L)))
    rows <- power_lapply(seq_along(data), function(i) {
      item <- data[[i]]
      vapply(gammas, function(gamma) {
        boundaries <- normalise_boundaries(
          svp_smuce_cpp(item$y, q = gamma, sigma2 = 1), n
        )
        paper_metrics(item$truth, boundaries, tolerance, terminal = n)[["F1"]]
      }, numeric(1))
    }, workers = workers)
    values <- do.call(rbind, rows)
    list(gammas = gammas, values = values)
  }

  upper_gamma <- round(upper_gamma, 6L)
  coarse <- sort(unique(c(seq(0, upper_gamma, by = 0.05), upper_gamma)))
  first <- evaluate_new_gammas(coarse)
  best_gamma <- function(gammas, values) {
    means <- colMeans(values)
    max(gammas[means == max(means)])
  }
  coarse_best <- best_gamma(first$gammas, first$values)
  middle_grid <- seq(max(0, coarse_best - 0.05),
                     min(upper_gamma, coarse_best + 0.05), by = 0.005)
  middle <- evaluate_new_gammas(middle_grid)
  all_gammas <- c(first$gammas, middle$gammas)
  all_values <- cbind(first$values, middle$values)
  middle_best <- {
    means <- colMeans(all_values)
    max(all_gammas[means == max(means)])
  }
  fine_grid <- seq(max(0, middle_best - 0.001),
                   min(upper_gamma, middle_best + 0.001), by = 0.0001)
  fine <- evaluate_new_gammas(fine_grid)
  all_gammas <- c(all_gammas, fine$gammas)
  all_values <- cbind(all_values, fine$values)
  keep <- !duplicated(all_gammas)
  all_gammas <- all_gammas[keep]
  all_values <- all_values[, keep, drop = FALSE]
  order <- order(all_gammas)
  all_gammas <- all_gammas[order]
  all_values <- all_values[, order, drop = FALSE]

  means <- colMeans(all_values)
  gamma <- max(all_gammas[means == max(means)])
  null_raw_path <- sub("\\.csv$", ".rds", GAUSSIAN_SMUCE_CALIBRATION_FILE)
  null_raw <- readRDS(null_raw_path)
  calibration_fpr <- vapply(all_gammas, function(q) {
    mean(null_raw$calibration > q)
  }, numeric(1))
  validation_fpr <- vapply(all_gammas, function(q) {
    mean(null_raw$validation > q)
  }, numeric(1))
  selected <- all_gammas == gamma
  summary <- data.frame(
    gamma = all_gammas,
    mean_up_f1 = means,
    null_calibration_fpr = calibration_fpr,
    null_validation_fpr = validation_fpr,
    selected = selected,
    n = n,
    reps_per_jump = reps,
    jump_sizes = length(jump_sizes),
    seed = seed,
    stringsAsFactors = FALSE
  )
  rownames(summary) <- NULL
  utils::write.csv(summary, GAUSSIAN_SMUCE_UP_TUNING_FILE, row.names = FALSE)
  raw <- list(summary = summary, f1_by_series = all_values,
              jump_by_series = vapply(data, `[[`, numeric(1), "jump"),
              selected_gamma = gamma, n = n, reps = reps,
              jump_sizes = jump_sizes, seed = seed,
              upper_gamma = upper_gamma,
              search = list(coarse_step = 0.05, middle_step = 0.005,
                            fine_step = 0.0001),
              null_raw_path = null_raw_path)
  saveRDS(raw, sub("\\.csv$", ".rds", GAUSSIAN_SMUCE_UP_TUNING_FILE))
  message("Selected gamma = ", format(gamma, digits = 8),
          "; mean up-curve F1 = ", format(max(means), digits = 6),
          "; null FPR validation = ",
          format(validation_fpr[which(selected)], digits = 4))
  summary
}

if (identical(tolower(Sys.getenv("SVP_TUNE_SMUCE_UP")), "true"))
  print(tune_gaussian_smuce_up())
