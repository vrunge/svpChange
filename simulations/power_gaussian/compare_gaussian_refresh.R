## Separate the effect of the corrected changepoint scoring from the effect of
## recalibrating the Gaussian SVP constants.

library(dplyr)

source(file.path("simulations", "power_common.R"))

GAUSSIAN_COMPARISON_ROOT <- file.path("simulations", "power_gaussian")

gaussian_configuration <- function(algorithm, refreshed = FALSE) {
  ifelse(
    algorithm == "PELT", "PELT",
    ifelse(
      grepl("right / c = 2\\.000|^SVP BIC$", algorithm), "right c=2",
      ifelse(
        grepl("right / c = 1\\.642|^SVP BIC calibrated$", algorithm),
        "right calibrated",
        ifelse(
          grepl("both / c = 1\\.987|^SVP BIC multiscale$", algorithm),
          "both calibrated", algorithm
        )
      )
    )
  )
}

rescore_gaussian_results <- function(results, n = 1000L,
                                     tolerance = round(0.0025 * n)) {
  metrics <- lapply(seq_len(nrow(results)), function(i) {
    mu <- generate_power_signal(n, results$pattern[i], results$jump[i])
    truth <- normalise_boundaries(which(diff(mu) != 0), n)
    paper_metrics(
      truth, results$changepoints[[i]], tolerance = tolerance, terminal = n
    )
  })
  metrics <- as.data.frame(do.call(rbind, metrics))
  for (name in names(metrics)) results[[name]] <- metrics[[name]]
  results
}

summarise_gaussian_stage <- function(results, stage) {
  results |>
    mutate(configuration = gaussian_configuration(algorithm)) |>
    group_by(configuration, pattern) |>
    summarise(
      stage = stage,
      mean_F1 = mean(F1),
      mean_precision = mean(Precision),
      mean_recall = mean(Recall),
      correct_number_probability = mean(CorrectNumCP),
      mean_segments = mean(NumSegments),
      mean_MSE = mean(MSE),
      .groups = "drop"
    ) |>
    select(stage, everything())
}

compare_gaussian_refresh <- function(root = GAUSSIAN_COMPARISON_ROOT) {
  historical <- readRDS(file.path(
    root, "reproduced_v020_matched", "results.rds"
  ))
  refreshed <- readRDS(file.path(root, "results.rds"))
  historical_rescored <- rescore_gaussian_results(historical)

  decomposition <- bind_rows(
    summarise_gaussian_stage(historical, "historical stored score"),
    summarise_gaussian_stage(historical_rescored, "historical fits, corrected score"),
    summarise_gaussian_stage(refreshed, "refreshed calibration and corrected score")
  )

  before <- summarise_gaussian_stage(
    historical_rescored, "historical fits, corrected score"
  ) |>
    rename_with(~paste0("before_", .x), -c(configuration, pattern, stage)) |>
    select(-stage)
  after <- summarise_gaussian_stage(
    refreshed, "refreshed calibration and corrected score"
  ) |>
    rename_with(~paste0("after_", .x), -c(configuration, pattern, stage)) |>
    select(-stage)
  changes <- full_join(before, after, by = c("configuration", "pattern")) |>
    mutate(
      old_constant = case_when(
        configuration == "right calibrated" ~ 1.5,
        configuration == "both calibrated" ~ 1.8,
        configuration == "right c=2" ~ 2,
        TRUE ~ NA_real_
      ),
      new_constant = case_when(
        configuration == "right calibrated" ~ 1.642,
        configuration == "both calibrated" ~ 1.987,
        configuration == "right c=2" ~ 2,
        TRUE ~ NA_real_
      ),
      delta_F1 = after_mean_F1 - before_mean_F1,
      delta_correct_number = after_correct_number_probability -
        before_correct_number_probability,
      segment_ratio = after_mean_segments / before_mean_segments,
      flag_material = abs(delta_F1) >= 0.03 |
        abs(delta_correct_number) >= 0.03 |
        abs(segment_ratio - 1) >= 0.10
    )

  write.csv(
    decomposition, file.path(root, "gaussian_performance_decomposition.csv"),
    row.names = FALSE
  )
  write.csv(
    changes, file.path(root, "gaussian_performance_changes.csv"),
    row.names = FALSE
  )
  saveRDS(
    historical_rescored,
    file.path(root, "reproduced_v020_matched", "results_corrected_score.rds")
  )
  invisible(list(decomposition = decomposition, changes = changes))
}

if (identical(tolower(Sys.getenv("SVP_COMPARE_GAUSSIAN")), "true")) {
  compare_gaussian_refresh()
}
