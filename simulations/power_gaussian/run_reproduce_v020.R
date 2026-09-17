## Reproduce Section 5.1 with the current svpChange source without
## overwriting the historical power_gaussian/results.* files.

source_path <- file.path("simulations", "power_gaussian", "run_power.R")
old_flag <- Sys.getenv("SVP_RUN_SIMULATIONS", unset = NA_character_)
Sys.setenv(SVP_RUN_SIMULATIONS = "false")
source(source_path)
if (is.na(old_flag)) Sys.unsetenv("SVP_RUN_SIMULATIONS") else
  Sys.setenv(SVP_RUN_SIMULATIONS = old_flag)

matched_root <- file.path(
  "simulations", "power_gaussian", "reproduced_v020_matched"
)

run_and_save_v020_matched <- function(workers = power_default_workers()) {
  run_and_save_gaussian(workers = workers, root = matched_root)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_v020_matched()
}
