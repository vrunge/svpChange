## Gaussian power study: PELT (calibrated), SVP none (calibrated), SVP right
## with a fixed, uncalibrated penalty of 2 * log(n), SVP right (calibrated),
## and SVP both (calibrated), all at a 1% null false-detection rate except
## the uncalibrated curve.
library(svpChange)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))
GAUSSIAN_ROOT <- file.path("simulations", "power_gaussian")

fit_gaussian_methods <- function(
    y, constants = read_gaussian_calibration(),
    pelt_constant = read_gaussian_pelt_calibration()) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = pelt_constant * log(n))
  methods <- list(
    normalise_boundaries(changepoint::cpts(pelt), n),
    normalise_boundaries(SVP(y, constants[["none"]] * log(n),
                            subtests = "none")$changepoints, n),
    normalise_boundaries(SVP(y, GAUSSIAN_UNCALIBRATED_RIGHT_CONSTANT * log(n),
                            subtests = "right")$changepoints, n),
    normalise_boundaries(SVP(y, constants[["right"]] * log(n),
                            subtests = "right")$changepoints, n),
    normalise_boundaries(SVP(y, constants[["both"]] * log(n),
                            subtests = "both")$changepoints, n))
  names(methods) <- gaussian_power_method_labels(constants, pelt_constant)
  methods
}

run_gaussian_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 4523167L,
    constants = read_gaussian_calibration(),
    pelt_constant = read_gaussian_pelt_calibration()) {
  run_power_grid(
    n, jump_sizes, reps, simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) {
      fit_gaussian_methods(y, constants, pelt_constant)
    },
    tolerance = round(n * 0.0025), workers = workers, seed = seed)
}

run_and_save_gaussian_all <- function(
    workers = power_default_workers(), root = GAUSSIAN_ROOT) {
  results <- run_gaussian_power(workers = workers)
  order <- validate_gaussian_power_results(results)
  plots <- save_power_outputs(
    results, root, power_scenario_plot(1000L, 0.6, stats::rnorm),
    plot_results = set_algorithm_order(results, order))
  invisible(results)
}

## Seven-method highlight (PELT, SVP both, SVP SMUCE, and four SVP variants
## explored in multiscale_methods_tested.pdf), kept independent of the main
## none/right/both comparison above. PELT/SVP both/SVP SMUCE are always
## first, since export_paper_figures.R's f1_three_methods.pdf subsets to
## exactly these three.
##
## "SVP both_correction" (subtests = "both_correction"): the multiscale simultaneous-
## CUSUM penalty from SVP_consistency_sharp.tex, q calibrated by simulation.
##
## "SVP both_calibrated Bonferroni" uses fitted constants (A=1.683,
## B=6.209, C=4.939 originally) unchanged for A and B, with C re-derived
## (calibrate_bonferroni_both_calibrated2.R) so that
## gamma(l)=A*log(l)+B*log(log(l))+C, applied once per distinct candidate
## length l within a "both"-pruned run, gives a true 1% null FPR: C=15.5,
## found by bisecting the real discrete FPR directly (two alternative
## correction shapes were also tried and rejected for lower power: a global
## multiplicative scale w=1.7109 on the whole fitted curve, and re-deriving
## A directly to A=6.2302 -- see multiscale_methods_tested.pdf for all
## three).
BOTH_CALIBRATED_BONFERRONI <- c(A = 1.683, B = 6.209, C = 15.5)

## "SVP both_rmin" (subtests = "both_rmin"): the minimum-arm restricted
## CUSUM, excluding splits whose shorter arm is below r_min. Calibrated by
## grid search over r_min (calibrate_rmin.R): r_min = 30, gamma bisected to
## a true 1% null FPR.
BOTH_RMIN <- c(r_min = 30, gamma = 25.3125)

## "SVP both_empirical" uses subtests = "both_empirical": a "both"-pruned
## test whose threshold gamma(t) = E_t + c*sqrt(V_t) varies with the CURRENT
## TOTAL observed length t (shared by every simultaneously alive candidate
## at that moment), with E_t, V_t, and c all estimated directly by
## simulation (calibrate_empirical_both.R) rather than any parametric or
## asymptotic formula.
##
## "SVP both_lengthonly" (subtests = "both_lengthonly"): the
## SVP_validity_rules.tex proposal, f_Gauss(y_{s:t}) <= [q +
## sqrt(2*log(e*n/m))]^2 with m the candidate's own current length and n
## the total series length -- a Dumbgen-Spokoiny/SMUCE-style multiscale
## correction, dispatched with "both" pruning like every other both-family
## method here. q calibrated by simulation (calibrate_lengthonly.R):
## q = 3.30, true 1% null FPR (0.90% at 8000 replicates).
LENGTHONLY_Q <- 3.30

## "SVP both_mosum" (subtests = "both_mosum"): the general MOSUM
## (moving-sum) test replaces the FOCUS/CUSUM statistic used by every
## other method with a comparison of the mean of G points to the right of
## a position k against G points to the left, maximized over every
## admissible k in [G, m-G], m = candidate's own current length. At its
## MAXIMAL bandwidth G = floor(m/2), that range collapses to a single
## point k = floor(m/2) (the segment's own midpoint), so this reduces to
## the single-anchored Gaussian LR statistic at a SINGLE, fixed split --
## the candidate's own midpoint -- instead of maximized over every
## internal split like f_Gauss. Segments shorter than 20 are always valid
## (hard floor, not asymptotic). gamma is a flat threshold calibrated by
## simulation (calibrate_mosum.R): gamma = 23.32, true 1% null FPR (0.96%
## at 8000 replicates).
##
## "SVP both_mosum_sqrt" (subtests = "both_mosum_sqrt"): the same general
## MOSUM test at a smaller, length-growing bandwidth G(m) = max(2,
## floor(sqrt(m))) instead of the maximal G = floor(m/2) used by
## "both_mosum". The threshold is the exact asymptotic boundary of
## Eichinger & Kirch (2018, Bernoulli 24(1), Theorem 2.1), evaluated at m
## in place of the total series length, plus a single additive correction
## gamma calibrated by simulation (calibrate_mosum_sqrt.R): gamma = 1.625,
## true 1% null FPR (0.76% at 8000 replicates). Substantially less
## powerful than the FOCUS-based methods at small-to-moderate jump sizes
## (a real property of this narrower-bandwidth MOSUM scan, not a
## calibration artifact -- see multiscale_methods_tested.pdf).
##
## "SVP both_mosum_adaptive" (subtests = "both_mosum_adaptive"): the
## length-adaptive version of both_mosum. The raw statistic never
## changes (still a single fixed-midpoint comparison, no internal
## maximum), but the threshold now grows as
## gamma(m) = [q + sqrt(2*log(log(m/20) + e))]^2, a double-log
## (iterated-logarithm) correction motivated by the fact that the raw
## statistic is already stationary across lengths (unlike f_Gauss), so a
## much gentler correction than the single-log f_Gauss-based methods need
## suffices -- see src/svp_tests_mosum_adaptive.h and
## multiscale_methods_tested.pdf for the full derivation. q calibrated by
## simulation (calibrate_mosum_adaptive.R): q = 3.304, true 1% null FPR
## (0.84% at 8000 replicates). See multiscale_methods_tested.pdf for all
## eleven methods.
gaussian_four_method_labels <- function(
    constants = read_gaussian_calibration(),
    pelt_constant = read_gaussian_pelt_calibration(),
    smuce_gamma = read_gaussian_smuce_calibration()) {
  c(paste0("PELT / c = ", formatC(pelt_constant, digits = 3L, format = "f")),
    paste0("SVP both / c = ", formatC(constants[["both"]], digits = 3L, format = "f")),
    paste0("SVP SMUCE / gamma = ", formatC(smuce_gamma, digits = 3L, format = "f")),
    "SVP both_correction (1% FPR)",
    "SVP both_calibrated Bonferroni (1% FPR)",
    "SVP both_rmin (1% FPR)",
    "SVP both_empirical (1% FPR)",
    "SVP both_lengthonly (1% FPR)",
    "SVP both_mosum (1% FPR)",
    "SVP both_mosum_sqrt (1% FPR)",
    "SVP both_mosum_adaptive (1% FPR)")
}

fit_gaussian_four_methods <- function(
    y, constants = read_gaussian_calibration(),
    pelt_constant = read_gaussian_pelt_calibration(),
    smuce_gamma = read_gaussian_smuce_calibration(n = length(y)),
    both_correction_q = read_gaussian_consistency_sharp_calibration(n = length(y))) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = pelt_constant * log(n))
  methods <- list(
    normalise_boundaries(changepoint::cpts(pelt), n),
    normalise_boundaries(SVP(y, constants[["both"]] * log(n),
                            subtests = "both")$changepoints, n),
    normalise_boundaries(svp_smuce_cpp(y, q = smuce_gamma, sigma2 = 1), n),
    normalise_boundaries(SVP(y, both_correction_q, subtests = "both_correction")$changepoints, n),
    normalise_boundaries(SVP(y, 0, subtests = "both_calibrated",
                            abc_A = BOTH_CALIBRATED_BONFERRONI[["A"]],
                            abc_B = BOTH_CALIBRATED_BONFERRONI[["B"]],
                            abc_C = BOTH_CALIBRATED_BONFERRONI[["C"]])$changepoints, n),
    normalise_boundaries(SVP(y, BOTH_RMIN[["gamma"]], subtests = "both_rmin",
                            r_min = BOTH_RMIN[["r_min"]])$changepoints, n),
    normalise_boundaries(SVP(y, 0, subtests = "both_empirical")$changepoints, n),
    normalise_boundaries(SVP(y, LENGTHONLY_Q, subtests = "both_lengthonly")$changepoints, n),
    normalise_boundaries(SVP(y, 23.32, subtests = "both_mosum")$changepoints, n),
    normalise_boundaries(SVP(y, 1.625, subtests = "both_mosum_sqrt")$changepoints, n),
    normalise_boundaries(SVP(y, 3.304, subtests = "both_mosum_adaptive")$changepoints, n))
  names(methods) <- gaussian_four_method_labels(
    constants, pelt_constant, smuce_gamma)
  methods
}

run_and_save_gaussian_four_methods <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 7791423L, root = GAUSSIAN_ROOT) {
  labels <- gaussian_four_method_labels()
  results <- run_power_grid(
    n, jump_sizes, reps, simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) fit_gaussian_four_methods(y),
    tolerance = round(n * 0.0025), workers = workers, seed = seed)
  results <- set_algorithm_order(results, labels)
  saveRDS(results, file.path(root, "f1_four_methods_results.rds"))
  utils::write.csv(dplyr::select(results, -changepoints),
                   file.path(root, "f1_four_methods_results.csv"),
                   row.names = FALSE)
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true"))
  run_and_save_gaussian_all()
