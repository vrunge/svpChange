##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' FOCuS Validity Test for Every Prefix
#'
#' @title FOCuS Validity Test for Every Prefix
#' @description Tests whether the FOCuS statistic remains below `gamma` for
#' every prefix of the segment. It is the prefix-wise validity rule available
#' for use with [svp0()].
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the FOCuS statistic.
#' @return `TRUE` if the statistic is strictly below `gamma` for every prefix
#' of `y`; otherwise, `FALSE`.
#' @details `valid_FOCUS()` updates the FOCuS statistic after each observation
#' and returns `FALSE` as soon as a prefix statistic reaches or exceeds
#' `gamma`. Use it with [svp0()] when every prefix must be valid. In contrast,
#' [SVP()] with `test = "gaussian_mean"` checks only the statistic at the
#' current segment endpoint, and [valid_FOCUS_last()] checks only the final
#' statistic of a supplied segment.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 40, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_FOCUS(data, gamma = 2 * log(length(data)))
#' @seealso [valid_FOCUS_last()]
#' @export
valid_FOCUS <- function(y, gamma) {
  .focus_valid_cpp(y, gamma, check_all_prefixes = TRUE)
}

################################################
#' FOCuS Validity Test for the Final Segment
#'
#' @title FOCuS Validity Test for the Final Segment
#' @description Tests whether the FOCuS statistic for the complete segment is
#' below `gamma`, without checking intermediate prefixes.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the FOCuS statistic.
#' @return `TRUE` if the final statistic is strictly below `gamma`; otherwise,
#' `FALSE`.
#' @details `valid_FOCUS_last()` processes the complete segment and checks only
#' the final statistic. An earlier prefix may have reached or exceeded `gamma`
#' without making the final segment invalid. Use [valid_FOCUS()] to require
#' every prefix to remain valid.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 40, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_FOCUS_last(data, gamma = 2 * log(length(data)))
#' @seealso [valid_FOCUS()]
#' @export
valid_FOCUS_last <- function(y, gamma) {
  .focus_valid_cpp(y, gamma, check_all_prefixes = FALSE)
}

################################################

# Estimate the AR(1) coefficient using the robust estimator used by the C++
# AR(1) implementation.
.valid_AR1_rho <- function(y) {
  if (length(y) < 3L) {
    stop("AR(1) estimation requires at least 3 observations")
  }
  lag1 <- abs(diff(y))
  lag2 <- abs(y[seq.int(3L, length(y))] - y[seq_len(length(y) - 2L)])
  med1 <- stats::median(lag1)
  med2 <- stats::median(lag2)
  if (!is.finite(med1) || !is.finite(med2) || med1 <= 0) {
    stop("Cannot estimate AR(1) correlation from constant or non-finite data")
  }
  max(-0.999, min(0.999, med2^2 / med1^2 - 1))
}

#' AR(1) Mean-Change Validity Test
#'
#' @title AR(1) Mean-Change Validity Test
#' @description Tests whether a segment is valid under a conditional Gaussian
#' AR(1) single-mean-change statistic implemented in R.
#' @param y A numeric AR(1) segment. At least four observations are needed for
#' a non-trivial single-change scan.
#' @param gamma A threshold for the AR(1) likelihood-ratio statistic.
#' @param rho An optional known AR(1) coefficient. If `NA`, it is estimated
#' robustly from `y`.
#' @param sigma2 A positive innovation variance used by the known-variance
#' statistic. In the AR(1) model, this is the variance of the new random shock
#' after accounting for the previous observation and the AR(1) mean structure;
#' it is not the marginal variance of the observed series. It is ignored when
#' `profile_sigma = TRUE`.
#' @param profile_sigma If `TRUE`, estimate the innovation variance separately
#' under the no-change and one-change models from their residual sums of
#' squares. This is useful when the innovation scale is unknown.
#' @return `TRUE` if the AR(1) statistic is strictly below `gamma`; otherwise,
#' `FALSE`.
#' @details The function conditions on the first observation, transforms the
#' remaining observations into AR(1) innovations, and scans change locations
#' that leave at least two observations on each side.
#' The statistic is one half of the Gaussian likelihood-ratio statistic, using
#' the same convention as [AR1_single_change()]. If `rho` is `NA`, the robust
#' estimator requires at least three observations. The function can be passed
#' to [svp0()] directly when `rho` is fixed, or through a closure when other
#' optional arguments are needed. The caller should enforce any minimum segment
#' length required by the chosen SVP procedure. When `profile_sigma = FALSE`,
#' the statistic uses the supplied `sigma2`; when `profile_sigma = TRUE`, the
#' innovation variance is profiled out and estimated under each competing
#' model.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 60, parameters = 0, sd_noise = 1,
#'   rho = 0.6, type = "gaussAR1"
#' )
#' valid_AR1(data, gamma = 10, rho = 0.6, sigma2 = 1)
#' @seealso [AR1_single_change()], [svp0()]
#' @export
valid_AR1 <- function(y, gamma, rho = NA_real_, sigma2 = 1,
                      profile_sigma = FALSE) {
  stopifnot(
    is.numeric(y),
    all(is.finite(y)),
    length(rho) == 1L,
    is.numeric(rho),
    length(sigma2) == 1L,
    is.numeric(sigma2),
    length(profile_sigma) == 1L,
    is.logical(profile_sigma),
    !is.na(profile_sigma)
  )

  if (!profile_sigma) {
    stopifnot(is.finite(sigma2), sigma2 > 0)
  }

  if (is.na(rho)) {
    rho <- .valid_AR1_rho(y)
  }
  stopifnot(is.finite(rho), abs(rho) < 1)

  n <- length(y)
  if (n < 4L) {
    return(0 < gamma)
  }

  innovations <- y[-1L] - rho * y[-n]
  prefix_sum <- c(0, cumsum(innovations))
  prefix_sq <- c(0, cumsum(innovations^2))
  transition_count <- n - 1L
  total_sum <- prefix_sum[transition_count + 1L]
  total_sq <- prefix_sq[transition_count + 1L]
  rss0 <- max(0, total_sq - total_sum^2 / transition_count)

  tau <- 2L:(n - 2L)
  pre_count <- tau - 1L
  post_count <- n - tau - 1L
  pre_sum <- prefix_sum[tau]
  transition <- innovations[tau]
  post_sum <- total_sum - prefix_sum[tau + 1L]
  one_minus_rho <- 1 - rho
  a <- pre_count * one_minus_rho^2 + rho^2
  b <- -rho
  d <- 1 + post_count * one_minus_rho^2
  v1 <- one_minus_rho * pre_sum - rho * transition
  v2 <- transition + one_minus_rho * post_sum
  determinant <- a * d - b^2
  eligible <- determinant > .Machine$double.eps
  if (!any(eligible)) {
    return(0 < gamma)
  }
  fitted <- (d[eligible] * v1[eligible]^2 -
    2 * b * v1[eligible] * v2[eligible] +
    a[eligible] * v2[eligible]^2) / determinant[eligible]
  rss1 <- min(pmax(0, total_sq - fitted))

  statistic <- if (profile_sigma) {
    if (rss0 > 0 && rss1 <= .Machine$double.eps * rss0) {
      Inf
    } else if (rss1 > 0 && rss0 > rss1) {
      0.5 * transition_count * log(rss0 / rss1)
    } else {
      0
    }
  } else {
    max(0, (rss0 - rss1) / (2 * sigma2))
  }
  statistic < gamma
}

################################################

#' Validity Test Based on Sum of Squared Errors
#'
#' @title Validity Test Based on Sum of Squared Errors
#' @description Checks if the sum of squared deviations from the segment mean
#' is lower than or equal to a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed sum of squared
#' errors.
#' @return TRUE if the sum of squared errors is less than or equal to gamma.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 30, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_SSE(data, gamma = 40)
#' @export
valid_SSE <- function(y, gamma) {
  sum((y - mean(y))^2) <= gamma
}

################################################

#' Validity Test Based on Range
#'
#' @title Validity Test Based on Range
#' @description Checks whether the range (max - min) of the segment is below a
#' threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed range.
#' @return TRUE if the range is less than or equal to gamma.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 30, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_RANGE(data, gamma = 5)
#' @export
valid_RANGE <- function(y, gamma) {
  max(y) - min(y) <= gamma
}

################################################

#' Validity Test Based on Trimmed Range
#'
#' @title Validity Test Based on Trimmed Range
#' @description Applies a slack range test using trimmed minimum and maximum
#' (ignores `trim` observations at each end).
#' @param y A numeric vector representing a segment with at least
#' `2 * trim + 1` observations.
#' @param gamma A numeric threshold for the trimmed range.
#' @param trim Number of observations to remove from each end. Defaults to 3
#' and must be a non-negative integer.
#' @return TRUE if the trimmed range is less than or equal to gamma.
#' @details The test discards `trim` smallest and `trim` largest observations
#' before computing the range. At least `2 * trim + 1`
#' observations are therefore required, leaving one observation after
#' trimming. The minimum segment length is enforced by the caller.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 30, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_RANGE_SLACK(data, gamma = 4, trim = 3)
#' @export
valid_RANGE_SLACK <- function(y, gamma, trim = 3) {
  stopifnot(
    length(trim) == 1L,
    is.numeric(trim),
    is.finite(trim),
    trim >= 0,
    trim == as.integer(trim)
  )
  trim <- as.integer(trim)
  sorted_y <- sort(y)
  sorted_y[length(y) - trim] - sorted_y[trim + 1L] <= gamma
}

################################################

#' Validity Test Based on Interquantile Range
#'
#' @title Validity Test Based on Interquantile Range
#' @description Checks whether the quantile range between `probs[1]` and
#' `probs[2]` is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the interquantile range.
#' @param probs A numeric vector of length two giving the lower and upper
#' quantile probabilities. Defaults to `c(0.05, 0.95)` and must be in
#' increasing order.
#' @return TRUE if the interquantile range is less than or equal to gamma.
#' @details The test computes the difference between the quantiles at
#' `probs[2]` and `probs[1]`.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 30, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_QUANTILE(data, gamma = 4)
#' @export
valid_QUANTILE <- function(y, gamma, probs = c(0.05, 0.95)) {
  stopifnot(
    length(probs) == 2L,
    is.numeric(probs),
    all(is.finite(probs)),
    all(probs >= 0, probs <= 1),
    probs[1] < probs[2]
  )
  unname(diff(quantile(y, probs = probs))) <= gamma
}

################################################

#' Validity Test Based on Robust Scale
#'
#' @title Validity Test Based on Robust Scale
#' @description Checks whether the Gaussian-normalized median absolute
#' deviation of a segment is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the robust scale of the segment.
#' @return `TRUE` if the normalized median absolute deviation is less than or
#' equal to `gamma`; otherwise, `FALSE`.
#' @details The scale is computed with [stats::mad()] using its segment median
#' as centre and the constant `1.4826`. This constant makes the MAD estimate
#' the standard deviation for Gaussian data. Because the test is based on the
#' complete segment, its result is not necessarily preserved when observations
#' are appended; use `subtests = "none"` and `PELT_pruning = FALSE` in [svp0()]
#' unless the required pruning properties have been established for the
#' intended application.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 30, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' valid_SCALE(data, gamma = 2)
#' @export
valid_SCALE <- function(y, gamma) {
  stats::mad(y, constant = 1.4826) <= gamma
}

################################################

#' Optimal Partitioning Cost Test
#'
#' @title Optimal Partitioning Cost Test
#' @description Tests whether the total cost of the segment is smaller than
#' the best two-part penalized segmentation cost.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for penalizing the introduction of a new
#' segment.
#' @return TRUE if the segment cannot be split into two parts with lower
#' (penalized) cost.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = c(20, 40), parameters = c(0, 2),
#'   sd_noise = 1, type = "gauss"
#' )
#' valid_OP(data, gamma = 2 * log(length(data)))
#' @export
valid_OP <- function(y, gamma) {
  len <- length(y)
  if (len == 1) {
    return(TRUE)
  }
  total <- sum((y - mean(y))^2)
  split <- seq_len(len - 1L)
  prefix_sum <- c(0, cumsum(y))
  prefix_sq_sum <- c(0, cumsum(y^2))
  left_sum <- prefix_sum[split + 1L]
  right_sum <- prefix_sum[len + 1L] - left_sum
  left_sse <- prefix_sq_sum[split + 1L] - left_sum^2 / split
  right_sse <- (prefix_sq_sum[len + 1L] - prefix_sq_sum[split + 1L]) -
    right_sum^2 / (len - split)
  val <- min(left_sse + right_sse)
  total < (val + 2 * gamma)
}
