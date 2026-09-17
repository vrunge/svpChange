# Validate the common fixed-variance Gaussian SMUCE arguments.
.smuce_validate <- function(y, gamma, sigma2, n) {
  stopifnot(
    is.numeric(y),
    length(y) > 0L,
    all(is.finite(y)),
    length(gamma) == 1L,
    is.numeric(gamma),
    is.finite(gamma),
    gamma >= 0,
    length(sigma2) == 1L,
    is.numeric(sigma2),
    is.finite(sigma2),
    sigma2 > 0,
    length(n) == 1L,
    is.numeric(n),
    is.finite(n),
    n >= length(y),
    n == floor(n)
  )
}

#' Internal helper: admissible interval of constant means
#' @param y Numeric observations.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @param n Total series length.
#' @return Numeric lower and upper admissible bounds.
#' @details
#' For the independent Gaussian model with known variance, the returned
#' interval is the set of all constant means satisfying the SMUCE constraint
#' on every subinterval of `y`. An interval containing `NA` values means that
#' no constant mean is admissible.
#' @references
#' Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point
#' Inference. *Journal of the Royal Statistical Society: Series B*, 76(3),
#' 495--580. doi:10.1111/rssb.12047.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 20, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' smuce_theta_interval(data, gamma = 2, sigma2 = 1)
#' @keywords internal
smuce_theta_interval <- function(y, gamma, sigma2 = 1, n = length(y)) {
  .smuce_validate(y, gamma, sigma2, n)
  # Center before prefix sums to avoid cancellation for translated data.
  offset <- y[1L]
  y_centered <- y - offset
  lo <- -Inf
  hi <- Inf
  cs <- c(0, cumsum(y_centered))
  for (u in seq_along(y_centered)) {
    for (v in u:length(y_centered)) {
      m <- v - u + 1
      centre <- (cs[v + 1] - cs[u]) / m
      r <- sqrt(sigma2 / m) *
        (gamma + sqrt(2 * (1 + log(n / m))))
      lo <- max(lo, centre - r)
      hi <- min(hi, centre + r)
      if (lo > hi) {
        return(c(NA_real_, NA_real_))
      }
    }
  }
  c(lo + offset, hi + offset)
}

#' Constrained Gaussian cost for a SMUCE-valid segment
#' @param y Numeric observations in the candidate segment.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @param n Total series length.
#' @return Numeric constrained Gaussian residual cost.
#' @details
#' The cost is the Gaussian residual sum of squares divided by `sigma2`,
#' minimized over the SMUCE-admissible constant means. It is therefore the
#' constrained Gaussian negative log-likelihood cost up to an additive and
#' multiplicative constant. `Inf` is returned when the segment has no
#' admissible constant mean.
#' @references
#' Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point
#' Inference. *Journal of the Royal Statistical Society: Series B*, 76(3),
#' 495--580. doi:10.1111/rssb.12047.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = 20, parameters = 0, sd_noise = 1, type = "gauss"
#' )
#' smuce_cost(data, gamma = 2, sigma2 = 1)
#' @export
smuce_cost <- function(y, gamma, sigma2 = 1, n = length(y)) {
  interval <- smuce_theta_interval(y, gamma, sigma2, n)
  if (anyNA(interval)) {
    return(Inf)
  }
  offset <- y[1L]
  y_centered <- y - offset
  theta <- min(max(mean(y), interval[1]), interval[2])
  sum((y_centered - (theta - offset))^2) / sigma2
}

#' SVP with SMUCE validity and constrained Gaussian cost
#' @param y Numeric observations.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @return Integer segment-end indices. Internal segment endpoints are the
#' estimated changepoints; the final endpoint is `length(y)`.
#' @details
#' For independent Gaussian observations with known variance, this function
#' computes the SMUCE dynamic program: it first minimizes the number of
#' constant segments satisfying the multiscale constraint, then minimizes the
#' constrained Gaussian residual cost among those partitions. The function
#' returns one optimal partition as its segment endpoints, not fitted levels
#' or confidence intervals. Exact changepoint locations are not unique when
#' multiple optimal partitions have the same objective; ties are resolved
#' deterministically by the dynamic program.
#' @references
#' Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point
#' Inference. *Journal of the Royal Statistical Society: Series B*, 76(3),
#' 495--580. doi:10.1111/rssb.12047.
#' @examples
#' set.seed(1)
#' data <- ts_generator(
#'   chpts = c(20, 40), parameters = c(0, 2),
#'   sd_noise = 1, type = "gauss"
#' )
#' svp_smuce(data, gamma = 1.5, sigma2 = 1)
#' @export
svp_smuce <- function(y, gamma, sigma2 = 1) {
  n <- length(y)
  .smuce_validate(y, gamma, sigma2, n)
  K <- rep(n + 1L, n + 1L)
  C <- rep(Inf, n + 1)
  prev <- rep(-1L, n + 1L)
  K[1] <- 0
  C[1] <- 0
  # Center before prefix sums to make the recurrence translation invariant.
  offset <- y[1L]
  y_centered <- y - offset
  cs <- c(0, cumsum(y_centered))
  cs2 <- c(0, cumsum(y_centered^2))

  # The admissible interval for a fixed start can only shrink as the segment
  # grows. Updating it incrementally reduces this implementation from O(n^4)
  # to O(n^3), while preserving the SMUCE dynamic program.
  for (s in 0:(n - 1L)) {
    if (!is.finite(C[s + 1L])) {
      next
    }
    lo <- -Inf
    hi <- Inf
    for (t in s:(n - 1L)) {
      for (u in s:t) {
        m <- t - u + 1L
        centre <- (cs[t + 2L] - cs[u + 1L]) / m
        radius <- sqrt(sigma2 / m) *
          (gamma + sqrt(2 * (1 + log(n / m))))
        lo <- max(lo, centre - radius)
        hi <- min(hi, centre + radius)
      }
      if (lo > hi) {
        break
      }
      m <- t - s + 1L
      mean_seg <- (cs[t + 2L] - cs[s + 1L]) / m
      theta <- min(max(mean_seg, lo), hi)
      rss <- (cs2[t + 2L] - cs2[s + 1L]) -
        2 * theta * (cs[t + 2L] - cs[s + 1L]) + m * theta^2
      cc <- max(0, rss) / sigma2
      kk <- K[s + 1L] + 1L
      val <- C[s + 1L] + cc
      if (kk < K[t + 2L] || (kk == K[t + 2L] && val < C[t + 2L])) {
        K[t + 2L] <- kk
        C[t + 2L] <- val
        prev[t + 2L] <- s
      }
    }
  }

  if (prev[n + 1L] < 0L) {
    stop("no valid SMUCE partition")
  }
  ends <- integer(0)
  t <- n
  while (t > 0L) {
    ends <- c(t, ends)
    t <- prev[t + 1L]
  }
  as.integer(ends)
}
