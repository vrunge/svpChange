.msl_interval_statistic <- function(y, start, end, sigma2 = 1) {
  values <- y[(start + 1L):end]
  m <- length(values)
  if (m < 2L) return(-Inf)
  prefix <- c(0, cumsum(values))
  max(vapply(seq_len(m - 1L), function(left) {
    right <- m - left
    z <- abs(prefix[left + 1L] / left -
               (prefix[m + 1L] - prefix[left + 1L]) / right) *
      sqrt(left * right / m) / sqrt(sigma2)
    z - sqrt(2 * log(exp(1) * m / min(left, right)))
  }, numeric(1)))
}

.msl_segment_valid <- function(y, start, end, q, sigma2 = 1) {
  if (end - start < 2L) return(TRUE)
  for (a in start:(end - 2L)) {
    for (b in (a + 2L):end) {
      if (.msl_interval_statistic(y, a, b, sigma2) > q) return(FALSE)
    }
  }
  TRUE
}

.msl_brute_partition <- function(y, q, sigma2 = 1) {
  n <- length(y)
  best <- NULL
  best_cost <- Inf
  for (mask in 0:(2^(n - 1L) - 1L)) {
    splits <- which(as.logical(intToBits(mask)[seq_len(n - 1L)]))
    ends <- c(splits, n)
    starts <- c(0L, head(ends, -1L))
    valid <- all(vapply(seq_along(ends), function(i) {
      .msl_segment_valid(y, starts[i], ends[i], q, sigma2)
    }, logical(1)))
    if (!valid) next
    cost <- sum(vapply(seq_along(ends), function(i) {
      x <- y[(starts[i] + 1L):ends[i]]
      sum((x - mean(x))^2)
    }, numeric(1)))
    if (is.null(best) || length(ends) < length(best) ||
        (length(ends) == length(best) && cost < best_cost)) {
      best <- ends
      best_cost <- cost
    }
  }
  best
}

test_that("MS-l SVP matches exhaustive hereditary validity and SSE ranking", {
  set.seed(943)
  for (i in seq_len(20L)) {
    y <- rnorm(7L) + rep(c(0, 2), c(3, 4))
    sigma2 <- if (i %% 2L) 1 else 2.25
    for (q in c(-0.5, 0.5, 1.5)) {
      expected <- .msl_brute_partition(y, q, sigma2)
      actual <- SVP(y, gamma = q, subtests = "both_msl",
                    sigma2 = sigma2)$changepoints
      expect_equal(actual, expected)
    }
  }
})

test_that("MS-l threshold uses the maximum over all contained intervals", {
  y <- c(0, 0, 4, 4, 0, 0)
  critical <- max(unlist(lapply(0:(length(y) - 2L), function(a) {
    vapply((a + 2L):length(y), function(b) {
      .msl_interval_statistic(y, a, b)
    }, numeric(1))
  })))
  expect_length(SVP(y, critical + 1e-7, subtests = "both_msl")$changepoints,
                1L)
  expect_gt(length(SVP(y, critical - 1e-7,
                       subtests = "both_msl")$changepoints), 1L)
})

test_that("MS-l handles known variance, shifts, signs, and inputs", {
  set.seed(944)
  y <- rnorm(30)
  expected <- SVP(y, 0.5, subtests = "both_msl")$changepoints
  expect_equal(SVP(y + 1e5, 0.5, subtests = "both_msl")$changepoints,
               expected)
  expect_equal(SVP(-y, 0.5, subtests = "both_msl")$changepoints,
               expected)
  expect_equal(SVP(2 * y, 0.5, subtests = "both_msl",
                   sigma2 = 4)$changepoints, expected)
  expect_equal(SVP(3, -10, subtests = "both_msl")$changepoints, 1)
  expect_error(SVP(1:4, Inf, subtests = "both_msl"), "finite")
  expect_error(SVP(1:4, 1, subtests = "both_msl", sigma2 = 0),
               "positive")
  expect_error(SVP(1:4, 1, test = "gamma_rate",
                   subtests = "both_msl"), "gaussian_mean")
})
