# Independent reference: enumerate all triples using direct sample means.
verzelen_oracle <- function(y, sigma2 = 1, N = length(y)) {
  best <- -Inf
  if (length(y) < 2L) return(best)
  for (a in 0:(length(y) - 2L)) {
    for (c in (a + 2L):length(y)) {
      for (b in (a + 1L):(c - 1L)) {
        cusum2 <- (b - a) * (c - b) / (c - a) *
          (mean(y[(a + 1L):b]) - mean(y[(b + 1L):c]))^2
        best <- max(best, cusum2 / (4 * sigma2) -
                      2 * log(N * (c - a) / ((b - a) * (c - b))))
      }
    }
  }
  best
}

test_that("Verzelen scans all triples with fixed global length", {
  stat <- svpChange:::.verzelen_statistic_cpp
  set.seed(246)
  series <- c(list(7, rep(2, 9), c(0, 0, 20, -20, 0, 0),
                   seq(-4, 4), rep(c(-2, 2), 5)),
              replicate(20, rnorm(12), simplify = FALSE))
  for (y in series) for (variance in c(0.4, 1, 3)) {
    expected <- verzelen_oracle(y, variance, N = 20)
    expect_equal(stat(y, variance, 20L), expected, tolerance = 1e-11)
    expect_equal(stat(y, variance, 20L, FALSE), expected, tolerance = 1e-11)
  }
  # The maximum here is strictly interior: neither endpoint is in its triple.
  y <- c(0, 0, 20, -20, 0, 0)
  expect_equal(stat(y), 200 - 2 * log(12))
  expect_equal(stat(rep(0, 8)), -2 * log(4))
  expect_equal(stat(y, N = 60L), stat(y) - 2 * log(10))
})

test_that("Verzelen bounds agree with exhaustive scans on longer series", {
  stat <- svpChange:::.verzelen_statistic_cpp
  set.seed(432)
  for (y in list(rnorm(140), c(rnorm(60), rnorm(80, 4)),
                 rep(0, 140), (1:140)^2 / 100, rep(c(-3, 3), 70))) {
    expect_equal(stat(y), stat(y, accelerated = FALSE), tolerance = 1e-10)
    expect_equal(stat(y + 1e6), stat(y), tolerance = 1e-7)
    expect_equal(stat(-y), stat(y), tolerance = 1e-10)
    expect_equal(stat(3 * y, sigma2 = 9), stat(y), tolerance = 1e-10)
  }
})

test_that("Verzelen SVP matches unpruned reference dynamic programming", {
  set.seed(52)
  for (y in c(list(rep(0, 8), c(0, 0, 7, 7, -5, -5)),
              replicate(10, rnorm(9, sd = 3), simplify = FALSE))) {
    for (gamma in c(-8, -3, 0, 2)) {
      valid <- function(segment, gamma)
        verzelen_oracle(segment, N = length(y)) <= gamma
      fit <- SVP(y, gamma, subtests = "verzelen")
      reference <- svp0(y, gamma, test = valid, subtests = "none",
                        PELT_pruning = FALSE)
      expect_equal(as.numeric(fit$changepoints),
                   as.numeric(reference$changepoints))
      expect_equal(unname(fit$R), unname(reference$R), tolerance = 1e-10)
      pruned <- svp0(y, gamma, test = valid, subtests = "both",
                    PELT_pruning = FALSE)
      expect_equal(as.numeric(fit$lastIndexSet), as.numeric(pruned$lastIndexSet))
      expect_equal(as.numeric(fit$nb), as.numeric(pruned$nb))
      expect_null(fit$costQ)
      critical <- svpChange:::.verzelen_statistic_cpp(y)
      expect_identical(length(fit$changepoints) == 1L, critical <= gamma)
    }
  }
})

test_that("Verzelen handles equality, variance and invalid arguments", {
  stat <- svpChange:::.verzelen_statistic_cpp
  y <- c(0, 8)
  threshold <- stat(y)
  expect_length(SVP(y, threshold, subtests = "verzelen")$changepoints, 1L)
  expect_length(SVP(y, threshold - 1e-10, subtests = "verzelen")$changepoints, 2L)
  expect_equal(as.numeric(SVP(4, -100, subtests = "verzelen")$changepoints), 1)
  expect_equal(SVP(3 * y, 0, sigma2 = 9, subtests = "verzelen")$changepoints,
               SVP(y, 0, subtests = "verzelen")$changepoints)
  for (g in c(NA_real_, Inf, -Inf))
    expect_error(SVP(y, g, subtests = "verzelen"), "finite")
  for (s in c(0, -1, NA_real_, Inf))
    expect_error(SVP(y, 1, sigma2 = s, subtests = "verzelen"), "sigma2")
  expect_error(SVP(y, -1), "positive")
  expect_error(SVP(y, 1, test = "AR1", subtests = "verzelen"), "gaussian_mean")
  expect_error(SVP(c(NA, 1), 1, subtests = "verzelen"), "finite")
  expect_error(stat(y, N = 1L), "at least")
})
