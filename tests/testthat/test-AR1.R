library(testthat)
library(svpChange)

simulate_ar1_change <- function(n, changepoint, means, rho, innovation_sd = 1) {
  innovations <- rnorm(n, sd = innovation_sd)
  mu <- rep(means, c(changepoint, n - changepoint))
  data <- numeric(n)
  data[1] <- mu[1] + innovations[1] / sqrt(1 - rho^2)
  for (i in 2:n) {
    data[i] <- mu[i] + rho * (data[i - 1] - mu[i - 1]) + innovations[i]
  }
  data
}

exact_ar1_reference <- function(data, rho, sigma2 = 1, profile_sigma = FALSE) {
  n <- length(data)
  z <- data[-1] - rho * data[-n]
  rss0 <- sum((z - mean(z))^2)
  taus <- 2:(n - 2)
  rss1 <- vapply(taus, function(tau) {
    design <- matrix(0, nrow = n - 1, ncol = 2)
    design[seq_len(tau - 1), 1] <- 1 - rho
    design[tau, ] <- c(-rho, 1)
    if (tau + 1 <= n - 1) {
      design[(tau + 1):(n - 1), 2] <- 1 - rho
    }
    sum(lm.fit(design, z)$residuals^2)
  }, numeric(1))
  best <- which.min(rss1)
  statistic <- if (profile_sigma) {
    0.5 * (n - 1) * log(rss0 / rss1[best])
  } else {
    (rss0 - rss1[best]) / (2 * sigma2)
  }
  list(
    statistic = statistic,
    changepoint = taus[best],
    rss0 = rss0,
    rss1 = rss1[best]
  )
}

test_that("robust AR1 estimate matches equation 3", {
  set.seed(1)
  data <- simulate_ar1_change(301, 150, c(0, 1.5), 0.6)
  lag1_median <- median(abs(diff(data)))
  lag2_median <- median(abs(data[3:length(data)] - data[1:(length(data) - 2)]))
  reference <- lag2_median^2 / lag1_median^2 - 1

  expect_equal(AR1_rho(data), max(-0.999, min(0.999, reference)), tolerance = 1e-12)
})

test_that("AR1 single-change test locates a simulated mean change", {
  set.seed(42)
  data <- simulate_ar1_change(400, 200, c(0, 2), 0.7)
  result <- AR1_single_change(data, gamma = 10, rho = 0.7)

  expect_false(result$valid)
  expect_lte(abs(result$changepoint - 200), 15)
})

test_that("AR1 statistic equals the exact fixed-rho reference", {
  set.seed(7)
  data <- simulate_ar1_change(120, 53, c(-0.5, 1.25), 0.65)
  reference <- exact_ar1_reference(data, rho = 0.65, sigma2 = 1)
  result <- AR1_single_change(data, gamma = Inf, rho = 0.65, sigma2 = 1)

  expect_equal(result$statistic, reference$statistic, tolerance = 1e-9)
  expect_equal(result$rss0, reference$rss0, tolerance = 1e-9)
  expect_equal(result$rss1, reference$rss1, tolerance = 1e-9)
  expect_equal(result$changepoint, reference$changepoint)
})

test_that("AR1 profiled statistic equals the exact reference", {
  set.seed(8)
  data <- simulate_ar1_change(100, 47, c(0, 1), 0.4, innovation_sd = 1.3)
  reference <- exact_ar1_reference(data, rho = 0.4, profile_sigma = TRUE)
  result <- AR1_single_change(
    data, gamma = Inf, rho = 0.4, profile_sigma = TRUE
  )

  expect_equal(result$statistic, reference$statistic, tolerance = 1e-9)
  expect_equal(result$changepoint, reference$changepoint)
})

test_that("profiled AR1 diagnostics and R validity ignore sigma2", {
  set.seed(9)
  data <- simulate_ar1_change(80, 40, c(0, 1), 0.4)

  default <- AR1_single_change(
    data, gamma = Inf, rho = 0.4, sigma2 = 1,
    profile_sigma = TRUE
  )
  zero <- AR1_single_change(
    data, gamma = Inf, rho = 0.4, sigma2 = 0,
    profile_sigma = TRUE
  )

  expect_equal(zero$statistic, default$statistic)
  expect_equal(zero$changepoint, default$changepoint)
  expect_identical(
    valid_AR1(data, gamma = 5, rho = 0.4, sigma2 = 0,
              profile_sigma = TRUE),
    valid_AR1(data, gamma = 5, rho = 0.4, sigma2 = 1,
              profile_sigma = TRUE)
  )
  expect_error(
    AR1_single_change(data, gamma = 5, rho = 0.4, sigma2 = 0),
    "sigma2"
  )
  expect_error(
    valid_AR1(data, gamma = 5, rho = 0.4, sigma2 = 0),
    "sigma2"
  )
})

test_that("AR1 validity test works inside SVP", {
  set.seed(42)
  data <- simulate_ar1_change(400, 200, c(0, 2), 0.7)
  result <- SVP(data, gamma = 10, test = "AR1", rho = 0.7)

  expect_equal(result$rho, 0.7)
  expect_equal(tail(result$changepoints, 1), length(data))
  expect_true(any(abs(head(result$changepoints, -1) - 200) <= 10))
})

test_that("AR1 validity test works inside main SVP with pruning options", {
  set.seed(43)
  data <- simulate_ar1_change(400, 200, c(0, 2), 0.7)
  result <- SVP(
    data,
    gamma = 10,
    test = "AR1",
    subtests = "both",
    rho = 0.7,
    sigma2 = 1
  )

  expect_equal(tail(result$changepoints, 1), length(data))
  expect_true(any(abs(head(result$changepoints, -1) - 200) <= 10))
})

test_that("AR1Focus matches the conditional AR1 SVP reference", {
  exact_validity <- function(rho, sigma2 = 1) {
    force(rho)
    force(sigma2)
    function(segment, gamma) {
      AR1_single_change(
        segment, gamma = gamma, rho = rho, sigma2 = sigma2,
        profile_sigma = FALSE
      )$valid
    }
  }

  for (rho in c(0, 0.3, 0.7, -0.4)) {
    set.seed(1)
    data <- simulate_ar1_change(300, 90, c(0, 1.5), rho)
    exact <- SVP(
    data, gamma = 8, test = "AR1", rho = rho, sigma2 = 1,
      subtests = "right"
    )
    focus <- SVP(
    data, gamma = 8, test = "AR1Focus", rho = rho, sigma2 = 1,
      subtests = "right"
    )
    reference <- svp0(
      data, gamma = 8, test = exact_validity(rho),
      subtests = "right",
      PELT_pruning = FALSE
    )

    # svp0() always ranks equal-length partitions with the Gaussian cost, so
    # the boundary comparison against it uses SVP() on the same cost scale.
    focus_gaussian <- SVP(
      data, gamma = 8, test = "AR1Focus", rho = rho, sigma2 = 1,
      subtests = "right", cost = "gaussian"
    )

    expect_equal(focus$changepoints, exact$changepoints)
    expect_equal(focus_gaussian$changepoints, reference$changepoints)

    # The cost never changes how many segments the validity test allows.
    expect_equal(
      length(focus$changepoints), length(reference$changepoints)
    )
  }
})

test_that("AR1Focus matches the conditional reference on randomized streams", {
  for (rho in c(-0.8, -0.4, 0, 0.4, 0.8, 0.95)) {
    for (seed in seq_len(5L)) {
      set.seed(1000 + seed)
      n <- 80L
      changepoint <- sample(15:65, 1L)
      data <- simulate_ar1_change(
        n, changepoint, c(0, 1.5), rho
      )
      gamma <- c(4, 6, 8, 10)[(seed %% 4L) + 1L]
      exact <- SVP(
        data, gamma = gamma, test = "AR1", rho = rho, sigma2 = 1,
        subtests = "right", cost = "gaussian"
      )
      focus <- SVP(
        data, gamma = gamma, test = "AR1Focus", rho = rho, sigma2 = 1,
        subtests = "right", cost = "gaussian"
      )
      expect_equal(focus$changepoints, exact$changepoints)
    }
  }
})

test_that("the AR(1) cost sums the interior innovations of each segment", {
  # Independent recomputation of the documented cost: within each segment
  # (s, t], the sum of squared deviations of z[u] = y[u] - rho * y[u - 1] from
  # their mean, taken over u = s + 2, ..., t.
  ar1_partition_cost <- function(y, boundaries, rho) {
    starts <- c(0L, head(boundaries, -1L))
    sum(vapply(seq_along(boundaries), function(k) {
      first <- starts[k] + 2L
      last <- boundaries[k]
      if (last < first) return(0)
      z <- y[first:last] - rho * y[(first - 1L):(last - 1L)]
      sum((z - mean(z))^2)
    }, numeric(1)))
  }

  for (rho in c(0, 0.5, 0.85, -0.6)) {
    set.seed(11)
    data <- simulate_ar1_change(240, 120, c(0, 2), rho)
    fit <- SVP(
      data, gamma = 8, test = "AR1Focus", rho = rho, sigma2 = 1,
      subtests = "right"
    )
    expect_equal(
      fit$R[length(data), 1L],
      ar1_partition_cost(data, fit$changepoints, rho)
    )
  }
})

test_that("cost = 'gaussian' keeps the previous squared-error ranking", {
  gaussian_partition_cost <- function(y, boundaries) {
    starts <- c(1L, head(boundaries, -1L) + 1L)
    sum(vapply(seq_along(boundaries), function(k) {
      segment <- y[starts[k]:boundaries[k]]
      sum((segment - mean(segment))^2)
    }, numeric(1)))
  }

  set.seed(12)
  data <- simulate_ar1_change(240, 120, c(0, 2), 0.7)
  fit <- SVP(
    data, gamma = 8, test = "AR1Focus", rho = 0.7, sigma2 = 1,
    subtests = "right", cost = "gaussian"
  )
  expect_equal(
    fit$R[length(data), 1L], gaussian_partition_cost(data, fit$changepoints)
  )

  # The default for a non-AR(1) test is still the Gaussian cost.
  set.seed(13)
  independent <- ts_generator(
    chpts = c(60, 120), parameters = c(0, 2), sd_noise = 1, type = "gauss"
  )
  expect_equal(
    SVP(independent, gamma = 8)$changepoints,
    SVP(independent, gamma = 8, cost = "gaussian")$changepoints
  )
})

test_that("the AR(1) cost localises a strongly correlated change better", {
  # Under strong serial dependence the Gaussian cost is misspecified, so it
  # places the boundary of an otherwise identical partition less accurately.
  rho <- 0.85
  errors <- vapply(seq_len(120), function(seed) {
    set.seed(2000 + seed)
    data <- simulate_ar1_change(400, 200, c(0, 4), rho)
    vapply(c("ar1", "gaussian"), function(cost) {
      changepoints <- SVP(
        data, gamma = 4 * log(400), test = "AR1Focus", rho = rho,
        sigma2 = 1, subtests = "right", cost = cost
      )$changepoints
      internal <- changepoints[changepoints < 400]
      if (length(internal) != 1L) return(NA_real_)
      abs(internal - 200)
    }, numeric(1))
  }, numeric(2))

  usable <- !is.na(errors[1L, ]) & !is.na(errors[2L, ])
  expect_gt(sum(usable), 40L)
  expect_lt(mean(errors[1L, usable]), mean(errors[2L, usable]))
  expect_gt(
    sum(errors[1L, usable] < errors[2L, usable]),
    sum(errors[1L, usable] > errors[2L, usable])
  )
})

test_that("cost accepts only the documented values", {
  set.seed(14)
  data <- simulate_ar1_change(120, 60, c(0, 2), 0.5)
  expect_error(SVP(data, gamma = 8, cost = "ar2"), "Unknown cost type")
})
