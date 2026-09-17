simulate_validity_ar1 <- function(n, changepoint, means, rho,
                                  innovation_sd = 1) {
  innovations <- rnorm(n, sd = innovation_sd)
  mu <- rep(means, c(changepoint, n - changepoint))
  data <- numeric(n)
  data[1] <- mu[1] + innovations[1] / sqrt(1 - rho^2)
  for (i in 2:n) {
    data[i] <- mu[i] + rho * (data[i - 1] - mu[i - 1]) + innovations[i]
  }
  data
}

test_that("valid_AR1 agrees with AR1_single_change", {
  set.seed(19)
  data <- simulate_validity_ar1(80, 35, c(0, 1.5), rho = 0.6)

  known_variance <- AR1_single_change(data, gamma = Inf, rho = 0.6)
  for (gamma in c(0, known_variance$statistic / 2,
                  known_variance$statistic + 1, Inf)) {
    expect_identical(
      valid_AR1(data, gamma = gamma, rho = 0.6),
      AR1_single_change(data, gamma = gamma, rho = 0.6)$valid
    )
  }

  profiled <- AR1_single_change(
    data, gamma = Inf, rho = 0.6, profile_sigma = TRUE
  )
  for (gamma in c(0, profiled$statistic / 2,
                  profiled$statistic + 1, Inf)) {
    expect_identical(
      valid_AR1(data, gamma = gamma, rho = 0.6, profile_sigma = TRUE),
      AR1_single_change(
        data, gamma = gamma, rho = 0.6, profile_sigma = TRUE
      )$valid
    )
  }
})

test_that("valid_AR1 supports robust rho estimation", {
  set.seed(20)
  data <- simulate_validity_ar1(80, 35, c(0, 1.5), rho = 0.6)

  expect_identical(
    valid_AR1(data, gamma = 5),
    AR1_single_change(data, gamma = 5)$valid
  )
})

test_that("valid_SCALE uses the normalized MAD", {
  data <- seq(-2, 2, length.out = 9)
  expected_scale <- stats::mad(data, constant = 1.4826)

  expect_true(valid_SCALE(data, gamma = expected_scale))
  expect_false(valid_SCALE(data, gamma = expected_scale - 0.001))
  expect_true(valid_SCALE(rep(0, 10), gamma = 0))
})
