test_that("singleton segments are always valid", {
  tested_lengths <- integer()
  reject_long_segments <- function(segment, gamma) {
    tested_lengths <<- c(tested_lengths, length(segment))
    FALSE
  }

  data <- c(1, 2, 3)
  result <- svp0(
    data,
    gamma = 0,
    test = reject_long_segments,
    subtests = "none",
    PELT_pruning = FALSE
  )

  expect_equal(as.numeric(result$changepoints), as.numeric(data))
  expect_false(any(tested_lengths == 1L))
  expect_equal(result$R[, 2], as.numeric(seq_along(data)))
  expect_equal(result$R[, 3], as.numeric(seq_along(data) - 1L))
  expect_null(result$costQ)
})

test_that("svp0 returns the same partition for every pruning option", {
  range_test <- function(segment, gamma) {
    diff(range(segment)) <= gamma
  }
  data <- c(rep(0, 3), rep(3, 3))
  expected <- c(3, 6)

  for (subtests in c("none", "right", "both")) {
    for (PELT_pruning in c(FALSE, TRUE)) {
      result <- svp0(
        data,
        gamma = 0,
        test = range_test,
        subtests = subtests,
        PELT_pruning = PELT_pruning
      )
      expect_equal(as.numeric(result$changepoints), expected)
    }
  }
})

test_that("svp0 validates subtests and rejects empty data", {
  valid <- function(segment, gamma) TRUE

  expect_error(
    svp0(1:3, gamma = 1, test = valid, subtests = "invalid"),
    "Invalid value for 'subtests'"
  )
  expect_error(
    svp0(1:3, gamma = 1, test = valid, subtests = "left"),
    "both.*right.*none"
  )
  expect_error(
    svp0(numeric(), gamma = 1, test = valid),
    "at least one observation"
  )
})
