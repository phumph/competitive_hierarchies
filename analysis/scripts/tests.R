library(testthat)
source("analysis/scripts/comp_utils.R")

test_that("Matrix calculations yield correct output", {
  input <- matrix(c(
    0, 0, 0, 1,
    1, 0, 1, 0,
    1, 0, 0, 0,
    0, 1, 1, 0
  ), nrow = 4)

  ivec <- c(0, 1, 0, 0)

  expected_matix <- matrix(c(
    0, 1, 0, 1,
    0, 0, 0, 0,
    1, 1, 0, 0,
    0, 1, 1, 0
  ), nrow = 4)

  output <- knockout_by_ivec(input, ivec)
  testthat::expect_equal(output, expected_matix)
})
