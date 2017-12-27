library(DEButilities)
context("numerical equivalenmce to DEBtool_M")
#crude test to compare example function calls with the corresponding DEBtool_M output

test_that("get_tb matches DEBtool", {
  expect_equal(get_tb(c(.1,.5,.03)), c(8.6038, 0.2658, 1), tolerance = 0.001)
})

