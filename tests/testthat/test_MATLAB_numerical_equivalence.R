library(DEButilities)
context("numerical equivalence to DEBtool_M")
#crude test to compare example function calls with the corresponding DEBtool_M output

test_that("get_tb matches DEBtool", {
  expect_equal(get_tb(c(.1,.5,.03)), c(8.6038, 0.2658, 1), tolerance = 0.001)
})

test_that("beta0 matches DEBtool_M", {
  expect_equal(beta0(0.1, 0.2), 0.062455351015153)
})

test_that("get_ue0 matches DEBtool_M", {
  g = 6; k = 6; kap = .8; uHb = .001; vHb = uHb/ (1 - kap);
  pars = c(g, k, vHb);
  expect_equal(get_ue0(pars), c(0.008171217468332, 0.188067711870661, 1))
})
