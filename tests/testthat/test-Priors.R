test_that("cal.prior uses the same values as cal.bayesian", {
  a <- cal.prior("Informative", n = 10)
  b <- cal.prior("Weak", n = 10)
  expect_equal(attr(a, "params")$sd, c(0.065, 0.004))
  expect_equal(attr(b, "params")$sd, c(1, 1))
})

test_that("cal.prior accepts custom priors", {
  a <- cal.prior(list(alpha = c(0.2, 0.01), beta = c(0.04, 0.001)), n = 10)
  expect_equal(attr(a, "params")$mean, c(0.2, 0.04))
  expect_error(cal.prior(list(alpha = c(0.2, -1), beta = c(0.04, 0.001))), "positive")
})

test_that("Flat priors can't be sampled", {
  expect_error(cal.prior("Flat"), "Flat")
  expect_error(cal.prior("Uninformative"), "Flat")
})

test_that("rec.prior adds temperatures", {
  a <- rec.prior("Informative", n = 10)
  expect_equal(colnames(a), c("alpha", "beta", "x", "xC"))
})
