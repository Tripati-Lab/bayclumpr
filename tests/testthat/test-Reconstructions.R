# Draws with fixed parameters, so results can be checked against known values
fixedDraws <- function(sigma) {
  data.frame(alpha = rep(0.268, 100), beta = rep(0.0369, 100), sigma = rep(sigma, 100))
}

# D47 of a sample at 15 C
x15 <- 10^6 / (273.15 + 15)^2
d47 <- 0.268 + 0.0369 * x15

test_that("rec.bayesian recovers a known temperature", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 1e-6)
  a <- rec.bayesian(fixedDraws(1e-6), recData, ndraws = 1000)
  expect_equal(a$meanTemp, 15, tolerance = 1e-3)
})

test_that("rec.bayesian SEs match the propagated uncertainty", {
  # First-order SD of temperature for sigma = 0.01: 3.24 C
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0)
  a <- rec.bayesian(fixedDraws(0.01), recData, ndraws = 1e5)
  expect_equal(a$SE, 3.24, tolerance = 0.02)
})

test_that("Medians and intervals are consistent", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0)
  a <- rec.bayesian(fixedDraws(0.01), recData, ndraws = 1e5)
  b <- rec.bayesian(fixedDraws(0.01), recData, ndraws = 1e5, level = 0.5)
  expect_equal(a$median, 15, tolerance = 0.01)
  expect_true(a$lower < a$median && a$median < a$upper)
  expect_true(b$upper - b$lower < a$upper - a$lower)
})

test_that("Bayesian SEs are larger than SEs without residual scatter", {
  set.seed(1)
  ds <- cal.dataset(error = "S1", nobs = 50)
  boot <- cal.ols(data = ds, replicates = 200)
  recData <- data.frame(Sample = c("A", "B"), D47 = c(0.6, 0.7), D47error = c(0.005, 0.02))
  bay <- rec.bayesian(boot, recData)
  freq <- rec.clumped(recData, boot)
  expect_true(all(bay$SE > freq$SE))
})

test_that("rec.clumped with prediction = TRUE includes sigma", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0.005)
  a <- rec.clumped(recData, fixedDraws(0.01), prediction = FALSE)
  b <- rec.clumped(recData, fixedDraws(0.01), prediction = TRUE)
  expect_gt(b$SE, a$SE)
  expect_equal(b$meanTemp, 15, tolerance = 1e-6)
  expect_error(rec.clumped(recData, fixedDraws(0.01)[, 1:2], prediction = TRUE), "sigma")
})

test_that("rec.bayesian works with a single sample", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0.01)
  expect_equal(nrow(rec.bayesian(fixedDraws(0.01), recData)), 1)
})

test_that("rec.bayesian matches materials in mixed models", {
  draws <- rbind(cbind(fixedDraws(1e-6), Material = "calcite"), cbind(fixedDraws(1e-6), Material = "aragonite"))
  draws$alpha[draws$Material == "aragonite"] <- 0.288
  recData <- data.frame(Sample = c("A", "B"), D47 = c(d47, d47 + 0.02), D47error = 1e-6, Material = c("calcite", "aragonite"))
  a <- rec.bayesian(draws, recData, ndraws = 1000)
  expect_equal(a$meanTemp, c(15, 15), tolerance = 1e-3)

  recData$Material <- c("calcite", "dolomite")
  expect_error(rec.bayesian(draws, recData), "Materials")
})

test_that("rec.bayesian stops when the calibration is not in 10^6/T^2", {
  draws <- fixedDraws(0.01)
  draws$beta <- -0.004
  recData <- data.frame(Sample = "A", D47 = 0.6, D47error = 0.01)
  expect_error(rec.bayesian(draws, recData), "10\\^6/T\\^2")
})

test_that("rec.bayesian warns about arguments from earlier versions", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0.01)
  expect_warning(rec.bayesian(fixedDraws(0.01), recData, iter = 1000), "iter")
})

test_that("rec.ion.bayesian removes the effect of the covariate", {
  draws <- cbind(fixedDraws(1e-6), gamma = 0.01)
  recData <- data.frame(Sample = "A", D47 = d47 + 0.01 * 2, D47error = 1e-6, Ion = 2, IonError = 1e-6)
  a <- rec.ion.bayesian(draws, recData, ndraws = 1000)
  b <- rec.ion.bayesian(draws, recData, ndraws = 1000, useIonError = TRUE)
  expect_equal(a$meanTemp, 15, tolerance = 1e-3)
  expect_equal(b$meanTemp, 15, tolerance = 1e-3)
})

test_that("rec.clumped SEs include the calibration spread", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0)
  boot <- data.frame(alpha = c(0.26, 0.27, 0.28), beta = c(0.0375, 0.0369, 0.0363))
  expect_gt(rec.clumped(recData, boot)$SE, 0)
})

test_that("SD uses the SD of a single replicate and SE uses the SE of the mean", {
  # With fixed parameters and no residual scatter, SD / SE is close to sqrt(N)
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0.005, N = 4)
  a <- rec.bayesian(fixedDraws(1e-6), recData, ndraws = 1e5)
  expect_equal(a$SD / a$SE, 2, tolerance = 0.02)

  # D47sd gives the same SD when N is missing
  recSD <- data.frame(Sample = "A", D47 = d47, D47error = 0.005, D47sd = 0.01)
  b <- rec.bayesian(fixedDraws(1e-6), recSD, ndraws = 1e5)
  expect_equal(b$SD, a$SD, tolerance = 0.02)
})

test_that("SD is NA without replicate information", {
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0.005)
  expect_true(is.na(rec.clumped(recData, fixedDraws(0.01))$SD))
})

test_that("SD and SE both include the calibration uncertainty", {
  # With one replicate, SD and SE describe the same quantity
  recData <- data.frame(Sample = "A", D47 = d47, D47error = 0.005, N = 1)
  a <- rec.clumped(recData, fixedDraws(0.01), ndraws = 1e5, prediction = TRUE)
  expect_equal(a$SD, a$SE, tolerance = 0.02)
})
