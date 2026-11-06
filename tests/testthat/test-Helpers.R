test_that("temp2x and x2temp convert both ways", {
  a <- temp2x(c(5, 25), tempError = c(0.5, 1))
  expect_equal(a$Temperature, 10^6 / (c(5, 25) + 273.15)^2)
  b <- x2temp(a$Temperature, xError = a$TempError)
  expect_equal(b$Temp, c(5, 25))
  expect_equal(b$TempSD, c(0.5, 1))
  expect_equal(temp2x(298.15, unit = "K")$Temperature, temp2x(25)$Temperature)
  expect_error(temp2x(-300), "absolute zero")
})

test_that("avg.replicates returns means and standard errors", {
  reps <- data.frame(Sample = rep(c("A", "B"), each = 3), D47 = c(0.601, 0.598, 0.605, 0.702, 0.695, 0.699), Material = "calcite", Run = 1:6)
  a <- avg.replicates(reps)
  expect_equal(a$D47, c(mean(reps$D47[1:3]), mean(reps$D47[4:6])))
  expect_equal(a$D47error[1], sd(reps$D47[1:3]) / sqrt(3))
  expect_equal(a$N, c(3, 3))
  expect_true("Material" %in% names(a))
  expect_false("Run" %in% names(a))
  expect_equal(avg.replicates(reps, pooledSD = 0.008)$D47error, rep(0.008 / sqrt(3), 2))
})

test_that("avg.replicates warns about single replicates", {
  reps <- data.frame(Sample = c("A", "A", "B"), D47 = c(0.6, 0.61, 0.7))
  expect_warning(avg.replicates(reps), "single replicate")
})

test_that("cal.ci returns a data.frame and respects the level", {
  ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
  ols <- cal.ols(ds, replicates = 200)
  a <- cal.ci(ols, from = 10, to = 14)
  b <- cal.ci(ols, from = 10, to = 14, level = 0.5)
  expect_s3_class(a, "data.frame")
  expect_equal(colnames(a), c("x", "median_est", "ci_lower_est", "ci_upper_est"))
  expect_true(all(b$ci_upper_est - b$ci_lower_est < a$ci_upper_est - a$ci_lower_est))
})

test_that("cal.summary reports estimates and standard errors", {
  ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
  ols <- cal.ols(ds, replicates = 200)
  a <- cal.summary(ols)
  expect_equal(a$parameter, c("alpha", "beta", "sigma"))
  expect_equal(colnames(a), c("parameter", "estimate", "SE", "median", "lower", "upper", "ESS", "Rhat"))
  expect_true(all(is.na(a$ESS)))
  expect_equal(a$SE[2], sd(ols$beta))

  # Mixed draws are summarized by material
  draws <- rbind(cbind(ols, Material = "calcite"), cbind(ols, Material = "aragonite"))
  b <- cal.summary(draws)
  expect_equal(nrow(b), 6)
  expect_equal(unique(b$Material), c("calcite", "aragonite"))
})
