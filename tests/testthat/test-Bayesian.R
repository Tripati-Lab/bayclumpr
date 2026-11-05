# These tests fit Stan models, so they are slow. Run them with Sys.setenv(RUN_STAN_TESTS = "true").
skipStan <- function() {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RUN_STAN_TESTS"), "true"), "Set RUN_STAN_TESTS = true to run Stan tests.")
}

test_that("The errors model separates residual scatter from D47 error", {
  skipStan()
  # Scenario S3: residual SD 0.0275 and D47 error 0.0125
  ds <- cal.dataset(error = "S3", nobs = 500, seed = 1)
  fit <- cal.bayesian(ds, models = c("errors", "noerrors"), priors = "Weak")
  s1 <- cal.draws(fit, "errors")$sigma
  s2 <- cal.draws(fit, "noerrors")$sigma
  expect_true(quantile(s1, 0.025) < 0.0275 && quantile(s1, 0.975) > 0.0275)
  expect_gt(mean(s2), mean(s1))
  expect_s3_class(cal.compare(fit), "compare.loo")
})

test_that("The mixed model recovers a shift between materials", {
  skipStan()
  ds1 <- cal.dataset(error = "S1", nobs = 200, seed = 1)
  ds2 <- cal.dataset(error = "S1", nobs = 200, seed = 2)
  ds2$D47 <- ds2$D47 + 0.02
  ds2$Material <- 2
  fit <- cal.bayesian(rbind(ds1, ds2), models = "mixed", priors = "Weak")
  d <- cal.draws(fit)
  shift <- d$alpha[d$Material == "2"] - d$alpha[d$Material == "1"]
  expect_true(quantile(shift, 0.025) < 0.02 && quantile(shift, 0.975) > 0.02)
})

test_that("The mixed model is skipped with one material", {
  skipStan()
  ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
  expect_message(fit <- cal.bayesian(ds, models = c("noerrors", "mixed")), "mixed")
  expect_equal(names(fit), "BLM1_fit_NoErrors")
})

test_that("The ion model runs and gives draws for gamma", {
  skipStan()
  ds <- cal.dataset(error = "S1", nobs = 100, seed = 1)
  ds$Ion <- rnorm(nrow(ds))
  fit <- cal.ion.bayesian(ds, IonError = 0.1)
  expect_true("gamma" %in% names(cal.draws(fit)))
})
