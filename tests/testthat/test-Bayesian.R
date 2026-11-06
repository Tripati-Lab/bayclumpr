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

test_that("Centering keeps the intercept on the usual scale", {
  skipStan()
  ds <- cal.dataset(error = "S1", nobs = 200, seed = 1)
  d <- cal.draws(cal.bayesian(ds, models = "noerrors", priors = "Weak"))
  expect_true(quantile(d$alpha, 0.025) < 0.268 && quantile(d$alpha, 0.975) > 0.268)
})

test_that("Both mixed models recover a shift between materials", {
  skipStan()
  ds1 <- cal.dataset(error = "S1", nobs = 200, seed = 1)
  ds2 <- cal.dataset(error = "S1", nobs = 200, seed = 2)
  ds2$D47 <- ds2$D47 + 0.02
  ds2$Material <- 2
  fit <- cal.bayesian(rbind(ds1, ds2), models = c("mixed_errors", "mixed_noerrors"), priors = "Weak")
  for (m in c("mixed_errors", "mixed_noerrors")) {
    d <- cal.draws(fit, m)
    shift <- d$alpha[d$Material == "2"] - d$alpha[d$Material == "1"]
    expect_true(quantile(shift, 0.025) < 0.02 && quantile(shift, 0.975) > 0.02)
  }

  # A new material has wider intervals than a known one
  recData <- data.frame(Sample = "A", D47 = 0.7, D47error = 0.01, Material = 1)
  known <- rec.bayesian(fit, recData, model = "mixed_errors")
  new <- rec.bayesian(fit, recData, model = "mixed_errors", newMaterial = TRUE)
  expect_gt(new$SE, known$SE)
})

test_that("Mixed models are skipped with one material", {
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

test_that("cal.summary reports ESS and Rhat for Bayesian fits", {
  skipStan()
  ds1 <- cal.dataset(error = "S1", nobs = 100, seed = 1)
  ds2 <- cal.dataset(error = "S1", nobs = 100, seed = 2)
  ds2$Material <- 2
  fit <- cal.bayesian(rbind(ds1, ds2), models = c("errors", "mixed_noerrors"), priors = "Weak")
  a <- cal.summary(fit, "errors")
  b <- cal.summary(fit, "mixed_noerrors")
  expect_true(all(a$ESS > 0) && all(a$Rhat > 0))
  expect_equal(nrow(b), 6)
  expect_false(anyNA(b$ESS))
})
