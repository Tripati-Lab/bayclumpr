ds <- sim.clumped(error = "S1", nobs = 50)

test_that("OLS works", {
  a <- simulateLM_measured(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})

test_that("weighted OLS works", {
  a <- simulateLM_inverseweights(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})

test_that("Deming works", {
  a <- simulateDeming(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})

test_that("York works", {
  a <- simulateYork_measured(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})
