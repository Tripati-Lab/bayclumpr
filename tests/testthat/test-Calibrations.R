ds <- cal.dataset(error = "S1", nobs = 50)

test_that("OLS works", {
  a <- cal.ols(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})

test_that("weighted OLS works", {
  a <- cal.wols(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})

test_that("Deming works", {
  a <- cal.deming(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})

test_that("York works", {
  a <- cal.york(data = ds, replicates = 10)
  expect_equal(nrow(a), 10)
})
