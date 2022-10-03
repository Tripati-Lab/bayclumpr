
test_that("The dataset is generated", {
  ds <- cal.dataset(error = "S1", nobs = 50)
  expect_equal(nrow(ds), 50)
})


test_that("The minimum set of colums are present", {
  ds <- cal.dataset(error = "S1", nobs = 50)
  targetCols <- c("Temperature", "TempError", "D47error", "D47", "Material")
  expect_true(all( targetCols %in% colnames(ds) ))
})
