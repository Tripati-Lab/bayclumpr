ds <- cal.dataset(error = "S1", nobs = 50)

paramdist <- cal.ols(data = ds, replicates = 10)
recData <- data.frame(Sample = c("Sample 1", "Sample 2"),
                      D47 = c(0.6, 0.7),
                      D47error = c(0.005, 0.005),
                      N = c(2,2),
                      Material = c(1,1))


test_that("Number of samples is correct", {
  a <- rec.clumped(recData = recData, obCal = paramdist)
  expect_equal(nrow(a), 2)
})

test_that("Number of columns is correct", {
  a <- rec.clumped(recData = recData, obCal = paramdist)
  expect_equal(ncol(a), 5)
})
