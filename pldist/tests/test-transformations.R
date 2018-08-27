context("transformations")

test_that("transformation stops when it should", {
  
  expect_error( runningmean(0, c(0,0)) )
  
})
