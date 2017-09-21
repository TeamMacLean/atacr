Sys.setenv("R_TESTS" = "")
library(atacr)


context("normalisation functions")


test_that("library_size_normalisation_internal() returns proper values", {

  expected_counts <- matrix(c(2.5, 4, 4.375, 5, 5, 5, 7.5, 6, 5.625), nrow=3)
  expected_se <-  SummarizedExperiment::SummarizedExperiment(assays=list(counts=expected_counts))
  in_se <-  SummarizedExperiment::SummarizedExperiment(assays=list(counts=matrix(1:9, nrow=3)))
  out_se <- library_size_normalisation_internal(in_se)

  expect_equal(out_se, expected_se)

})

test_that("get_scaling_factors() gets proper values", {
  expected_factors <- c(2,1,0.666666667)
  tm <- matrix(c(rep(1,3), rep(2,3), rep(3,3)), nrow=3)
  expect_equal(get_scaling_factors(tm), expected_factors)
})


test_that("scale_normalise() returns proper values", {
  expected_mat <- matrix(rep(2,9), nrow=3)
  tm <- matrix(c(rep(1,3), rep(2,3), rep(3,3)), nrow=3)
  expect_equal(scale_normalise(tm, c(2,1,0.66667)), expected_mat, tolerance = 1e-05)
})
