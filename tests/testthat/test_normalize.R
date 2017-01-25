library(atacr)
context("Raw to Normalized Data")

test_that("normalisation gives expected values", {
  df <- atacr::quick_normalize(tenbaits_threesamples, "S1")
  expect_equal(df$log_normalized_value[1], 0)
})

test_that("control_ratio_intensity_plot returns ggplot", {
  df <- atacr::quick_normalize(tenbaits_threesamples, "S1")
  expect_is(atacr::control_ratio_intensity_plot(df), "ggplot")
})
