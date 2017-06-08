library(atacr)
Sys.setenv("R_TESTS" = "")

context("differential count functions")

test_that("get_t() returns proper value", {
  expect_equal( unname(get_t(1:100, c(1:10, 90:100))), -64.6472, tolerance = 0.0000001)
})


