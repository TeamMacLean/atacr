Sys.setenv("R_TESTS" = "")
library(atacr)

context("methods")

test_that("as.data.frame.atacr() returns proper dataframe",{
  d <- as.data.frame(sim_counts)
  expect_vectors_equal(names(d), c("chromosome", "start", "stop", "strand", "sample", "count", "window_type"))
  expect_is(d$chromosome, "factor")
  expect_is(d$start, "integer")
  expect_is(d$stop, "integer")
  expect_is(d$sample, "factor")
  expect_is(d$count, "numeric")
  expect_is(d$window_type, "factor")
  expect_equal(levels(d$chromosome), c("synth_chrom"))

})
