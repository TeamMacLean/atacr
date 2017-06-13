library(atacr)
Sys.setenv("R_TESTS" = "")

context("summary and count functions")

test_that("target_count_summary() returns proper dataframe", {

  smry <- target_count_summary(sim_counts)
  expect_is(smry, "data.frame") #right class
  expect_vectors_equal(colnames(smry), c("sample","percent_on_target", "on_target", "off_target"))

})

test_that("coverage_count_summary() returns proper dataframe", {

  smry <- coverage_count_summary(sim_counts)
  expect_is(smry, "data.frame")
  expect_vectors_equal(colnames(smry), c("on_target", "off_target", "sample"))

})

test_that("target_count_coverage() returns proper dataframe", {
  cov <- target_count_coverage(sim_counts)
  expect_is(cov, "data.frame")
  expect_vectors_equal(colnames(cov), c("sample", "target", "count_sum", "mean_coverage"))
})

test_that("target_count_coverage() returns proper sized dataframe", {
  cov <- target_count_coverage(sim_counts)
  expect_length(cov$count_sum, 12)
  expect_equal(nrow(cov[cov$target == "on_target",]), 6)
  expect_equal(nrow(cov[cov$target == "off_target",]), 6)
})

test_that("target_count_coverage() returns proper values in dataframe", {
  cov <- target_count_coverage(sim_counts)
  expect_vectors_equal(cov$count_sum,c(15001.00,15170.00,14976.00,16665.77,16755.63,16640.31,355.00,359.00,360.00,364.00,405.00,376.00))
})

test_that("sample_kmeans() returns proper dataframe and proper values in the dataframe", {
  k <- sample_kmeans_cluster(sim_counts)
  expect_vectors_equal(colnames(k), c("cluster_id", "sample"))
  expect_vectors_equal(k$cluster_id, c(1,1,1,2,2,2))
})

test_that("count_windows_under_threshold() returns proper dataframe with proper values", {
  th <- count_windows_under_threshold(sim_counts, threshold=15)
  expect_vectors_equal(colnames(th), c("count", "threshold", "sample"))
  expect_vectors_equal(th$count, c(0,0,0,4,4,4))
})

test_that("calc_quantiles() returns list() when threshold == NULL", {
  l <- calc_quantiles(sim_counts)
  expect_is(l, "list")
  expect_has_all_and_only_these_members(l, c("bait_windows", "non_bait_windows"))
})

test_that("get_fit() gets fit", {
  f <- get_fit("norm", small_counts)
  expect_equal(f$chisqpvalue[1], 0.3154778 )
})

test_that("get_fits() gets fits", {
  f <- get_fits(small_counts)
  expect_equal(unique(f$distribution), c("norm", "pois", "nbinom"))
})

test_that("get_expected_values() returns right random numbers",{
  set.seed(1234)
  expect_vectors_equal(get_expected_values(c(1,2,3,4,1,2,3,4),dist="norm"),c(1.057280,2.831591,3.796155,-0.303645,3.012902,3.104852,1.813054,1.846650))
})

test_that("observed_expected_bins() gives right values", {
  l <- observed_expected_bins(c(1,2,3,4,1,2,3,4))
  expect_has_all_and_only_these_members(l, c("observed", "expected"))
  expect_vectors_equal(l$observed, c(8))
  expect_vectors_equal(l$expected, c(8))
})
