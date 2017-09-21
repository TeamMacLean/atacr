library(atacr)
Sys.setenv("R_TESTS" = "")

context("differential count functions")

test_that("get_t() returns proper value", {
  expect_equal( unname(get_t(1:100, c(1:10, 90:100))), -64.6472, tolerance = 0.0000001)
})

test_that("select_comparisons() extracts proper columns", {
  l <- select_comparisons(sim_counts, "treatment", "control")
  expect_has_all_and_only_these_members(l, c("treatment_a_data", "treatment_b_data"))
  expect_vectors_equal(l$treatment_a_data, c("treatment_001", "treatment_002", "treatment_003"))
  expect_vectors_equal(l$treatment_b_data, c("control_001", "control_002", "control_003"))
})

test_that("estimate_fdr() returns proper dataframe", {
  expect_vectors_equal(names(estimate_fdr(sim_counts, "control", "treatment")), c("window", "t", "p_value", "fdr", "mean_count_a", "mean_count_b", "sd_a", "sd_b", "log2fc", "is_sig"))
})
