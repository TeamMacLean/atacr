library(atacr)

context("loading BAM files")

test_that("get_bait_regions_from_text() gets correct bait regions", {
  regions <- atacr::get_bait_regions_from_text('../../inst/extdata/individual_bait_regions.txt')
  expect_is(regions, "GRanges") #right class
  expect_that(levels(regions@seqnames), equals(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"))) #right seqnames
  expect_that(regions[1]@ranges@NAMES, equals("AT1G01680_1")) #right first name
  expect_that(regions[1]@ranges@start, equals(249021)) #right start
  expect_that(regions[1]@ranges@width, equals(120)) #right calculated width
  expect_that(length(regions), equals(2219)) #right number of regions
})

test_that("when loading bam files, region names load correctly", {

  data <- atacr::make_counts('../../inst/extdata/individual_bait_regions.txt', 'sample_treatment_bam_mappings_for_test.csv')
  expect_is(data$whole_genome, "RangedSummarizedExperiment")
  expect_is(data$bait_windows, "RangedSummarizedExperiment")
  expect_is(data$non_bait_windows, "RangedSummarizedExperiment")

  ## check that first range in each set of windows (genome, bait, none_baits) has right rownames - presumably is parsed correctly... this tests whether the windows are loaded correctly

  expect_that(names(data$whole_genome)[1], equals("Chr1:1-50"))
  expect_that(names(data$non_bait_windows)[1], equals("Chr1:1-50"))
  expect_that(names(data$bait_windows)[1], equals("Chr1:245951-246000"))



})
