library(atacr)
Sys.setenv("R_TESTS" = "")

context("loading BAM files")

test_that("get_bait_regions_from_text() gets correct bait regions", {
  regions <- get_bait_regions_from_text('individual_bait_regions.txt')

  expect_is(regions, "GRanges") #right class
  expect_that(levels(regions@seqnames), equals(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"))) #right seqnames
  expect_that(regions[1]@ranges@NAMES, equals("AT1G01680_1")) #right first name
  expect_that(regions[1]@ranges@start, equals(249021)) #right start
  expect_that(regions[1]@ranges@width, equals(120)) #right calculated width
  expect_that(length(regions), equals(2219)) #right number of regions

})

atac_data <-
  make_counts('individual_bait_regions.txt',
    'sample_treatment_bam_mappings_for_test.csv')

rnaseq_data <-
  make_counts('bait_genes.gff',
    'sample_treatment_bam_mappings_for_test.csv',
    is_rnaseq = TRUE)

test_that("when loading RNASeq, genome subsections are RangedSummarizedExperiments",
  {
    expect_is(rnaseq_data$whole_genome, "RangedSummarizedExperiment")
    expect_is(rnaseq_data$bait_windows, "RangedSummarizedExperiment")
    expect_is(rnaseq_data$non_bait_windows,
      "RangedSummarizedExperiment")

  })

test_that("when loading ATACSeq, genome subsections are RangedSummarizedExperiments",
  {
    expect_is(atac_data$whole_genome, "RangedSummarizedExperiment")
    expect_is(atac_data$bait_windows, "RangedSummarizedExperiment")
    expect_is(atac_data$non_bait_windows, "RangedSummarizedExperiment")

  })

test_that("when loading bam files for ATACSeq, region names load correctly", {
  # check that first range in each set of windows (genome, bait, none_baits) has right rownames - presumably is parsed correctly... this tests whether the windows are loaded correctly

  expect_that(names(atac_data$whole_genome)[1], equals("Chr1:1-50"))
  expect_that(names(atac_data$non_bait_windows)[1], equals("Chr1:1-50"))
  expect_that(names(atac_data$bait_windows)[1],
    equals("Chr1:245951-246000"))

})

test_that("when loading BAM files for RNAseq, region names are computed correctly",
  {
    expect_that(names(rnaseq_data$non_bait_windows)[1],
      equals("Chr1:1-246000"))
    expect_that(names(rnaseq_data$non_bait_windows)[2],
      equals("Chr1:246201-246700"))
    expect_that(names(rnaseq_data$bait_windows)[1], equals("FakeGeneA"))
    expect_that(names(rnaseq_data$bait_windows)[2], equals("FakeGeneB"))

  })
