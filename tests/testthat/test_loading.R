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



all_atac <-   make_counts('individual_bait_regions.gff',
  'sample_treatment_bam_mappings_for_test.csv',
  filter_params = make_params(paired_map = FALSE, minq=1, dedup = F))

filtered_atac <-
  make_counts('individual_bait_regions.gff',
    'sample_treatment_bam_mappings_for_test.csv')

filtered_rnaseq <-
  make_counts('bait_genes.gff',
    'sample_treatment_bam_mappings_for_test.csv',
    is_rnaseq = TRUE
    )

all_rnaseq <-
  make_counts('bait_genes.gff',
    'sample_treatment_bam_mappings_for_test.csv',
    is_rnaseq = TRUE,
    filter_params = NULL
  )



test_that("when loading RNASeq, genome subsections are RangedSummarizedExperiments",
  {
    expect_is(all_rnaseq$whole_genome, "RangedSummarizedExperiment")
    expect_is(all_rnaseq$bait_windows, "RangedSummarizedExperiment")
    expect_is(all_rnaseq$non_bait_windows,
      "RangedSummarizedExperiment")

  })

test_that("when loading ATACSeq, genome subsections are RangedSummarizedExperiments",
  {
    expect_is(all_atac$whole_genome, "RangedSummarizedExperiment")
    expect_is(all_atac$bait_windows, "RangedSummarizedExperiment")
    expect_is(all_atac$non_bait_windows, "RangedSummarizedExperiment")

  })

test_that("when loading bam files for ATACSeq, region names load correctly", {
  # check that first range in each set of windows (genome, bait, none_baits) has right rownames - presumably is parsed correctly... this tests whether the windows are loaded correctly

  expect_that(names(all_atac$whole_genome)[1], equals("Chr1:1-50"))
  expect_that(names(all_atac$non_bait_windows)[1], equals("Chr1:1-50"))
  expect_that(names(all_atac$bait_windows)[1],
    equals("Chr1:245951-246000"))

})

test_that("when loading BAM files for RNAseq, region names are computed correctly",
  {
    expect_that(names(all_rnaseq$non_bait_windows)[1],
      equals("Chr1:1-246000"))
    expect_that(names(all_rnaseq$non_bait_windows)[2],
      equals("Chr1:246201-246700"))
    expect_that(names(all_rnaseq$bait_windows)[1], equals("FakeGeneA"))
    expect_that(names(all_rnaseq$bait_windows)[2], equals("FakeGeneB"))

  })

test_that("when loading BAM files for RNAseq, poor reads are filtered properly", {

  expect_that(unname(SummarizedExperiment::assay(filtered_rnaseq$bait_windows)["FakeGeneA",]), equals(c(0,4,3,70)))
  expect_that(unname(SummarizedExperiment::assay(filtered_rnaseq$bait_windows)["FakeGeneB",]), equals( c(3,7,7,186)))

})

test_that("when loading BAM files for ATACseq, poor reads are filtered properly", {

  expect_that(unname(SummarizedExperiment::assay(filtered_atac$bait_windows)["Chr1:246251-246300",]), equals(c(0,1,1,8)))

})

p <- make_csaw_params(make_params())

test_that("make_csaw_params() returns properly populated object", {
  expect_is(p, "readParam")
  expect_that(unname(p@pe), equals(c("both")))
  expect_that(p@max.frag, equals(500))
  expect_that(p@dedup, equals(TRUE))
  expect_that(p@minq, equals(30))
  expect_that(p@forward, equals(NA))
})

q <- make_scanBamParam(make_params(), filtered_rnaseq$bam_files[1])
test_that("make_scanBamParam() returns properly populated object", {

  expect_is(q, "ScanBamParam")
  expect_that(unname(q@flag), equals(c(2045,1023)))
  expect_that(q@simpleCigar, equals(FALSE))
  expect_that(q@reverseComplement, equals(FALSE))
  expect_that(q@mapqFilter, equals(30))
})
