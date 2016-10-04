library(pbbamr)

# setwd("/Users/nigel/git/pbbamr/tests/testthat/")
# ffile = "/Users/nigel/pacbio/NewChemTests/ecoliK12_pbi_March2013.fasta"
ofile = "loadedAln.Rd"
refdset = "referenceset.xml"
bamdset = "AlignmentSet/m54006_160504_020705.alignmentset.xml"
sbamdset = "SubreadSet/m54006_160504_020705.tiny.subreadset.xml"
#v = pbbamr::getReadReport(bamdset, ffile)

test_that("misMatchReports", {
  # Simple combine with a name vector
  #v = pbbamr::getReadReport(bamdset, ffile)
  #a = data.frame(var1 = 1:3, var2 = 2:4)
})

test_that("ClippingReport", {
  ibam = system.file("extdata", "internalsample.bam", package = "pbbamr")
  fasta = system.file("extdata", "All4Mer.V2.11.fna", package = "pbbamr")
  data = getReadReport(ibam, fasta)
  clips = data[["clipping"]]
  # This just verifies that Clipped << Unclipped and that current numbers match
  # what they were before, I haven't verified exactly.
  expect_equal(clips[1, 2], 1229)
  expect_equal(clips[2, 2], 18568)
})
