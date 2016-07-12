library(pbbamr)
#bfile = "/Users/nigel/git/pbbamr/tests/testthat//test.aligned.bam"
#ffile = "/Users/nigel/git/pbbamr/tests/testthat/lambdaNEB.fa"
#ofile = "/Users/nigel/git/pbbamr/tests/testthat/loadedAln.Rd"

bfile = "test.aligned.bam"
ffile = "lambdaNEB.fa"
ofile = "loadedAln.Rd"
refdset = "referenceset.xml"
bamdset = "alignmentset.xml"
sbamdset = "subreadset.xml"
test_that("dataMatches", {
  d = loadPBI(bfile)
  head(d)
  aln = loadDataAtOffsets(d$offset, bfile, ffile)
  #org_aln = aln
  #save(org_aln, file=ofile)
  load(ofile)
  for (i in 1:length(org_aln)) {
    expect_equal(any(aln[[i]] != org_aln[[i]]), FALSE)
  }
})

test_that("loadFromDataFrame", {
  d = loadPBI(bfile)
  aln = loadAlnsFromIndex(d, ffile)
  expect_equal(length(aln), 10)
  aln = loadAlnsFromIndex(d, ffile, 2:4)
  expect_equal(length(aln), 3)
})

test_that("headerMatches", {
  d = loadHeader(bfile)
  expect_equal(as.numeric(as.character(d$sequences$length[1])), 48502)
})


test_that("fastafound", {
  fasta = getFastaFileNameFromDatasetFile(refdset)
  expect_equal("./sequence/lambdaNEB.fasta", fasta)
})

test_that("alignedBAMfound", {
  bam = getBAMNameFromDatasetFile(bamdset)
  expect_equal("/home/nechols/src/PacBioSampleData/data/AlignmentSet/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.aligned_subreads.bam", bam)
})

test_that("subreadsBAMfound", {
  bam = getBAMNameFromDatasetFile(sbamdset)
  expect_equal("/home/nechols/src/PacBioSampleData/data/SubreadSet/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.bam", bam)
})
