library(pbbamr)
#bfile = "/Users/nigel/git/pbbamr/tests/testthat//test.aligned.bam"
#ffile = "/Users/nigel/git/pbbamr/tests/testthat/lambdaNEB.fa"
#ofile = "/Users/nigel/git/pbbamr/tests/testthat/loadedAln.Rd"

bfile = "test.aligned.bam"
ffile = "lambdaNEB.fa"
ofile = "loadedAln.Rd"

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

test_that("headerMatches", {
  d = loadHeader(bfile)
  expect_equal(as.numeric(as.character(d$sequences$length[1])), 48502)
})
