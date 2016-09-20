library(pbbamr)
#bfile = "/Users/nigel/git/pbbamr/tests/testthat//test.aligned.bam"
#ffile = "/Users/nigel/git/pbbamr/tests/testthat/lambdaNEB.fa"
#ofile = "/Users/nigel/git/pbbamr/tests/testthat/loadedAln.Rd"
#sbamdset = "/Users/nigel/git/pbbamr/tests/testthat/SubreadSet/m54006_160504_020705.tiny.subreadset.xml"
#setwd("/Users/nigel/git/pbbamr/tests/testthat/")
#setwd("/Users/ytian/Documents/Git/pbbamr/tests/testthat/")
bfile = "test.aligned.bam"
ffile = "lambdaNEB.fa"
ofile = "loadedAln.Rd"
refdset = "referenceset.xml"
bamdset = "AlignmentSet/m54006_160504_020705.alignmentset.xml"
sbamdset = "SubreadSet/m54006_160504_020705.tiny.subreadset.xml"

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

test_that("snrLoads", {
  d = loadPBI(bamdset, loadSNR = TRUE, loadRQ = TRUE)
  expect_equal(abs(4.73512 - d$snrA[1]) < 1e-5, TRUE)
})

test_that("headerMatches", {
  d = loadHeader(bfile)
  expect_equal(as.numeric(as.character(d$sequences$length[1])), 48502)
})


test_that("fastafound", {
  fasta = getFastaFileNameFromDatasetFile(refdset)
  expect_equal("./sequence/lambdaNEB.fasta", fasta)
})

test_that("SC Loads", {
  ibam = "SubreadSet/m54006_160504_020705.tiny.scraps.bam"
  ind = loadPBI(ibam, loadSC = TRUE)
  expect_equal(ind$sc[1], factor("LQRegion", levels = levels(ind$sc)))
})


test_that("Start Frame Loads", {
  ibam = "SubreadSet/m54006_160504_020705.tiny.scraps.bam"
  ind = loadPBI(ibam, loadSC = TRUE)
  expect_equal(ind$sc[1], factor("LQRegion", levels = levels(ind$sc)))
})

test_that("Unaligned Data Loads", {
  # Load an internal sample with name/read/ipd/pw/pkmid/sf
  ibam = system.file("extdata", "internalsample.bam", package="pbbamr")
  ind = loadPBI(ibam)
  data = loadReadsFromIndex(ind, 1:10)
})


test_that("alignedBAMfound", {
  bam = getBAMNamesFromDatasetFile(bamdset)
  expected = c("AlignmentSet/m54006_160504_020705.tiny_mapped.1.subreads.bam",
               "AlignmentSet/m54006_160504_020705.tiny_mapped.2.subreads.bam")
  expect_equal(expected, bam)
})

test_that("subreadsBAMfound", {
  bam = getBAMNamesFromDatasetFile("SubreadSet/m54006_160504_020705.tiny.subreadset.xml")
  expect_equal("SubreadSet/m54006_160504_020705.tiny.subreads.bam", bam[1])
})

test_that("FilteringWorks", {
  # Filter is RQ > 0.901
  index = loadPBI("SubreadSet/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreadset.xml", loadRQ = TRUE)
  index2 = loadPBI("SubreadSet/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0_filtered.subreadset.xml", loadRQ = TRUE)
  expect_true(min(index2$rq) > 0.901)
  expect_true(min(index$rq) < 0.901)

  })

test_that("IndexLoadAsInteger", {
  d = loadPBI(bfile)
  index <- c("tstart", "tend", "astart", "aend")
  for (i in index) {
    expect_equal(class(d[,i]), "integer")
  }
})
