library(pbbamr)
library(dplyr)
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

test_that("loadRefName", {
  d = loadPBI(bfile)
  expect_equal(as.character(d$ref[1]), "lambda_NEB3011")
  b = loadPBI(bamdset)
  expect_equal(as.character(b$ref[1]), "ecoliK12_pbi_March2013")
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

# This will test that pluses are correctly removed from the alignment
test_that("Pulses excluded", {
  ibam = "goat.bam"
  ind = loadPBI(ibam)
  fastaname = "All4mers_circular_72x_l50256.fasta"
  data = loadAlnsFromIndex(ind, fastaname)[[1]]
  source("goat.values.R")
  expect_equal(length(pm), nchar(fulseq))
  readnogaps = data$read[data$read!="-"]
  expect_equal(length(readnogaps)+1, nchar(seqb))# Plus 1 to the read as the soft-clipping
  readAsStr <- paste(as.character(readnogaps), sep="", collapse="")
  alnSeq <- substr(seqb, 0, nchar(seqb)-1)
  expect_equal(readAsStr, alnSeq)
  fulseqVec <- strsplit(fulseq, "")[[1]]
  upperCaseVec <- grepl("[A-Z]",fulseqVec)
  expect_equal(sum(upperCaseVec), nchar(seqb))
  # Check pkmid value matches the sequence
  selectedPkmid <- data$pkmid[data$read!="-"]
  selectedPm <- pm[upperCaseVec]
  selectedPm = head(selectedPm, length(selectedPm)-1)
  # pbbam has "float" pkmid, and pbbamr uses "double", leaving a small
  # difference. Here we just check this small difference.
  pkmidDiff = abs(selectedPm/10 - selectedPkmid)
  expect_false(any(pkmidDiff > 1e-4))
  # Check sf value matches the sequence
  selectedSf <- data$sf[data$read!="-"]
  selectedBamsf <- sf[upperCaseVec]
  selectedBamsf = head(selectedBamsf, length(selectedBamsf)-1)
  expect_equal(selectedSf,selectedBamsf)
})

test_that("Unaligned Data Loads", {
  # Load an internal sample with name/read/ipd/pw/pkmid/sf
  ibam = system.file("extdata", "internalsample.bam", package="pbbamr")
  ind = loadPBI(ibam)
  data = loadReadsFromIndex(ind, 1:10)
  nms = colnames(data[[1]])
  expect_equal(nms, c("name", "read", "ipd", "pw", "pkmid", "sf", "snrA", "snrC", "snrG", "snrT" ))
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


test_that("loadPBI2 & loadExtras work as loadPBI did", {

    df <- (loadPBI("AlignmentSet/m54006_160504_020705.alignmentset.xml",
                  loadSNR=T, loadRQ=T, loadNumPasses=T)
           %>% select(-starts_with("bc")))
    attr(df, "bam.file") <- NULL

    df2 <- loadPBI2("AlignmentSet/m54006_160504_020705.alignmentset.xml")
    ex <- loadExtras(df2, loadSNR=T, loadRQ=T, loadNumPasses=T)
    df2 <- cbind(df2, ex)

    expect_equal(df, df2)
})
