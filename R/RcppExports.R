# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Get Reference Fasta File Path From Dataset
#'
#' This function parses a referenceset.xml file and assuming it follows a "Standard"
#' format returns the path to the file name inside the file.
#'
#' @param filename The BAM file name
#' @export
#' @examples getFastaFileNameFromDatasetFile("referenceset.xml")
getFastaFileNameFromDatasetFile <- function(dataset_name) {
    .Call('pbbamr_getFastaFileNameFromDatasetFile', PACKAGE = 'pbbamr', dataset_name)
}

#' Get Aligned BAM Name
#'
#' This function parses a alignmentset.xml file and assuming it follows a "Standard"
#' format returns the path to the file name inside the file.
#'
#' @param filename The BAM file name
#' @export
#' @examples getBAMNameFromDatasetFile("alignmentset.xml")
getBAMNamesFromDatasetFile <- function(dataset_name) {
    .Call('pbbamr_getBAMNamesFromDatasetFile', PACKAGE = 'pbbamr', dataset_name)
}

#' Load BAM header
#'
#' This function loads the header of a BAM file into a List.  It will attempt to load the
#' list of sequences, program used and read groups as well as the version of the
#' BAM file, skipping over any information that is not present.
#'
#' @param filename The BAM file name
#' @export
#' @examples loadHeader("~git/pbbam/tests/data/dataset/bam_mapping_1.bam")
loadHeader <- function(filename) {
    .Call('pbbamr_loadHeader', PACKAGE = 'pbbamr', filename)
}

#' Load PBI BAM index file
#'
#' This function loads a pbi index file into a dataframe.  Depending on the
#' number of attributes present, it will either load just the basic data, or optionally
#' the mapping and barcode data.
#'
#' The original BAM file can also be read to gather additional covariates such as the SNR, read quality
#' and number of passes, though this may take longer.
#'
#' @param filename The BAM file name (without .pbi)
#' @param loadSNR Should we load the four channel SNR data? (Default = FALSE)
#' @param loadNumPasses Should we load the number of passes data? (Default = FALSE)
#' @param loadRQ Should we load the read quality? (Default = FALSE)
#' @param loadSC Load the SC tag for a scraps.bam file? (Only possible if file ends with '.scraps.bam')
#' @export
#' @examples loadPBI("~git/pbbam/tests/data/dataset/bam_mapping_1.bam")
loadPBI <- function(filename, loadSNR = FALSE, loadNumPasses = FALSE, loadRQ = FALSE, loadSC = FALSE) {
    .Call('pbbamr_loadPBI', PACKAGE = 'pbbamr', filename, loadSNR, loadNumPasses, loadRQ, loadSC)
}

#' Load BAM alignments as a list of data frames.
#'
#' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
#' @param bamName The BAM file name to grab
#' @param indexedFastaName The name of the indexed fasta file this should come from.
#'
#' @return Returns a list of alignments as data frames.  If the IPD and Pulse Width are available, they will be columns in the returned data as well.
#'
#' @export
loadDataAtOffsets <- function(offsets, bamName, indexedFastaName) {
    .Call('pbbamr_loadDataAtOffsets', PACKAGE = 'pbbamr', offsets, bamName, indexedFastaName)
}

#' Load BAM subreads as a list of data frames.
#'
#' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
#' @param bamName The BAM file name to grab
#'
#' @return Returns a list of subreads as a list of character strings.
#' @export
loadSubreadsAtOffsets <- function(offsets, bamName) {
    .Call('pbbamr_loadSubreadsAtOffsets', PACKAGE = 'pbbamr', offsets, bamName)
}

#' Load a section of the reference genome as a character string.
#'
#' @param id The name of the sequence
#' @param start The start of the sequence
#' @param end The end of the sequence
#' @param indexedFastaName The name of the indexed fasta file this should come from.
#'
#' @return Returns a character vector of the sequence
#' @export
loadReferenceWindow <- function(id, start, end, indexedFastaName) {
    .Call('pbbamr_loadReferenceWindow', PACKAGE = 'pbbamr', id, start, end, indexedFastaName)
}

#' Load BAM alignments as a list of list for the HMM model.
#'
#' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
#' @param bamName The BAM file name to grab
#' @param indexedFastaName The name of the indexed fasta file this should come from.
#' @param trimToLength How much should we subsample the alignments?
#'
#' @return Returns a list of phase2datasets as data frames.
#' @export
loadHMMfromBAM <- function(offsets, bamName, indexedFastaName, trimToLength = 140L) {
    .Call('pbbamr_loadHMMfromBAM', PACKAGE = 'pbbamr', offsets, bamName, indexedFastaName, trimToLength)
}

#' Load BAM alignment from a single ZMW as a list of list for the HMM model.
#'
#' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
#' @param bamName The BAM file name to grab
#' @param indexedFastaName The name of the indexed fasta file this should come from.
#' @param windowBreakSize We generate a new "window" every time 2 basepairs are matching in a particular gap.
#' @param minSize What is the minimum window size necessary for return, if the window at the end is less than this, we drop it. (NOT CURRENTLY USED).
#' @return Returns a list of phase2datasets as data frames.
#' @export
loadSingleZmwHMMfromBAM <- function(offsets, bamName, indexedFastaName, windowBreakSize = 140L, minSize = 50L) {
    .Call('pbbamr_loadSingleZmwHMMfromBAM', PACKAGE = 'pbbamr', offsets, bamName, indexedFastaName, windowBreakSize, minSize)
}

#' Load regions table
#'
#' Load a "regions table" for a subreads.bam (and optionally, a
#' corresponding `aligned_subreads.bam`).
#'
#' Presently, this needs to traverse the entire BAM files---the PBI
#' doesn't contain any "region type" information.  If this becomes
#' crucial we might ' want to consider including it in the PBI.
#'
#' @param subreadsBamName the .subreads.bam file name.  The corresponding .scraps.bam must be available in the same directory
#' @return a data frame with columns HoleNumber, RegionType, RegionStart, RegionEnd, indicating the base extent of regions of each type.
#' @export
loadRegionsTable <- function(subreadsBamName) {
    .Call('pbbamr_loadRegionsTable', PACKAGE = 'pbbamr', subreadsBamName)
}

#' Load some extra data values that are not in the pbi file, only in
#' the BAM.  This method is intended to be used to fetch information
#' about a significantly reduced subset of the entire BAM file
#' contents.
#'
#' This function returns a dataframe with same nrow dimension as df but augmented with extra columns.
#'
#' @param df the result of a call to loadPBI, or a subset of the rows of such a result
#' @param loadSNR Should we load the four channel SNR data? (Default = FALSE)
#' @param loadNumPasses Should we load the number of passes data? (Default = FALSE)
#' @param loadRQ Should we load the read quality? (Default = FALSE)
#' @param loadSC Load the SC tag for a scraps.bam file? (Only possible if file ends with '.scraps.bam')
#' @export
loadExtras <- function(df, loadSNR = FALSE, loadNumPasses = FALSE, loadRQ = FALSE, loadSC = FALSE) {
    .Call('pbbamr_loadExtras', PACKAGE = 'pbbamr', df, loadSNR, loadNumPasses, loadRQ, loadSC)
}

#' Load PBI BAM index file
#'
#' This function loads a pbi index file into a dataframe.  Depending on the
#' number of attributes present, it will either load just the basic data, or optionally
#' the mapping and barcode data.
#'
#' Input can be a BAM filename, or an XML dataset file---not a
#' FOFN.
#'
#' @export
loadPBI2 <- function(filename) {
    .Call('pbbamr_loadPBI2', PACKAGE = 'pbbamr', filename)
}

#' Get Per Read Metrics
#'
#' This function loads a dataset and parses each read, passing it to a class
#' which collects metrics on each read.  It returns a list of data frames, one
#' for each metric analyzed.
#'
#' @param datasetname The dataset/BAM file name.
#' @param indexedFastaName The fasta file used in the alignment.
#' @export
getReadReport <- function(datasetname, indexedFastaName) {
    .Call('pbbamr_getReadReport', PACKAGE = 'pbbamr', datasetname, indexedFastaName)
}

