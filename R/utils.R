### Collection of helper functions to make manipulating data easier

#' Helper function for verifying arguments
#' @export
chkClass <- function(var, classname, msg) {
  if (charmatch(classname, class(var), 0) == 0) {
    stop(msg)
  }
}

#' Helper function to combine loaded data frames.
#'
#' Takes a list of data frames, and names for those data frames, and combines
#' them using the data.table packages' rbindlist command and then adds a factor
#' variable (called `Condition`) with the name used to identify the data frame.
#' Useful for combining different datasets and analyzing them.
#'
#' @param dfs A list of data frames to combine
#' @param cond_names A character vector of strings that will be used to identify each
#' dataset.  If missing, the names of the list elements are used
#' @examples
#' # Using two arguments and a separate vector of names
#' ad = list(df1, df2)
#' nms = c("myFirst", "mySecond")
#' all = combineConditions(ad, nms)
#'
#' # Using one argument and the names of variables
#' ad = list(myFirst = df1, mySecond = df2)
#' all = combineConditions(ad)
#' @export
combineConditions <- function(dfs, cond_names = names(dfs)) {
  # Verify args
  if (length(dfs) != length(cond_names)) {
    stop("Names must be the same size as the number of data frames.")
  }
  chkClass(cond_names, "character", "Names for data frames must be a character frame.")
  chkClass(dfs, "list", "Data frames must be a list")
  ncols = sapply(dfs, ncol)
  if (length(unique(ncols)) != 1) {
    allcolnames = unique(unlist(sapply(dfs, colnames)))
    for (i in 1:length(dfs)) {
      Missing <- setdiff(allcolnames, names(dfs[[i]]))  # Find names of missing columns
      dfs[[i]][Missing] <- NA  # Add them, filled with 'NA's
      dfs
    }
  }
  if (any(colnames(dfs[[1]]) == "Condition")) {
    stop("Can't create a condition column as it already exists.")
  }
  n = length(dfs)
  Condition = factor(as.vector(
    unlist(sapply(1:n, function(i) rep(cond_names[i], nrow(dfs[[i]]))))
  ), levels = unique(cond_names))
  nd = data.table::rbindlist(dfs)
  nd$Condition = Condition
  nd
}

#' Function converts to Phred Scale, capping the error at a specified minimum
#'
#' @param errRate The error rate observed
#' @param minError The minimum allowed error rate.
#' @return A vector of phred scaled error rates
#' @examples
#' toPhred(10 ^ (-(1:4)), 1e-3)
#' @export
toPhred <- function(errRates, minError = -Inf) {
  clampedRates = sapply(errRates, max, minError)
  (-10 * log10(clampedRates))
}

#' Function converts from Phred to error rate
#'
#' @param qv The QVs associated with each value
#' @param minError The minimum allowed error rate.
#' @return A vector of error rates
#' @examples
#' fromPhred(seq(10, 70, 10), 1e-5)
#' @export
fromPhred <- function(qv, minError) {
  tmp = 10 ^ (-qv / 10)
  sapply(tmp, max, minError)
}

#' Convert a holeNumber into the X coordinate on the chip.
#'
#' Only applies to Sequel data where the hole number is
#' X << 16 & Y
#' @param holeNumber the ZMW hole number
#' @export
getHoleX <- function(holeNumber)  {
  holeNumber %/% 65536
}

#' Convert a hole number into the Y coordinate on the chip.
#'
#' Only applies to Sequel data where the hole number is
#' X << 16 & Y
#' @param holeNumber the ZMW hole number
getHoleY = function(holeNumber) {
  holeNumber %% 65536
}

#' Given a data frame returned by loadPBI and a fasta file, load a list of alignments
#'
#' @param index A data frame returned by loadPBI
#' @param indexed_fasta_name A character vector for the returned fasta file.
#' @param rows A list of rows to select the data from (default = use all data)
#' @seealso loadReadsFromIndex
#' @return A list of alignments as data frames
#' @export
loadAlnsFromIndex <- function(index, indexed_fasta_name, rows = NULL ) {
  chkClass(index, "data.frame", "Index is expected to be a data frame.")
  chkClass(indexed_fasta_name, "character", "Filename needs to be a character pointing to a Fasta file.")

  if(sum(colnames(index) %in% c("file", "offset")) != 2) {
    stop("Index was expected to have a file and offset column.")
  }

  # DO NOT SUBSET THE DATA.FRAME
  # R data.frames must have a row.names attribute, and under certain encodings
  # R will generate a new seq(1, nrow(df)) integer vector whenever this attribute
  # is asked for (see the `SEXP getAttrib(SEXP vec, SEXP name)` function)
  # For large data.frames this triggers an Allocation/GC/Fill routine that destroys
  # performance, so we avoid it by not subsetting the data.frame directly
  rng = if (is.null(rows)) 1:nrow(index) else rows

  lapply(rng, function(i) {
    bam_name = as.character(index$file[i])
    loadDataAtOffsets(index$offset[i], bamName = bam_name, indexedFastaName = indexed_fasta_name)[[1]]
  })
}

#' Get the full path to a reference FASTA from a dataset file.  For
#' XML reference sets, this is determined by the FASTA file referenced in
#' the XML; for FASTA reference sets, it is just the FASTA file itself.
##'
#' @param  p A referenceset XML file or reference FASTA
#' @export
getReferencePath <- function(p) {
  p = normalizePath(p)
  if (grepl(".fasta$", p)) { return(p) }
  dp = pbbamr::getFastaFileNameFromDatasetFile(p)
  split_path <- function(x) if (dirname(x) == x) x else c(basename(x),split_path(dirname(x)))
  ap = rev(split_path(dp))
  if (ap[1] == ".") {
    return(do.call(file.path, as.list(c(dirname(p), ap[2:length(ap)]))))
  } else {
    return(dp)
  }
}

#' Given a data frame returned by loadPBI  load a list of the unaligned sequencing data.
#' (this is the unaligned equivalent of loadAlnsFromIndex)
#'
#' @param index A data frame returned by loadPBI
#' @param rows A list of rows to select the data from (default = use all data)
#' @seealso loadAlnsFromIndex
#' @return A list of alignments as data frames
#' @export
loadReadsFromIndex <- function(index, rows = NULL ) {
  chkClass(index, "data.frame", "Index is expected to be a data frame.")

  if (sum(colnames(index) %in% c("file", "offset")) != 2) {
    stop("Index was expected to have a file and offset column.")
  }

  rng = if (is.null(rows)) 1:nrow(index) else rows

  lapply(rng, function(i) {
    bam_name = as.character(index$file[i])
    loadSubreadsAtOffsets(index$offset[i], bamName = bam_name)[[1]]
  })
}



