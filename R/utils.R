### Collection of helper functions to make manipulating data easier

#' Helper function for verifying arguments
chkClass <- function(var, classname, msg) {
  if (class(var) != classname) {
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
    stop("All data frames must have the same number of columns.")
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
