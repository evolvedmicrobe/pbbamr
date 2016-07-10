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
#' @param names A character vector of strings that will be used to identify each
#' dataset.
#' @export
combineConditions <- function(dfs, names) {
  # Verify args
  if (length(dfs) != length(names)) {
    stop("Names must be the same size as the number of data frames.")
  }
  chkClass(names, "character", "Names for data frames must be a character frame.")
  chkClass(dfs, "list", "Data frames must be a list")
  ncols = sapply(dfs, ncol)
  if (length(unique(ncols)) != 1) {
    stop("All data frames must have the same number of columns.")
  }
  if (any(colnames(dfs[[1]]) == "Condition")) {
    stop("Can't create a condition column as it already exists.")
  }
  n = length(dfs)
  Condition = Condition = factor(as.vector(
    unlist(sapply(1:n, function(i) rep(names[i], nrow(dfs[[i]]))))
  ), levels = unique(names))
  nd = data.table::rbindlist(dfs)
  nd$Condition = Condition
  nd
}
