% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{loadSingleZmwHMMfromBAM}
\alias{loadSingleZmwHMMfromBAM}
\title{Load BAM alignment from a single ZMW as a list of list for the HMM model.}
\usage{
loadSingleZmwHMMfromBAM(offsets, bamName, indexedFastaName,
  windowBreakSize = 140L, minSize = 50L)
}
\arguments{
\item{offsets}{The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).}

\item{bamName}{The BAM file name to grab}

\item{indexedFastaName}{The name of the indexed fasta file this should come from.}

\item{windowBreakSize}{We generate a new "window" every time 2 basepairs are matching in a particular gap.}

\item{minSize}{What is the minimum window size necessary for return, if the window at the end is less than this, we drop it. (NOT CURRENTLY USED).}
}
\value{
Returns a list of phase2datasets as data frames.
}
\description{
Load BAM alignment from a single ZMW as a list of list for the HMM model.
}

