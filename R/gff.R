
## GFF parsing code imported from legacy pbutils codebase
## (//depot/software/R/common/pbutils)

readGFF <- function(gffFile, attrClasses = c(), keepAttributes = FALSE,
                    sep = "\t", fill = TRUE, flush = TRUE) {
  ##
  ## Based on the description here:
  ## http://www.sanger.ac.uk/resources/software/gff/spec.html
  ##
  d <- read.table(gffFile, comment.char = '#', stringsAsFactors = FALSE,
                  sep = sep, fill = fill, flush = flush,
                  col.names=c('seqid','source','type','start','end','score','strand',
                              'frame','attributes'))

  if (nrow(d) == 0) { d }
  else {
    attrs <- strsplit(d$attributes, split = ";")
    keyValues <- lapply(attrs, function(a) strsplit(a, split = "="))
    attrNames <- unique(unlist(lapply(keyValues, function(a) sapply(a, "[[", 1))))
    flat <- lapply(keyValues, function(kv) {
      l <- vector("list", length(attrNames))
      l[1:length(l)] <- NA
      names(l) <- attrNames
      l[match(sapply(kv, "[[", 1), names(l))] <- sapply(kv, function(z) {
        if (length(z) >= 2) z[2] else ''
      })

      return(l)
    })
    w <- if(keepAttributes) {
      w <- 1:ncol(d)
    } else {
      -which(colnames(d) == "attributes")
    }
    nd <- cbind(d[,w], do.call(rbind, lapply(flat, unlist)),
                stringsAsFactors = FALSE)

    for (i in seq.int(attrClasses)) {
      class(nd[,names(attrClasses)[i]]) <- attrClasses[i]
    }
    return(nd)
  }
}

readVariantsGFF <- function(gffFile, keepAttributes = TRUE) {
  a <- readGFF(gffFile, attrClasses = c("coverage" = "integer",
                          "confidence" = "numeric",
                          "length" = "integer"), keepAttributes = keepAttributes)
  a <- a[, !(colnames(a) %in% c("source", "frame"))]
  a$x <- (a$start + a$end)/2
  a$length[a$type == "SNV"] <- 1
  a$genotype[a$type == "deletion"] <- '-'
  return(a)
}

readAlignmentSummaryGFF <- function(gffFile, keepAttributes = TRUE) {
  a <- readGFF(gffFile, attrClasses = c("ins" = "integer", "del" = "numeric",
                          "snv" = "integer"),
               keepAttributes = keepAttributes)
  a <- a[, !(colnames(a) %in% c("source", "frame"))]

  pCommaSep <- function(nm, colnames) {
    x <- as.data.frame(do.call(rbind, strsplit(a[,nm], ",")))
    x[,] <- as.numeric(as.matrix(x))
    colnames(x) <- colnames
    return(x)
  }
  a <- cbind(a, pCommaSep("cov2", c("v.mean", "v.sd")))
  a <- cbind(a, pCommaSep("cov", c("v.min", "v.med", "v.max")))
  a <- cbind(a, pCommaSep("gaps", c("v.ngaps", "v.gapsize")))
  a$x <- (a$start + a$end)/2
  return(a)
}
