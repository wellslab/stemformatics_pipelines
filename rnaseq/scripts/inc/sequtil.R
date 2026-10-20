## Author: O.Korn
##

## ===========================================================================
##                                  HELPERS
## ===========================================================================

## Convert CPM to RPKM (agnostic as to input frags or reads).
## Must supply gene / feature length sizes to calculate CPKM.
cpm2rpkm <- function (counts, isLogged=F, geneLengths=list()) {
  if (! is.matrix(counts) && ! is.data.frame(counts)) {
    write("Error: cpm2rpkm(): Input counts must be matrix or data frame!", stderr())
    return()
  }
  if (length(geneLengths) != nrow(counts)) {
    write("Error: cpm2rpkm(): Gene size list dimension mismatch to input count data!", stderr())
    return()
  }
  if (isLogged) {
    counts <- 2^counts
  }
  newCounts <- apply(base::cbind(counts, geneLengths), 1, function(x) {
    x[1:length(x)-1] / (tail(x,1) / 1000)
  })
  return(t(newCounts))
}
