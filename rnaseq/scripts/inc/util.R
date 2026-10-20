## ===========================================================================
##                             ARGS PROCESSING
## ===========================================================================
## Captures args invoked like:
##   R -f command_args_processing.R --args -input="input.txt" -output="output.txt"
##
## Author: O.Korn

## Only capture the part given on the shell after the "--args" bit
BH_args <- commandArgs(trailingOnly=TRUE)

## Split into key/vals
BH_argskv <- strsplit(BH_args, "=")

## Create a dictionary lookup of arguments
## e.g. structure like:
##
## argsLookup [
##     input : input.txt
##     output : output.txt
##     ...
## ]
##
ARGV <- list()
for (kv in BH_argskv) {
  ## Strip any single surrounding double quotes from key=value string,
  ## strip the leading "-" from arg name while we're here, and
  ## any single or double quotes wrapping the value itself
  k <- gsub("^[\"]|[\"]$", "", kv[1])
  ARGV[gsub("^-", "", k)] <- gsub("^[\"\']+|[\"\']+$", "", kv[2])
}

ARGV.get <- function (arg, required=FALSE) {
  if (! arg %in% names(ARGV)) {
    if (required) {
      stop(paste("No such argument '",arg,"'!",sep=""))
    }
    warning(paste("No such argument '",arg,"'!",sep=""))
    NA
  } else {
    ARGV[[arg]]
  }
}
