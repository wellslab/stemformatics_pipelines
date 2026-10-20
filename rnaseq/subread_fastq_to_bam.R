## ==========================================================================
## File:    subread_fastq_to_bam
## Author:  O.Korn
## Created: 2014-09-03
## Modified: 2017-04-13
##
## Align reads (single end or paired end) to a given genome index.
## (The genome index must already exist).
##
## ==========================================================================
library("Rsubread")
library("limma")
library("tools")
## Include util script for args processing (relative to module root path)
source("scripts/inc/util.R")

## NOTE: All input paths expected to be full ABSOLUTE paths as we are
##   required to set our working directory to within the genome index path.
inputdir <- ARGV.get("fastqdir", T)
outdir <- ARGV.get("output", T)
indexname <- ARGV.get("indexname", T)
indexdir <- ARGV.get("indexdir", T)
targetsFile <- ARGV.get("targets", T)
ncpu <- ARGV.get("ncpu")
mismatches <- ARGV.get("mismatches")
indels <- ARGV.get("indels")
## Phred offset - 64 (old Illumina sequencers etc) or 33 (default)
## NOTE: Subread tells you when loaded phred scores don't look right. The default is
##       fine 99% of the time for fairly modern datasets.
## See: https://en.wikipedia.org/wiki/FASTQ_format
phredoffset <- ARGV.get("phredoffset")
## boolean to disable unique-only mappings
uniqueonly <- ARGV.get("uniqueonly")
## nBestLocations flag of align(). If not given, default = 1
## Maximum of 16 apparently? UPDATE: Can specify more but maybe ignored beyond 16?
nbestlocations <- ARGV.get("nbestlocations")
## nsubreads flag of align(). Default if not given is 10.
## Subread guide recommends nsubreads=4 for small RNA (e.g. for >=19bp miRNA)
## (Formula: nsubreads = read_length - 15)
nconsensus <- ARGV.get("nconsensus")
## Alignment type - DNA or RNA (default)
aligntype <- ARGV.get("type")
## Mininum fragment length for DNA alignment. Default 50
minfrag <- ARGV.get("minfrag")
## Maximum fragment length for DNA alignment. Default 600
maxfrag <- ARGV.get("maxfrag")


## Default nthreads=2 if not given
if (is.na(ncpu)) {
  ncpu <- 2
} else {
  ncpu <- abs(as.integer(ncpu))
}

## Default mismatches (SNPs) = 5 if not given
## NOTE: Subread default = 3
if (is.na(mismatches)) {
  mismatches <- 5
} else {
  mismatches <- abs(as.integer(mismatches))
}
## Default indels = 5 (Subread default) if not given
if (is.na(indels)) {
  indels <- 5
} else {
  indels <- abs(as.integer(indels))
}
## Default phred offset = 33
if (is.na(phredoffset)) {
  phredoffset <- 33
} else {
  phredoffset <- as.integer(phredoffset)
}
## Default unique mapping == TRUE
if (! is.na(uniqueonly)) {
  uniqueonly <- as.logical(uniqueonly)
} else {
  uniqueonly <- TRUE
}
## Max number of best locations to report for multimapping reads
if (! is.na(nbestlocations)) {
  nbestlocations <- as.numeric(nbestlocations)
} else {
  nbestlocations <- 1
}
## Min number of consensus subreads. Default is 10.
if (! is.na(nconsensus)) {
  nconsensus <- as.numeric(nconsensus)
} else {
  nconsensus <- 10
}
## If nconsensus less than 10, calculate "TH1" param based on 30% of value
## since TH1 could be set too high otherwise. Default TH1 value is 3 for nconsensus==10
## NOTE: We don't have to worry about 'TH2' (read pair 2) because default is 1.
TH1 <- 3
if (nconsensus < 10) {
  TH1 <- ceiling(nconsensus * 0.3)
  print(paste("NOTE: Got 'nconsensus' < 10, setting 'TH1' (minimum required hits out of 'nconsensus' reads for read pair 1) to ",TH1," accordingly",sep=""))
}
if (is.na(aligntype)) {
  aligntype <- "RNA"
} else {
  aligntype <- as.character(aligntype)
  if (! aligntype %in% c("RNA","rna","DNA","dna")) {
    stop("Error: Parameter '-type' value must be 'RNA' or 'DNA' (case insensitive)")
  }
}
if (is.na(minfrag)) {
  minfrag <- 50
} else {
  minfrag <- as.numeric(minfrag)
}
if (is.na(maxfrag)) {
  maxfrag <- 600
} else {
  maxfrag <- as.numeric(maxfrag)
}


## Targets file usage as per example:
## http://bioinf.wehi.edu.au/RNAseqCaseStudy/
##
## Must have columns: 'SampleID', 'ReadFile1', and if requiring paired-end, 'ReadFile2'
targets <- readTargets(targetsFile)

## Check for existence of column header "ReadFile1" (always required)
if (! length(grep("ReadFile1", colnames(targets), fixed=T))) {
  stop("Error: Targets file MUST specify 'ReadFile1' column")
}

## Subread expects index files to exist in current working directory.
## We don't need to worry about relative paths for 'inputdir', 'outdir' etc
## because this script expects full absolute paths.
oldwd <- setwd(indexdir)

## Align paired end reads (output BAM)
for (i in 1:nrow(targets)) {
  readfile1 <- paste(inputdir,"/",targets$ReadFile1[i],sep="")
  if (! file.exists(readfile1)) {
    stop(paste("Error: Targets file missing (or invalid) file name for read file 1 in line #", i, sep=""))
  }
  if ("ReadFile2" %in% colnames(targets) && ! is.na(targets$ReadFile2[i])) {
    readfile2 <- paste(inputdir,"/",targets$ReadFile2[i],sep="")
  } else {
    readfile2 <- NULL
  }
  input_format <- "gzFASTQ"
  ## Work out if input is gzipped FASTQ or not
  ext <- file_ext(readfile1)
  if (tolower(ext) != "gz") {
    input_format <- "FASTQ"
  }
  ## Align single-end or paired-end reads (if readfile2 specified) to target genome index and output aligned BAM
  ## NOTE: By default, best unique mapping per read retained only (no multi-mapping)
  align(index=indexname,
        readfile1=readfile1,
        readfile2=readfile2,
        input_format=input_format,
        output_format="BAM",
        output_file=paste(outdir,"/",targets$SampleID[i],".bam",sep=""),
        unique=uniqueonly,
        nBestLocations=nbestlocations,
        nsubreads=nconsensus,
        TH1=TH1,
        maxMismatches=mismatches,
        indels=indels,
        phredOffset=phredoffset,
        nthreads=ncpu,
        type=aligntype,
        minFragLength=minfrag,
        maxFragLength=maxfrag)
}
setwd(oldwd)
