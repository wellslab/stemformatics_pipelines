## ==========================================================================
## File: subread_feature_count_gene.R
##
## Annotate read counts to user supplied GTF feature annotation file.
##
## If paired end reads, must supply "-ispaired=TRUE" arg.
##
## NOTE: Annotates to highest level meta-feature ID in GTF file
##
## O.Korn
## 2014-09-05
## modified by S.Butcher 2021-03-02
## ==========================================================================


## example usage

## R -f /scratch/yf0/sb0894/7387_test/subread_feature_count_gene.R --args -inputdir=/scratch/yf0/sb0894/7387_test/output/ -outdir=/scratch/yf0/sb0894/7387_test/output/feature_counts/ -gtf=/scratch/yf0/sb0894/Homo_sapiens.GRCh38.91.gtf -ncpu=8 -ispaired=TRUE -targets=/scratch/yf0/sb0894/7387_test/FASTQ_Files/7387_targets_file_all_samples.txt -multimap=FALSE -ignoredup=TRUE -dsid=7387

## Note: This script must be run from the directory containing "scripts/inc/util.R" and "scripts/inc/sequtil.R". Required arguments include inputdir, outdir, gtf, dsid. Rsubread and edgeR packages are also required. 


## load libraries
library("Rsubread")
library("edgeR")

## Include util scripts (relative to module root path)
source("scripts/inc/util.R")
source("scripts/inc/sequtil.R")

## samdir or bamdir: Must be full absolute path
samBamDir <- ARGV.get("inputdir", T)

## Path to write summarise count tables
outdir <- ARGV.get("outdir", T)

## gtf: Must be full absolute path
gtf <- ARGV.get("gtf", T)

## Dataset ID
dsid <- ARGV.get("dsid", T)

## Number of threads, default = 2
ncpu <- ARGV.get("ncpu")

## Are input reads paired end? Default = FALSE
ispaired <- ARGV.get("ispaired")

## Limit input BAMs to those in targets file, if supplied
targetsfile <- ARGV.get("targets")

## Count multi-mapping reads?
do_multimap <- ARGV.get("multimap")

## Ignore duplicate reads? Default is TRUE
do_ignoredup <- ARGV.get("ignoredup")


if (is.na(dsid)) {
  stop("Error: dsid must be given!")
}

if (is.na(outdir)) {
  stop("Error: Output path must be given!")
}
if (is.na(samBamDir)) {
  stop("Error: Must supply directory containing SAM or BAM files!")
}

## Default nthreads=2 if not given
if (is.na(ncpu)) {
  ncpu <- 2
} else {
  ncpu <- as.integer(ncpu)
}

## Settings for counting multi mapping and fractional counting
do_fractional <- FALSE

if (is.na(do_multimap)) {
  do_multimap <- FALSE
} else {
  do_multimap <- as.logical(do_multimap)
  ## fractional counting enabled if multimap == TRUE
  ## NOTE: Output raw counts may not be whole numbers
  if (do_multimap) {
    do_fractional <- TRUE
  }
}
print(paste("DEBUG: do_multimap ==",do_multimap))

if (is.na(do_ignoredup)) {
  do_ignoredup <- TRUE
} else {
  do_ignoredup <- as.logical(do_ignoredup)
}
print(paste("DEBUG: do_ignoredup ==",do_ignoredup))


## Set SAM or BAM input files
if (! is.na(samBamDir)) {
  inputfiles <- list.files(samBamDir, pattern=".[bs]am$", ignore.case=T)
  inputnames <- basename(inputfiles)
} else {
  stop("Error: No SAM or BAM files found in input directory!")
}

## Input reads are paired end? Default = FALSE
if (! is.na(ispaired)) {
  ispaired <- as.logical(ispaired)
} else {
  ispaired <- F
}

oldwd <- setwd(samBamDir)


## Load targets file if supplied
targets <- NULL
if ( ! is.na(targetsfile)) {
  ## Targets file usage as per example:
  ## http://bioinf.wehi.edu.au/RNAseqCaseStudy/
  ##
  ## Must have columns: 'SampleID', 'ReadFile1', and if requiring paired-end, 'ReadFile2'
  targets <- readTargets(targetsfile)

  ## Check for existence of column header "ReadFile1" (always required)
  if (! length(grep("ReadFile1", colnames(targets), fixed=T))) {
    stop("Error: Targets file MUST specify 'ReadFile1' column")
  }
  bamFileTargets <- paste(targets$SampleID, ".bam", sep="")

  ## Override input files with those that intersect with sample IDs given in targets
  inputfiles <- intersect(inputfiles, bamFileTargets)
  inputnames <- gsub(".bam", "", basename(inputfiles), fixed=T)
}


## Transcript/exon level summarisation
featureCounts_gene <- featureCounts(
  inputfiles,
  annot.ext=gtf,
  isGTFAnnotationFile=T,
  useMetaFeatures=T,
  isPairedEnd=ispaired,
  allowMultiOverlap=T,
  reportReads=NULL,
  countMultiMappingReads=do_multimap,
  fraction=do_fractional,
  minMQS=30,
  nthreads=ncpu,
  ignoreDup=do_ignoredup)

## Reset sample names in counts table e.g. if were numeric, might have had 'X' prepended to sample names by featureCounts()
featureCounts_gene$targets <- inputnames
colnames(featureCounts_gene$counts) <- inputnames
save(featureCounts_gene, file=paste(outdir,"featureCounts_gene.RData",sep="/"))


## Write raw and transformed counts
write.table(round(featureCounts_gene$counts), paste(outdir, "/", dsid, ".raw.tsv",sep=""), sep="\t", quote=F, na="")


setwd(oldwd)
