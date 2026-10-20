########################################################################
######### Description and example usage
########################################################################

## This script takes an SraRunTable from NCBI SRA Run Selector and 
## constructs the targets file required for downstream processing scripts 
## including subread_fastq_to_bam.R and subread_feature_count_gene.R.

## example usage: 
# R -f /scratch/yf0/sb0894/rnaseq/rnaseq_scripts/prepare_targets_file_from_sra_metadata.R --args -input=/scratch/yf0/sb0894/rnaseq/datasets/7387_2.0/sample_information/SraRunTable.txt -outdir=/scratch/yf0/sb0894/rnaseq/datasets/7387_2.0/sample_information/ -dsid=7387 -paired=TRUE


########################################################################
######## Source util scripts (relative to module root path) ############
########################################################################

source("scripts/inc/util.R")


########################################################################
######## Source arguments from command line ############################
########################################################################


## input argument must supply the full absolute path of input file 
## (sra metadata file). This file is assumed to be an SraRunTable.txt file 
## downloaded from SRA Run Selector and must contain the following columns: 
## GEO_Accession (exp), source_name, Run.

input.file <- ARGV.get("input", T)

## dsid: stemformatics dataset id

dsid <- ARGV.get("dsid", T)

## output directory: must be full absolute path

outdir <- ARGV.get("outdir", T)

## paired: is the data paired end? Logical (TRUE or FALSE)

paired <- ARGV.get("paired", T)

########################################################################
######### Read in metadata #############################################
########################################################################

sra.metadata <- read.table(input.file, sep=",", stringsAsFactors = FALSE, header=TRUE)



########################################################################
############# Extract and construct required columns ###################
########################################################################


if(paired == TRUE){
	
    # extract relevant columns    
    required.cols <- sra.metadata[ , c("GEO_Accession..exp.", "source_name", "Run", "Run")]
    
    # change column names to those required by downstream scripts
    names(required.cols) <- c("SampleID", "SampleType", "ReadFile1", "ReadFile2")
    
    # construct ReadFile1 column
    required.cols$ReadFile1 <- paste(required.cols$ReadFile1, "_1.fastq", sep="")
    
    # construct Readfile2 column
    required.cols$ReadFile2 <- paste(required.cols$ReadFile2, "_2.fastq", sep="")
    
    # update SampleID column to have dsid prepended
    required.cols$SampleID <- paste(dsid, "_", required.cols$SampleID, sep="")
    
    
} else if(paired !=TRUE){

    # extract relevant columns
    required.cols <- sra.metadata[,c("GEO_Accession..exp.", "source_name", "Run")]

    # change column names to those required by downstream scripts
    names(required.cols) <- c("SampleID", "SampleType", "ReadFile1")

    # construct ReadFile1 column
    required.cols$ReadFile1 <- paste(required.cols$ReadFile1, ".fastq", sep="")
    
    # update SampleID column to have dsid prepended
    required.cols$SampleID <- paste(dsid, "_", required.cols$SampleID, sep="")

}

## write targets file
write.table(required.cols, paste(outdir, dsid, "_targets.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
