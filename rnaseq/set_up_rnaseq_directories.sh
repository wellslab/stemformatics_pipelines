#!/bin/bash
# This script creates the directory structure required for RNAseq dataset processing

# DSHOME should be full path to directory where all datasets are to be stored
# DSID_VERSION should be in the format StemformaticsDatasetID_version. Eg. 1234_1.0

# example useage: ./set_up_rnaseq_directories.sh "~/Documents/rnaseq_datasets 1234_1.0

########################################################################
# begin
########################################################################
## obtain relevant arguments from command line

SCRIPT=$0

DSHOME=$1

DSID_VERSION=$2


echo -e "\n Checking if DSHOME directory exists \n"


# check if DSHOME exists. If not, create it.
if [ -d $DSHOME ]; then

    echo -e "DSHOME directory already exists. \n"

else

    echo -e "DSHOME directory does not exist. Creating directory DSHOME. \n"

    mkdir $DSHOME;

fi;



echo -e "Checking if DSID_VERSION directory exists. \n"

# check if DSID_VERSION diretory exists. If not, create it and its subdirectories.

if [ -d $DSHOME/$DSID_VERSION ]; then


    echo -e  "DSID_VERSION directory already exists. \n"

else 

    echo -e "Creating DSID_VERSION directory and subdirectories. \n"

    mkdir $DSHOME/$DSID_VERSION;
  
    cd $DSHOME/$DSID_VERSION
    
    mkdir bams
    
    mkdir feature_counts
    
    mkdir sample_information
   
    mkdir pbs_scripts
   
    mkdir FASTQ
   
    mkdir job_reports
   
    mkdir qc
    
    mkdir qc/picard_rna_metrics

    mkdir qc/fastqc_out

    mkdir qc/multiqc_data

    mkdir pcas;
    
fi;

echo "DONE"
