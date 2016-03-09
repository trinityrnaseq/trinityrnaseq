#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 4 --no_cleanup --long_reads longReads.fa --output test_trinity_long_reads

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi


