#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity \
    --seqType fq \
    --max_memory 2G \
    --left reads.left.fq.gz \
    --right reads.right.fq.gz \
    --SS_lib_type RF \
    --CPU 4 \
    --no_cleanup \
    --use_bowtie2 \
    --output test_trinity_bowtie2


##### Done Running Trinity #####


