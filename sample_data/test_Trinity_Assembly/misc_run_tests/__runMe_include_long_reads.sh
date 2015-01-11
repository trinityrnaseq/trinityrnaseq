#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity.pl --seqType fq --JM 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 4 --no_cleanup --long_reads longReads.fa

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi


