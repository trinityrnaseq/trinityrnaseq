#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity.pl --seqType fq --JM 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 2 --no_cleanup

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi

