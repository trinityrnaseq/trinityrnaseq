#!/bin/bash -ve



#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz,reads2.left.fq.gz --right reads.right.fq.gz,reads2.right.fq.gz --SS_lib_type RF --CPU 4 --normalize_reads --no_cleanup --output trinity_test_no_qtrim

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi




