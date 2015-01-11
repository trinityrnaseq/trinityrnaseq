#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity.pl --seqType fq --JM 2G --left reads.left.fq.gz --right reads.right.fq.gz  --CPU 4 

##### Done Running Trinity #####

exit 0
