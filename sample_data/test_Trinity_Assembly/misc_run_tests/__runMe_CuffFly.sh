#!/bin/bash -ve


if [ -e reads.right.fq.gz ] && [ ! -e reads.right.fq ]; then
    gunzip -c reads.right.fq.gz > reads.right.fq
fi

if [ -e reads.left.fq.gz ] && [ ! -e reads.left.fq ]; then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi



#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity.pl --seqType fq --JM 2G --left reads.left.fq --right reads.right.fq --SS_lib_type RF --CPU 4 --CuffFly

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi


sleep 2


###########################################
# use RSEM to estimate read abundance  ####
###########################################

sleep 2

../../util/RSEM_util/run_RSEM_align_n_estimate.pl  --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF 

###### Done running RSEM ########



