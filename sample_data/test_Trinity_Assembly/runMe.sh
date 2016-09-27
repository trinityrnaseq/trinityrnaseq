#!/bin/bash -ve


if [ -e reads.right.fq.gz ] && [ ! -e reads.right.fq ]; then
    gunzip -c reads.right.fq.gz > reads.right.fq
fi

if [ -e reads.left.fq.gz ] && [ ! -e reads.left.fq ]; then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi

if [ -e reads2.right.fq.gz ] && [ ! -e reads2.right.fq ]; then
    gunzip -c reads2.right.fq.gz > reads2.right.fq
fi

if [ -e reads2.left.fq.gz ] && [ ! -e reads2.left.fq ]; then
    gunzip -c reads2.left.fq.gz > reads2.left.fq
fi



#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

../../Trinity --seqType fq --max_memory 2G \
              --left reads.left.fq.gz,reads2.left.fq.gz \
              --right reads.right.fq.gz,reads2.right.fq.gz \
              --SS_lib_type RF \
              --CPU 4 

##### Done Running Trinity #####

