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


${TRINITY_HOME}/Trinity --seqType fq \
              --max_memory 2G \
              --left reads.left.fq,reads2.left.fq \
              --right reads.right.fq,reads2.right.fq \
              --SS_lib_type RF \
              --CPU 4 \
              --trimmomatic \
              --normalize_reads \
              --normalize_by_read_set \
              --output trinity_trim_and_norm_outdir



exit 0
