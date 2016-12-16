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

# monitoring at 1 second intervals because this test runs very quick.  You might monitor on the order of minutes rather than seconds for 'regular' runs.
../../Trinity --seqType fq --max_memory 1G --left reads.left.fq --right reads.right.fq --SS_lib_type RF --CPU 4 --monitoring --monitor_sec 1

