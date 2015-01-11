#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################



../../Trinity --seqType fq --JM 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 4 --grid_conf_file ../../htc_conf/BroadInst_LSF.test.conf


