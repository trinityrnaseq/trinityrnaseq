#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 4 --grid_exec /home/unix/bhaas/utilities/trinity_uger_cmd_processor.py --output trinity_outdir_uger


