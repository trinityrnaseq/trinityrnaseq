#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

${TRINITY_HOME}/Trinity --seqType fq --max_memory 2G \
              --left reads.left.fq.gz \
              --right reads.right.fq.gz \
              --SS_lib_type RF \
              --incl_salmon \
              --CPU 4 --include_supertranscripts --output trinity_incl_salmon

exit 0

