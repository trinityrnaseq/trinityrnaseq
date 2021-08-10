#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

${TRINITY_HOME}/Trinity --seqType fq --max_memory 2G \
              --left reads.left.fq.gz \
              --right reads.right.fq.gz \
              --SS_lib_type RF \
              --no_symlink \
              --CPU 4 --output trinity_no_symlinks

exit 0

