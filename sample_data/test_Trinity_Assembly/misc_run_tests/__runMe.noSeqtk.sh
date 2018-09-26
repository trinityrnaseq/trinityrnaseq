#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

if [ -z ${TRINITY_HOME} ]; then
    echo "Must set env var TRINITY_HOME"
    exit 1
fi


${TRINITY_HOME}/Trinity --seqType fq --max_memory 2G \
              --left reads.left.fq.gz \
              --right reads.right.fq.gz \
              --SS_lib_type RF \
              --CPU 4 \
              --output trinity_out_dir_noseqtk --NO_SEQTK

