#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################


docker run --rm -v `pwd`:`pwd` trinityrnaseq/trinityrnaseq \
              Trinity --seqType fq --max_memory 2G \
              --left `pwd`/reads.left.fq.gz \
              --right `pwd`/reads.right.fq.gz \
              --SS_lib_type RF \
              --CPU 4 \
              --output `pwd`/trinity_docker_outdir

