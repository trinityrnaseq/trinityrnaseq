#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

# stop before inchworm (just in silico norm)
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz  --CPU 4 --output trinity_piecemeal --no_run_inchworm

# stop before chrysalis
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz  --CPU 4 --output trinity_piecemeal --no_run_chrysalis

# stop before phase 2
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz  --CPU 4 --output trinity_piecemeal --no_distributed_trinity_exec

# finish it up
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz  --CPU 4 --output trinity_piecemeal 




##### Done Running Trinity #####

exit 0
