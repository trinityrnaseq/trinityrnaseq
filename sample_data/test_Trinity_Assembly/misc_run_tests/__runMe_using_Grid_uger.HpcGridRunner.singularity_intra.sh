#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

cmd="../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 4 --grid_exec '/home/unix/bhaas/GITHUB/HpcGridRunner/hpc_cmds_GridRunner.pl --grid_conf /home/unix/bhaas/GITHUB/HpcGridRunner/hpc_conf/BroadInst_UGER.conf -c '  --output trinity_outdir_uger_simg_intra --singularity_img /seq/RNASEQ/TOOLS/TRINITY/SINGULARITY/trinityrnaseq.v2.15.0-predev.simg"


echo $cmd

eval $cmd



