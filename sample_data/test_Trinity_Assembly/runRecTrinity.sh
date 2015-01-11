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

../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz,reads2.left.fq.gz --right reads.right.fq.gz,reads2.right.fq.gz --SS_lib_type RF --CPU 4 --recursive_trinity --grid_node_max_memory 1G

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi


sleep 2


###########################################
# use RSEM to estimate read abundance  ####
###########################################

sleep 2

# first try RSEM
../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference

# try eXpress
../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference


###### Done running RSEM ########



