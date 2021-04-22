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
              --CPU 1 

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa trinity_out_dir/Trinity.fasta 90

    ./test_FL.sh --query trinity_out_dir/Trinity.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse
fi

