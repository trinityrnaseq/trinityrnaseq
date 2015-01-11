#!/bin/bash -ve

if [ ! -e reads.single.fq ]; then
    gunzip -c ../test_Trinity_Assembly/reads.left.fq.gz > reads.single.fq
fi


# just for testing purposes, use --max_cov 30 or higher for real applications.
../../util/insilico_read_normalization.pl --JM 2G --single reads.single.fq --seqType fq --max_cov 5 --no_cleanup --tmp_dir_name single_tmp_norm_reads



