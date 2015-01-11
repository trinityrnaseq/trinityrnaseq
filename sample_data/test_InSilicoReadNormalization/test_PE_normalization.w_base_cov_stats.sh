#!/bin/bash -ve

if [ ! -e reads.left.fq ]; then
    gunzip -c ../test_Trinity_Assembly/reads.left.fq.gz > reads.left.fq
fi

if [ ! -e reads.right.fq ]; then
    gunzip -c ../test_Trinity_Assembly/reads.right.fq.gz > reads.right.fq
fi

# just for testing purposes, use --max_cov 30 or higher for real applications.
../../util/insilico_read_normalization.pl --JM 2G --left reads.left.fq --right reads.right.fq --seqType fq --max_cov 5 --pairs_together --__devel_report_kmer_cov_stats



