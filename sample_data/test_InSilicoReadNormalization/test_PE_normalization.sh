#!/bin/bash -ve

# just for testing purposes, use --max_cov 30 or higher for real applications.
../../util/insilico_read_normalization.pl --JM 2G --left ../test_Trinity_Assembly/reads.left.fq.gz --right ../test_Trinity_Assembly/reads.right.fq.gz --seqType fq --max_cov 5 --pairs_together --no_cleanup --CPU 4 --tmp_dir_name tmp_PE_norm_dir



