#!/bin/bash -ve

../../Trinity --seqType fq --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF  --CPU 2 --jaccard_clip --max_memory 1G --output trinity_w_jaccard



