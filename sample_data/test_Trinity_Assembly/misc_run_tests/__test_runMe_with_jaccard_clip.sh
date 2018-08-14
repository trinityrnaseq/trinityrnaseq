#!/bin/bash -ve

${TRINITY_HOME}/Trinity --seqType fq --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF  --CPU 1 --jaccard_clip --max_memory 1G --output __test_trinity_w_jaccard --grid_node_max_memory 2G --grid_node_CPU 2


exit 0

