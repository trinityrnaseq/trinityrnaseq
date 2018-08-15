#!/bin/bash -ve

${TRINITY_HOME}/Trinity --samples_file samples.PE.txt \
              --seqType fq \
              --max_memory 1G \
              --min_kmer_cov 3 \
              --output trinity_test_samples_PE_min_kmer_cov_3


exit 0

