#!/bin/bash -ve

../../Trinity --samples_file samples.PE.txt \
              --seqType fq \
              --max_memory 1G \
              --output trinity_test_samples_PE
