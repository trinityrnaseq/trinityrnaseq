#!/bin/bash

#source /broad/software/scripts/useuse

#reuse Perl-5.8
#reuse .samtools-0.1.19
#reuse GCC-4.9



CMD="`dirname $0`/run_Trinity_eval.pl $*"

eval $CMD

exit $?

