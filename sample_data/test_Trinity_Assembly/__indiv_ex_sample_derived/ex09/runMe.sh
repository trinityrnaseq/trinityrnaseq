#!/bin/sh
../../Trinity.pl --seqType fq --left ex9.reads.left.fq --right ex9.reads.right.fq --SS_lib_type RF --bfly_opts "--edge-thr=0.05 --stderr -V 18" --run_butterfly  --output trinity_outdir

