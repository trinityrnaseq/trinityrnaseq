#!/bin/bash
# Assemble one benchmark data set.
# usage: run_trinity.sh <Trinity executable> <pombe|mouse> <data_dir> <out_dir> [CPU=6] [max_memory=20G]
# Run the same command for the previous release and the release candidate (e.g. inside each
# version's Docker/Apptainer image) so that the two assemblies differ only by Trinity version.
set -u
TRIN=$1; DS=$2; D=$3; OUT=$4; CPU=${5:-6}; MEM=${6:-20G}
mkdir -p $OUT && cd $OUT
start=$(date +%s)
$TRIN --version | head -1 > version.txt
$TRIN --seqType fq --left $D/reads/${DS}_10M.left.fq --right $D/reads/${DS}_10M.right.fq \
      --SS_lib_type RF --max_memory $MEM --CPU $CPU --output trinity_out_dir > run.log 2>&1
rc=$?
echo -e "$(cat version.txt)\t$DS\texit=$rc\telapsed_s=$(( $(date +%s) - start ))" | tee status.txt
exit $rc
