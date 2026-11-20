#!/bin/bash
# Count full-length reconstructed reference transcripts (FL_trans_analysis_pipeline.pl, default
# thresholds: >=99% of reference length, >=99% identity, <=1% gaps), with BLAT split over N chunks.
# usage: fl_eval.sh <TRINITY_HOME> <Trinity.fasta> <ref.fa> <out_dir> [nchunks=6]
# Requires blat on PATH. Use the SAME TRINITY_HOME (evaluator code) for every assembly being compared.
set -eu
TH=$1; Q=$(readlink -f $2); T=$(readlink -f $3); OD=$4; N=${5:-6}
export PATH=$PATH:$TH/trinity-plugins/slclust/bin
mkdir -p $OD/chunks && cd $OD
rm -f chunks/q.*
awk -v n=$N -v d=chunks 'BEGIN{i=-1} /^>/{i++} {print > (d "/q." (i%n) ".fa")}' $Q
for f in chunks/q.*.fa; do blat -t=dna -q=dna -out=pslx -noHead $T $f $f.pslx > /dev/null & done
wait
cat chunks/q.*.fa.pslx > FL.pslx
$TH/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $T --query $Q --out_prefix FL --reuse > FL.log 2>&1
$TH/util/TrinityStats.pl $Q > TrinityStats.txt
# FL.pslx.maps.summary columns: label, FL genes, FL transcripts, FL genes (multi-iso), FL transcripts (multi-iso)
cat FL.pslx.maps.summary
