#!/bin/sh
../../../Chrysalis/GraphFromFasta -i inchworm.K25.L25.fa -r both.fa -min_contig_length 100 -min_glue 0 -glue_factor 0 -min_iso_ratio 0 -t 4 -k 24 -kk 48
