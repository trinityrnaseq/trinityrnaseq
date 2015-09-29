#!/bin/bash -ve

../../Trinity --seqType fq --left top100k.Left.fq.gz --right top100k.Right.fq.gz --genome top100k.genome.gz --genome_guided_use_bam SP2.chr.bam --JM 1G --CPU 4 --genome_guided_max_intron 1000 --SS_lib_type RF --bfly_opts "--no_path_merging --triplet_strict -R 1 -O 10 --edge-thr=0.04" --genome_guided_sort_buffer 6G --genome_guided_CPU 2 --output test_Schizo_GG_outdir
# adding some extra params just to make sure all opts get forwarded on to where they're needed.

