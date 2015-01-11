#!/bin/bash -ve

../../Trinity --genome_guided_bam SP2.chr.bam --max_memory 1G --CPU 2 --genome_guided_max_intron 1000 --jaccard_clip --SS_lib_type RF --output test_Schizo_trinityGG_jaccard_RF_outdir
