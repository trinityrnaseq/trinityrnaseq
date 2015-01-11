#!/bin/bash -ve

if [ ! -e top100k.genome ]; then
    gunzip -c top100k.genome.gz > top100k.genome
fi


../../Trinity --seqType fq --single top100k.Left.fq.gz --genome top100k.genome --JM 1G --CPU 2 --genome_guided_max_intron 1000 --SS_lib_type R --genome_guided_sort_buffer 2G --genome_guided_CPU 2 --output test_Schizo_GG_SE_outdir


