#!/bin/bash -ve


../../Trinity --genome_guided_max_intron 1000 --genome_guided_bam transAligns.cSorted.bam --max_memory 2G  --output test_GG_use_small_multiscaff_bam_trinity_outdir
