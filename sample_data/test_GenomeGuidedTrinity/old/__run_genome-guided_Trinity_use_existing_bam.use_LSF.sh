#!/bin/bash -ve

if [ -e top100k.Left.fq.gz ] && ! [ -e top100k.Left.fq ]
then
	gunzip -c top100k.Left.fq.gz > top100k.Left.fq
fi

if [ -e top100k.Right.fq.gz ] && ! [ -e top100k.Right.fq ]
then
	gunzip -c top100k.Right.fq.gz > top100k.Right.fq
fi



if [ -e top100k.genome.gz ] && ! [ -e top100k.genome ]
then
	gunzip -c top100k.genome.gz > top100k.genome
fi



../../Trinity --left top100k.Left.fq --right top100k.Right.fq --jaccard_clip --genome top100k.genome --genome_guided_max_intron 1000 --genome_guided_use_bam SP2.chr.bam --JM 1G --seqType fq --grid_conf ../../htc_conf/BroadInst_LSF.test.conf --genome_guided_just_prep

../../Trinity --left top100k.Left.fq --right top100k.Right.fq --jaccard_clip --genome top100k.genome --genome_guided_max_intron 1000 --genome_guided_use_bam SP2.chr.bam --JM 1G --seqType fq --grid_conf ../../htc_conf/BroadInst_LSF.test.conf 

