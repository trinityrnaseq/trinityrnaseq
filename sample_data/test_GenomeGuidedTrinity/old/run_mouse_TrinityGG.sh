#!/bin/bash -ve

if [ -e mm9chr17.fasta.gz ] && [ ! -e mm9chr17.fasta ]; then
    gunzip -c mm9chr17.fasta.gz > mm9chr17.fasta
fi

if [ -e mm9chr17.annotation.bed.gz ] && [ ! -e mm9chr17.annotation.bed ]; then
    gunzip -c mm9chr17.annotation.bed.gz > mm9chr17.annotation.bed
fi

if [ -e mm9chr17.tophat.bam ] && [ ! -e mm9chr17.tophat.sam ]; then
    samtools view mm9chr17.tophat.bam > mm9chr17.tophat.sam
fi


../../util/support_scripts/prep_rnaseq_alignments_for_genome_assisted_assembly.pl --coord_sorted_SAM mm9chr17.tophat.sam -I 100000 --SS_lib_type RF 

find Dir_mm9chr17.tophat.* -regex ".*reads" > read_files.list

../../util/support_scripts/GG_write_trinity_cmds.pl --reads_list_file read_files.list --paired --SS  > trinity_GG.cmds


../../trinity-plugins/parafly/bin/ParaFly -c trinity_GG.cmds -CPU 2 -vv

## execute the trinity commands, and then pull together the aggregate fasta file like so:
find Dir_mm9chr17.tophat.*  -name "*inity.fasta"  | ../../util/GG_trinity_accession_incrementer.pl > mm9.Trinity-GG.fasta
# which will ensure that each trinity accession is unique.

