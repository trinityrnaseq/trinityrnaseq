#!/bin/bash -ve


## include sense and antisense versions of the target transcripts
../../util/fasta_write_sense_n_anti.pl trinity_out_dir/Trinity.fasta > trinity_out_dir/Trinity.sense_n_anti.fasta


## align reads back to the transcripts

../../util/alignReads.pl --left reads.left.fq.gz --right reads.right.fq.gz --target trinity_out_dir/Trinity.sense_n_anti.fasta --aligner bowtie --seqType fq --SS_lib_type RF

# use RSEM to estimate read abundance
../../util/RSEM_util/run_RSEM.pl --transcripts trinity_out_dir/Trinity.sense_n_anti.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.sam.+.sam.PropMapPairsForRSEM.bam --paired


# compute FPKM values
../../util/RSEM_util/summarize_RSEM_fpkm.pl --transcripts trinity_out_dir/Trinity.sense_n_anti.fasta --RSEM RSEM.isoforms.results --fragment_length 300 --no_group_by_component | tee Trinity.sense_n_anti.RSEM.fpkm


