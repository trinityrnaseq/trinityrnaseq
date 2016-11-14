#!/bin/sh -ve

# run the pipeline
../../util/run_Trinity_edgeR_pipeline.pl  --samples_file `pwd`/samples_n_reads_decribed.txt 

# get assembly stats
../../util/TrinityStats.pl trinity_out_dir/Trinity.fasta

# filter based on fpkm
../../util/filter_fasta_by_rsem_values.pl  --rsem_output=ds_rep1.isoforms.results,hs_rep1.isoforms.results,log_rep1.isoforms.results,plat_rep1.isoforms.results --fasta=trinity_out_dir/Trinity.fasta --fpkm_cutoff=0.1 --output=filtered.fasta

