#!/bin/sh -ve

# run the pipeline
../../util/run_Trinity_edgeR_pipeline.pl  --samples_file `pwd`/samples_n_reads_decribed.txt 

# get assembly stats
../../util/TrinityStats.pl trinity_out_dir/Trinity.fasta

# filter based on fpkm
../../util/filter_fasta_by_rsem_values.pl  --rsem_output=ds_rep1.isoforms.results,hs_rep1.isoforms.results,log_rep1.isoforms.results,plat_rep1.isoforms.results --fasta=trinity_out_dir/Trinity.fasta --fpkm_cutoff=0.1 --output=filtered.fasta

# run bowtie alignment to examine read mapping stats:
../../util/bowtie_PE_separate_then_join.pl  --seqType fq --left reads.ALL.left.fq --right reads.ALL.right.fq --aligner bowtie --output read_content_analysis --target trinity_out_dir/Trinity.fasta -- -p 4 --all --best --strata -m 300

# examine read mapping stats:
../../util/SAM_nameSorted_to_uniq_count_stats.pl read_content_analysis/read_content_analysis.nameSorted.bam



