test:
	../../Trinity --samples_file samples_n_reads_decribed.txt --seqType fq --SS_lib_type RF --max_memory 1G --trimmomatic
	../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --samples_file samples_n_reads_decribed.txt --seqType fq
	find . -maxdepth 2  -name '*RSEM.genes.results' > rsem.genes.files
	../../util/abundance_estimates_to_matrix.pl  --est_method RSEM --name_sample_by_basedir --quant_files rsem.genes.files --out_prefix genes
	find . -maxdepth 2 -name "*RSEM.isoforms.results" > rsem.isoform.files
	../../util/abundance_estimates_to_matrix.pl  --est_method RSEM --name_sample_by_basedir --quant_files rsem.isoform.files --out_prefix isoforms
	../../Analysis/DifferentialExpression/run_DE_analysis.pl -m isoforms.counts.matrix -s samples_n_reads_decribed.txt --method edgeR -o edgeR_isoforms --dispersion 0.1
	../../Analysis/DifferentialExpression/run_DE_analysis.pl -m genes.counts.matrix -s samples_n_reads_decribed.txt --method edgeR -o edgeR_genes --dispersion 0.1
	../../util/misc/fasta_seq_length.pl trinity_out_dir/Trinity.fasta > Trinity.fasta.seq_lens
	../../util/misc/TPM_weighted_gene_length.py  --gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map  --trans_lengths Trinity.fasta.seq_lens --TPM_matrix isoforms.TMM.EXPR.matrix > Trinity.gene_lengths.txt

clean:
	./cleanme.pl

