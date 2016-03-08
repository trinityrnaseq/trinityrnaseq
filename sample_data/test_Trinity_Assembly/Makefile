
test_trinity:
	./runMe.sh


test: test_assembly test_abundance_estimation test_bowtie_PE_read_estimates


test_assembly:
	./runMe.sh
	./misc_run_tests/__test_runMe_include_normalization.sh
	./misc_run_tests/__runMe_with_qual_trimming.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalization.sh
	./misc_run_tests/__test_runMe_with_jaccard_clip.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalize_libs_separately.sh

test_abundance_estimation: test_RSEM test_eXpress test_kallisto

test_RSEM:
	../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir RSEM_PE


test_eXpress:
	../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF --est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference --output_dir eXpress_PE

test_kallisto:
	ifeq "" "$(shell which kallisto)"
		echo skipping kallisto test
	else
		../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --est_method kallisto --trinity_mode --prep_reference --output_dir kallisto_PE
	endif

test_bowtie_PE_read_estimates:
	../../util/bowtie_PE_separate_then_join.pl  --seqType fq --left reads.left.fq --right reads.right.fq --target trinity_out_dir/Trinity.fasta --aligner bowtie
	../../util/SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.bam


test_full:
	./runMe.sh 1

clean:
	./cleanme.pl

