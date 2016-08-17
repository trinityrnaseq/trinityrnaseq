
test_trinity:
	./runMe.sh


test: test_assembly test_bowtie_PE_read_estimates


test_assembly:
	./runMe.sh
	./misc_run_tests/__test_runMe_include_normalization.sh
	./misc_run_tests/__runMe_with_qual_trimming.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalization.sh
	./misc_run_tests/__test_runMe_with_jaccard_clip.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalize_libs_separately.sh
	./misc_run_tests/__runMe_include_long_reads.sh
	./misc_run_tests/__runMe_use_bowtie2.sh

test_bowtie_PE_read_estimates:
	../../util/bowtie_PE_separate_then_join.pl  --seqType fq --left reads.left.fq --right reads.right.fq --target trinity_out_dir/Trinity.fasta --aligner bowtie
	../../util/SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.bam


test_full:
	./runMe.sh 1

clean:
	./cleanme.pl

