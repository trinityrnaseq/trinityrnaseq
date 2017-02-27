
test:
	./runMe.sh
	touch test

test_all:
	./runMe.sh
	./misc_run_tests/__test_runMe_NO_normalization.sh
	./misc_run_tests/__runMe_with_qual_trimming.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalization.sh
	./misc_run_tests/__test_runMe_with_jaccard_clip.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalize_libs_separately.sh
	./misc_run_tests/__runMe_include_long_reads.sh
	./misc_run_tests/__run_PE_samples_file.sh
	./misc_run_tests/__run_SE_samples_file.sh
	./misc_run_tests/__runMe_use_workdir.sh
	touch test_all

clean:
	./cleanme.pl

