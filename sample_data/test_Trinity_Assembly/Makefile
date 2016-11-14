
test_trinity:
	./runMe.sh

test_full: test

test:
	./runMe.sh
	./misc_run_tests/__test_runMe_NO_normalization.sh
	./misc_run_tests/__runMe_with_qual_trimming.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalization.sh
	./misc_run_tests/__test_runMe_with_jaccard_clip.sh
	./misc_run_tests/__runMe_with_qual_trimming_and_normalize_libs_separately.sh
	./misc_run_tests/__runMe_include_long_reads.sh



clean:
	./cleanme.pl

