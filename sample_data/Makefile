
DIRS = test_GenomeGuidedTrinity \
test_InSilicoReadNormalization \
test_Trinity_Assembly \
test_DE_analysis \
test_align_and_estimate_abundance \
test_full_edgeR_pipeline \
test_GOSeq_trinotate_pipe \
test_profiling_report



test:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) test) || exit $$?; done


clean:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done
