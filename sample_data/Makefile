
DIRS = test_GenomeGuidedTrinity \
test_InSilicoReadNormalization \
test_DE_analysis \
test_align_and_estimate_abundance \
test_full_edgeR_pipeline \
test_GOSeq_trinotate_pipe \
test_profiling_report \
test_PtR



test_all: test_trin_assembly
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) test) || exit $$?; done

test_trin_assembly:
	cd test_Trinity_Assembly && make test_all


clean:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done
