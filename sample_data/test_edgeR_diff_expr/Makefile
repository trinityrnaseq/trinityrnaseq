
DIRS = example_edgeR_analysis_pipeline testEdgeRfuncs 
#example_DESeq2_analysis_pipeline


test:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) test); done


clean:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) clean); done


