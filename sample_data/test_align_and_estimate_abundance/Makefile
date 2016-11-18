
DIRS = PAIRED_END_ABUNDANCE_ESTIMATION \
       SINGLE_END_ABUNDANCE_ESTIMATION \
       PAIRED_END_ABUNDANCE_ESTIMATION_via_samples_file_direct \
       SINGLE_END_ABUNDANCE_ESTIMATION_via_samples_file_direct

test:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) test) || exit $$?; done


clean:
	@for i in $(DIRS); do \
	echo "Running clean in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done

