DIRS = test_Trinity_Assembly

test_all: test_trin_assembly

test_trin_assembly:
	cd test_Trinity_Assembly && make test_all


clean:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done
