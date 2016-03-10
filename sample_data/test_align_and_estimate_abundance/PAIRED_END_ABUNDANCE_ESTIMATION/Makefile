test: test_RSEM test_eXpress test_kallisto
	../../Analysis/DifferentialExpression/merge_matrices.pl   RSEM-trans.TPM.not_cross_norm eXpress-trans.TPM.not_cross_norm kallisto-trans.TPM.not_cross_norm > combined.TPM.not_cross_norm.matrix
	./plot_paired_comparisons.Rscript

test_RSEM:
	./align_and_estimate_tester_PAIRED_END.pl RSEM

test_eXpress:
ifeq ("", "$(shell which eXpress)")
	echo WARNING: skipping eXpress test, cant locate eXpress tool in PATH	
else
	./align_and_estimate_tester_PAIRED_END.pl eXpress
endif

test_kallisto:
ifeq ("", "$(shell which kallisto)")
	echo WARNING: cannot find kallisto installed in PATH, skipping kallisto test
else
	./align_and_estimate_tester_PAIRED_END.pl kallisto
endif


clean:
	./cleanme.pl


