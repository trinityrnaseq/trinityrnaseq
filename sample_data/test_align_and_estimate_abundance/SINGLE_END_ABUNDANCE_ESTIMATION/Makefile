TARGET="samples.txt"

test: test_RSEM test_eXpress test_kallisto test_salmon
	../../../Analysis/DifferentialExpression/merge_matrices.pl RSEM-trans.TPM.not_cross_norm express-trans.TPM.not_cross_norm kallisto-trans.TPM.not_cross_norm salmon-quasi-trans.TPM.not_cross_norm salmon-fmd-trans.TPM.not_cross_norm > combined.TPM.not_cross_norm.matrix
	../plot_paired_comparisons.Rscript

test_RSEM:
ifeq ("", "$(shell which rsem-calculate-expression)")
	echo WARNING: skipping RSEM test, cannot locate rsem-calculate-expression tool in PATH
else
	../align_and_estimate_tester.pl RSEM ${TARGET}
endif

test_eXpress:
ifeq ("", "$(shell which express)")
	echo WARNING: skipping express test, cant locate express tool in PATH	
else
	../align_and_estimate_tester.pl express ${TARGET}
endif

test_kallisto:
ifeq ("", "$(shell which kallisto)")
	echo WARNING: cannot find kallisto installed in PATH, skipping kallisto test
else
	../align_and_estimate_tester.pl kallisto ${TARGET}
endif

test_salmon:
ifeq ("", "$(shell which salmon)")
	echo WARNING: cannot find salmon installed in PATH, skipping salmon test
else
	../align_and_estimate_tester.pl salmon-quasi ${TARGET}
	../align_and_estimate_tester.pl salmon-fmd ${TARGET}
endif


clean:
	./cleanme.pl


