test: test_RSEM test_eXpress test_kallisto
	../../Analysis/DifferentialExpression/merge_matrices.pl   RSEM-trans.TPM.not_cross_norm eXpress-trans.TPM.not_cross_norm kallisto-trans.TPM.not_cross_norm > combined.TPM.not_cross_norm.matrix
	./plot_paired_comparisons.Rscript

test_RSEM:
	./align_and_estimate_tester.pl RSEM

test_eXpress:
	./align_and_estimate_tester.pl eXpress

test_kallisto:
ifeq "" "$(shell which kallisto)"
	echo skipping kallisto test
else
	./align_and_estimate_tester.pl kallisto
endif


clean:
	./cleanme.pl

