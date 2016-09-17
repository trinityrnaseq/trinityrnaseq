

test: test_edgeR test_DESeq2 test_voom test_ROTS test_PtR_PCA

test_edgeR:
	../../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method edgeR -s samples.txt -o edgeR_outdir
	cd edgeR_outdir && ../../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  \
	--matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt && \
	../../../../Analysis/DifferentialExpression/validate_UP_subset.Rscript

test_DESeq2:
	../../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method DESeq2 -s samples.txt -o DESeq2_outdir
	cd DESeq2_outdir && ../../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  \
	--matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt && \
	../../../../Analysis/DifferentialExpression/validate_UP_subset.Rscript

test_voom:
	../../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method voom -s samples.txt -o voom_outdir
	cd voom_outdir && ../../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  \
	--matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt && \
	../../../../Analysis/DifferentialExpression/validate_UP_subset.Rscript

test_ROTS:

	../../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method ROTS -s samples.txt -o ROTS_outdir --ROTS_B 50 # make it go faster, use more in practice
	cd ROTS_outdir && ../../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  \
	--matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt && \
	../../../../Analysis/DifferentialExpression/validate_UP_subset.Rscript


test_PtR_PCA:
	../../../Analysis/DifferentialExpression/PtR -m Trinity_trans.TMM.EXPR.matrix -s samples.txt --prin_comp 3 --add_prin_comp_heatmaps 30

clean:
	rm -rf ./edgeR_outdir
	rm -rf ./DESeq2_outdir
	rm -rf ./voom_outdir
	rm -rf ./ROTS_outdir
	rm -f ./Trinity_trans.counts.matrix.*
	rm -f ./Trinity_trans.TMM.EXPR.matrix.*
	rm -f ./*.pdf



# in case we need to restrict tests here
#rversion := $(shell R --version | head -1)
#ifeq (,$(findstring 3.1.1, $(rversion)))
#  ...
#endif
