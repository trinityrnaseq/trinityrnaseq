

test: test_edgeR test_DESeq2 test_voom test_ROTS

test_edgeR:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method edgeR -s samples.txt -o edgeR_outdir
	cd edgeR_outdir && ../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  --matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt

test_DESeq2:
ifeq "" "$(shell which DESeq2)"
	echo skipping DESeq2 test
else
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method DESeq2 -s samples.txt -o DESeq2_outdir
	cd DESeq2_outdir && ../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  --matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt 
endif

test_voom:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method voom -s samples.txt -o voom_outdir
	cd voom_outdir && ../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  --matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt 

test_ROTS:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method ROTS -s samples.txt -o ROTS_outdir --ROTS_B 50 # make it go faster, use more in practice
	cd ROTS_outdir && ../../../Analysis/DifferentialExpression/analyze_diff_expr.pl  --matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt 

clean:
	rm -rf ./edgeR_outdir
	rm -rf ./DESeq2_outdir
	rm -rf ./voom_outdir
	rm -rf ./ROTS_outdir
	rm -f ./Trinity_trans.counts.matrix.*
	rm -f ./Trinity_trans.TMM.EXPR.matrix.*
	rm -f ./*.pdf


