

test: test_edgeR test_DESeq2 test_voom

test_edgeR:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method edgeR -s samples.txt -o edgeR_outdir

test_DESeq2:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method DESeq2 -s samples.txt -o DESeq2_outdir

test_voom:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method voom -s samples.txt -o voom_outdir

clean:
	rm -rf ./edgeR_outdir
	rm -rf ./DESeq2_outdir
	rm -rf ./voom_outdir



