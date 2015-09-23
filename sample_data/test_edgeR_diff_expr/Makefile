

test:
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method edgeR -s samples.txt -o edgeR_outdir
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method DESeq2 -s samples.txt -o DESeq2_outdir
	../../Analysis/DifferentialExpression/run_DE_analysis.pl  -m Trinity_trans.counts.matrix --method voom -s samples.txt -o voom_outdir

clean:
	rm -rf ./edgeR_outdir
	rm -rf ./DESeq2_outdir
	rm -rf ./voom_outdir



