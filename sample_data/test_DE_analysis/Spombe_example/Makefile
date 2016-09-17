test:
	../../../Analysis/DifferentialExpression/run_DE_analysis.pl -m Trinity_genes.counts.matrix --method edgeR -s samples.txt --dispersion 0.01 --output edgeR_genes

	cd edgeR_genes && ../../../../Analysis/DifferentialExpression/analyze_diff_expr.pl -m ../Trinity_genes.TMM.EXPR.matrix \
        -s ../samples.txt  --examine_GO_enrichment --GO_annots ../Trinotate_report.xls.gene_ontology --gene_lengths ../Trinity_genes.lengths && \
           ../../../../Analysis/DifferentialExpression/validate_UP_subset.Rscript ../samples.txt

clean:
	rm -rf ./edgeR_genes
