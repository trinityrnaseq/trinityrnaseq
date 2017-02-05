test:
	tar xvf DE_results.tar.gz
	gunzip -c transcripts.TMM.fpkm.avg_reps.matrix.gz > transcripts.TMM.fpkm.avg_reps.matrix
	../../Analysis/DifferentialExpression/TissueEnrichment/DE_results_to_pairwise_summary.pl transcripts.TMM.fpkm.avg_reps.matrix . > DE_pairwise_summary.txt
	../../Analysis/DifferentialExpression/TissueEnrichment/pairwise_DE_summary_to_DE_classification.pl  DE_pairwise_summary.txt
	touch test

clean:
	rm -f ./*.DE_results
	rm -f ./DE_pairwise_summary.txt*
	rm -f ./transcripts.TMM.fpkm.avg_reps.matrix
	rm -f ./test
