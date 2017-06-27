test:
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --heatmap --heatmap_colorscheme 'black,yellow'
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --sample_cor_matrix
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --prin_comp 3 --add_prin_comp_heatmaps 5 --Zscale
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --compare_replicates
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --imgfmt svg --boxplot_log2_dist 0.8
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --top_variable_genes 300 --min_across_ANY_samples_gene_expr_val 5 --gene_cor_matrix
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --mean_vs_sd
	head -n11 candidate.Trinity_trans.TMM.EXPR.matrix > test.matrix && ../../Analysis/DifferentialExpression/PtR -m test.matrix -s candidate.samples.txt --per_gene_plots --per_gene_plots_incl_vioplot
	touch ./test


clean:
	rm -f ./candidate.Trinity_trans.TMM.EXPR.matrix.*
	rm -f ./test
	rm -f ./test.matrix ./test.matrix.*
	rm -f ./*.pdf
	rm -f ./*.svg
