test:
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --heatmap --heatmap_colorscheme 'black,yellow'
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --sample_cor_matrix
	../../Analysis/DifferentialExpression/PtR -m candidate.Trinity_trans.TMM.EXPR.matrix -s candidate.samples.txt --log2 --prin_comp 3 --add_prin_comp_heatmaps 5 --Zscale

clean:
	rm -f ./candidate.Trinity_trans.TMM.EXPR.matrix.*

