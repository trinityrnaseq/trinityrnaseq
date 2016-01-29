	
source("../../../Analysis/DifferentialExpression/R/edgeR_funcs.R")
myDGEList = target_list_file_to_DGEList("targets.dat")
postscript("MA_plot.eps")
top_entries_table = edgeR_DE_analysis(myDGEList, dispersion=0.1, FDR=0.05); # set dispersion to zero for Poisson-equivalent results.
dev.off()
write.table(top_entries_table, quote=F, file="diff_expr.results.txt", sep="\t")
