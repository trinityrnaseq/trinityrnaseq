source("/seq/bhaas/SVN/trinityrnaseq/trunk/Analysis/DifferentialExpression/R/edgeR_funcs.R")
myDGEList = target_list_file_to_DGEList("__tmp_targets.1820.dat")
postscript("myc_vs_cysts.eps")
top_entries_table = edgeR_DE_analysis(myDGEList, dispersion=0.1, FDR=0.05); # set dispersion to zero for Poisson-equivalent results.
dev.off()
write.table(top_entries_table, quote=F, file="myc_vs_cysts.results.txt", sep="\t")
