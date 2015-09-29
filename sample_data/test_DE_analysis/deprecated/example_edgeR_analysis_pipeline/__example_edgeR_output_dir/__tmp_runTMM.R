source("/seq/bhaas/SVN/trinityrnaseq/trunk/Analysis/DifferentialExpression/R/edgeR_funcs.R")
myDGEList = target_list_file_to_DGEList("all_data_files.list")
myDGEList$samples$eff.lib.size = myDGEList$samples$lib.size * myDGEList$samples$norm.factors
write.table(myDGEList$samples, file="TMM_info.txt", quote=F, sep="\t", row.names=F)
