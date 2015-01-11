manually_define_clusters = function(hclust_obj, data_matrix) {

    plot(as.hclust(hclust_obj), labels=F, sub="", xlab="")
	
	myClusters <- identify(as.hclust(hclust_obj), N=200, MAXCLUSTER=length(hclust_obj$order))

	dirname = paste("manually_defined_clusters_", length(myClusters), sep='')
	dir.create(dirname)

	for (i in 1:length(myClusters)) {
		cluster_i = as.data.frame(myClusters[[i]])
		gene_names_i = rownames(cluster_i)
		data_i = data_matrix[gene_names_i,]
		
		outfile = paste(dirname, "/cluster_", i, sep='')
		write.table(data_i, file=outfile, sep="\t", quote=F)
	}

}

