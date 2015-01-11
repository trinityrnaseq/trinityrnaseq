
library("edgeR")

target_list_file_to_DGEList = function(targetsFile, minSumTagCounts = 5) {
	
	# format should be like so, including header:
#   files	group	description
#   fileA.txt  grpA    my grpA replicate
#   fileB.txt  grpB    my grpB replicate
#   ...

	targets = read.delim(file=targetsFile, stringsAsFactors = FALSE)
	d = readDGE(targets)

    # only use those rows that have a minimum number of tags (min to be DE)
    d = d[rowSums(d$counts) >= minSumTagCounts, ];
    
    # reset library sizes
    d$samples$lib.size = colSums(d$counts)
    
    # run TMM normalization
    d = calcNormFactors(d)
    

	return(d)
}


edgeR_DE_analysis = function (DGEList, dispersion=0.1, FDR=0.05) {

	de.tgw = exactTest(DGEList, dispersion=dispersion)

    topTagCount = -1;
    if (! is.null(de.tgw$table$p.value)) {  # the old edgeR way (R-2.12 and earlier)	
    	topTagCount = sum(de.tgw$table$p.value <= FDR);
    }
    else {
        topTagCount = sum(de.tgw$table$PValue <= FDR); # the latest way of doing it.
    }
    
	toptgw = topTags(de.tgw, n = topTagCount);
	
	detags_toptgw = toptgw$table[toptgw$table[,4] <= FDR,];
	
	detags_names = rownames(detags_toptgw)

	plotSmear(DGEList, de.tags=detags_names)

	return(detags_toptgw); # table with entries of interest.
}


	
