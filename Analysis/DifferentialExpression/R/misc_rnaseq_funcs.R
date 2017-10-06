
plot_counts_matrix_log2_dist = function(matrix_file) {

	
	data = read.table(file=matrix_file, com='', row.names=1, header=T)

	conditions = colnames(data)
	colors = rainbow(length(conditions))


	plot(density(log2(data[,1])), col=colors[1], main=matrix_file, xlab='log2(frag_counts)', ylab='density')

	for (i in 2:length(data[1,])) {

		points(density(log2(data[,i])), type='l', col=colors[i])

	}

	legend('topright', conditions, col=colors, pch=15)

}


matrix_to_color_assignments = function(matrix_m, col=NULL, by=c("matrix", "row", "col")) {

	if (! is.matrix(matrix_m))
		stop("Error, matrix_to_color_assignments() requires a matrix as parameter.")
	num_colors = 0
    
    if (is.null(col)) {
        num_colors = min(nrow(matrix_m), ncol(matrix_m))
        col = rainbow(num_colors)
    }
    else {
        num_colors = length(col)
    }
    
    by = match.arg(by)
    
    if (by == "matrix") {

        min_val = min(matrix_m, na.rm=T)
	    matrix_m = matrix_m - min_val
	    max_val = max(matrix_m, na.rm=T)
	    matrix_m = matrix_m / max_val * num_colors
        #print(matrix_m)
   	    matrix_m = apply(matrix_m, 1:2, function(x) ifelse (x<1, as.character(col[1]), as.character(col[x])));
		
        matrix_m = matrix(as.character(matrix_m), nrow=dim(matrix_m)[1])
	}
	else {

		row_or_col_only_color_selector_func = function(x) { 
				a = min(x, na.rm=T); 
				b = max(x, na.rm=T); 
				c = (x-a)/(b-a) * num_colors;
                c = round(c);
				c = ifelse (c<1, 1, c); 
                #print(paste(c("color selection: (a)", a, " (b)", b, " (c)", paste(c, sep=',')))); 
                colors = as.character(col[c]);
                return(colors);
		}
	
		if (by == "row") {
            matrix_m = apply(matrix_m, 1, row_or_col_only_color_selector_func);
            print(matrix_m)
            print("dim matrix_m after apply"); print(dim(matrix_m))
            matrix_m = t(matrix_m);
            print("dim matrix_m after transpose: "); print(dim(matrix_m))
		}
		else {
			# by column
			matrix_m = apply(matrix_m, 2, row_or_col_only_color_selector_func);
		}
	}
    
	#print(matrix_m)
	return(matrix_m)
}

sample_matrix_to_color_assignments = function(sampleAnnotationsMatrix, colors) {

	if (missing(colors))
		colors = rainbow(nrow(sampleAnnotationsMatrix))

	nsamples = nrow(sampleAnnotationsMatrix);

	if (length(colors) < nrow(sampleAnnotationsMatrix))
		stop("Error, only ", length(colors), " colors specified, but have ", nsamples, " samples");

	for (i in 1:nrow(sampleAnnotationsMatrix)) {
		c = colors[i]
		sampleAnnotationsMatrix[i,] = sapply(sampleAnnotationsMatrix[i,], function(x) ifelse( x, as.character(c), 'white'))
	}

	return(sampleAnnotationsMatrix);

}

