
plot_log2_fpkm_dist = function (fpkm_files) {

    num_files = length(fpkm_files);

	data_list = list();
	
	max_y = 0
	xlim = c(0,0)

	for (i in 1:num_files) {

		file = fpkm_files[i]
      	data = read.table(file, header=F)
		data = data[,6];
		data = log2(data+1)
		
		den = density(data);
	
		data_list[[i]] = den

		y = max(den$y);
		if (y > max_y) {
			max_y = y;
		}
	
		x = min(den$x)
		if (x < xlim[1]) {
			xlim[2] = x
		}	

		x = max(den$x);
		if (x > xlim[2]) {
			xlim[2] = x
		}


	}
	
	colors = rainbow(num_files);

	for (i in 1:num_files) {

		if (i == 1) {

			plot(data_list[[1]], col=colors[1], xlim=xlim, ylim=c(0,max_y), xlab="log2(fpkm+1)")
		}
		else {
			points(data_list[[i]], col=colors[i], type='l')
		}
	}

	return;
}





plot_expressed_gene_counts = function(fpkm_file, 
                                      title="expressed transcript counts vs. min fpkm", 
                                      fpkm_range=seq(0,10,0.2), 
                                      total=0, 
                                      outfile="count_summary.txt") {

    data = read.table(fpkm_file, header=T, row.names=1);
	data = data[,5]
	
	
	counts_expressed = c();
    counts_not_expressed = c();

    print_not_expressed_flag = 1;
    if (total == 0) {
        total = length(data);
        print_not_expressed_flag = 0;
    }

    count_expressed = total;
	for (i in fpkm_range) {
   
        if (i > 0) {
       		count_expressed = sum(data>=i)
		}
        count_not_expressed = total - count_expressed;
	
		counts_expressed[length(counts_expressed)+1] = count_expressed;
		counts_not_expressed[length(counts_not_expressed)+1] = count_not_expressed;
		
	}


    orig_settings = par(mfrow=c(1,2))

	plot(fpkm_range, counts_expressed, type='o', col='black', xlab="min(fpkm)", 
        main=title, ylab="count of transcripts", ylim=c(0,max(counts_expressed, counts_not_expressed)))
    
    if (print_not_expressed_flag) {
    	points(fpkm_range, counts_not_expressed, type='o', col='blue')
        legend('topright', c('expressed', 'not expressed'), col=c('black', 'blue'), pch=15);
    }   
	else {
		legend('bottomright', c('expressed transcripts'), col=c('black'), pch=15);
	}
   
	data_table = data.frame(fpkm_range=fpkm_range, expressed=counts_expressed, not_expressed=counts_not_expressed);
	
	write.table(data_table, file=outfile, quote=F, sep='\t', row.names=F);
	

    ## make density plot for non-zero FPKM values
    data = data[data>0]
    plot(density(log2(data)), xlab="log2(fpkm)")


    par(orig_settings)
    
}
