
get_Poisson_conf_intervals = function(seq_range, quantile_vec = c(0.05,0.95), plot=T) {
	
	num_rows = length(seq_range)
	num_cols = length(quantile_vec)
	m = matrix(nrow=num_rows, ncol=num_cols)
	
    colnames(m) = quantile_vec	
	rownames(m) = seq_range

	row_count = 0
	for (i in seq_range) {
		q = qpois(quantile_vec, i)
		row_count = row_count + 1
		m[row_count,] = q			
	
	}

    results = list()
    results$mat = m
	
	if (plot) {
		percents = plot_conf_intervals(m)
 	    
        results$pct = percents
    }

	return(results)

}


get_NB_conf_intervals = function(seq_range, dispersion, quantile_vec = c(0.05,0.95), plot=T) {
	
    size = 1/dispersion  #according to the mu-definition of size in nbinom of R, where var = mean + (1/size)mean^2
    
	num_rows = length(seq_range)
	num_cols = length(quantile_vec)
	m = matrix(nrow=num_rows, ncol=num_cols)
	
    colnames(m) = quantile_vec	
	rownames(m) = seq_range

	row_count = 0
	for (i in seq_range) {
		q = qnbinom(quantile_vec, mu=i, size=size)
		row_count = row_count + 1
		m[row_count,] = q			
	
	}
	
    results = list()
    results$mat = m
    

	if (plot) {
		percents = plot_conf_intervals(m)
 	
        results$pct = percents
    }

	return(results)

}






plot_conf_intervals = function(m) {

	par(mfrow=c(1,2))

	c_names = colnames(m)
	r_vals = as.numeric(rownames(m))

	max_val = max(m)
	
	plot(r_vals, r_vals, xlab="known read counts", ylab="Poisson read counts dist", t='l')

	# plot the confidence levels
	line_colors = rainbow(length(c_names))
	for(i in 1:length(c_names)) {
		points(r_vals,m[,i], t='l', col=line_colors[i])
	}

	# plot the max percentage of value for 95% conf level.
	percents=c()
	for (i in 1:length(r_vals)) {

		max_delta = max(abs(m[i,]-r_vals[i]))
		percent = max_delta/r_vals[i]*100
		
		percents[i]=percent;
	}
	
	plot(r_vals, percents, ylim=c(0,100), xlab="read counts", ylab="percent of value for 95% conf interval")

	percents

}
