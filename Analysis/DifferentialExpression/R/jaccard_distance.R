jaccard_dist_matrix = function(x, by.row = F) {
    
    x = as.matrix(x);
    
    if (by.row == F) {
        x = t(x)
    }
    
    m = matrix(nrow=nrow(x), ncol=nrow(x));

    diag(m) = 0;
    
    colnames(m) = rownames(x);
    rownames(m) = rownames(x);
    
    # operates by rows.
    for (i in 1:(nrow(x)-1)) {

        message("row (",i,")");
                
        for (j in (i+1):nrow(x)) {
    
            num_same = sum(na.omit(x[i,] == x[j,]));
            
            num_diff = sum(na.omit(x[i,] != x[j,]));

            total = num_same + num_diff;
            jaccard = ifelse (total > 0, 1 - (num_same / (num_same + num_diff)), 1);
            m[i,j] = jaccard;
            m[j,i] = jaccard;
            
            
        }
    }
        
    return(m)
    
}
