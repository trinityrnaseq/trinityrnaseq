plot_MA = function(featureNames, logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20, top_gene_labels_show=20) {

    plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);;

    text(logCounts[1:top_gene_labels_show], logFoldChange[1:top_gene_labels_show], labels=featureNames[1:top_gene_labels_show], cex= 0.7, pos=3)
    
}


plot_Volcano = function(featureNames, logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20, top_gene_labels_show=20) {

    plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
    text(logFoldChange[1:top_gene_labels_show], (-1*log10(FDR))[1:top_gene_labels_show], labels=featureNames[1:top_gene_labels_show], cex= 0.7, pos=3)
    
}


plot_MA_and_Volcano = function(featureNames, logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {

    def.par = par(no.readonly = TRUE) # save default, for resetting...

    #gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    #layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 

    plot_MA(featureNames, logCounts, logFoldChange, FDR);
    plot_Volcano(featureNames, logFoldChange, FDR);

    # draw again, but use a smaller dot for data points
    #plot_MA(logCounts, logFoldChange, FDR, pch='.');
    #plot_Volcano(logFoldChange, FDR, pch='.');
    

    par(def.par)   
        
    
}
