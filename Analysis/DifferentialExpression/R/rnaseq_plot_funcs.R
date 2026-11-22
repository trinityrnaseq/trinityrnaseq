## MA and volcano plots for DE results, drawn with ggplot2 (#1153).
## Points are semi-transparent so dense regions stay readable; features with FDR <= 0.05 are red.
## The first top_gene_labels_show features are labeled (DE result tables are sorted by significance).
## Labels use ggrepel when it's installed, otherwise plain text labels.

suppressPackageStartupMessages(library(ggplot2))


.de_plot_df = function(featureNames, x, y, FDR, top_gene_labels_show) {

    FDR[is.na(FDR)] = 1  # eg. DESeq2 reports NA padj for filtered features

    df = data.frame(feature=featureNames, x=x, y=y, FDR=FDR, stringsAsFactors=FALSE)
    df$significant = factor(ifelse(df$FDR <= 0.05, "FDR <= 0.05", "not significant"),
                            levels=c("not significant", "FDR <= 0.05"))

    n_label = min(top_gene_labels_show, nrow(df))
    df$label = ""
    if (n_label > 0) {
        df$label[1:n_label] = df$feature[1:n_label]
    }

    return(df)
}


.de_plot_draw = function(df, xlab, ylab, title, pch, hline) {

    df = df[is.finite(df$x) & is.finite(df$y), , drop=FALSE]

    p = ggplot(df, aes(x=x, y=y, color=significant)) +
        geom_hline(yintercept=hline, linetype="dashed", color="grey50") +
        geom_point(alpha=0.4, size=1.2, shape=pch) +
        scale_color_manual(values=c("not significant"="black", "FDR <= 0.05"="red"), drop=FALSE) +
        labs(x=xlab, y=ylab, title=title, color=NULL) +
        theme_bw() +
        theme(legend.position="bottom")

    labeled = df[df$label != "", , drop=FALSE]
    if (nrow(labeled) > 0) {
        if (requireNamespace("ggrepel", quietly=TRUE)) {
            p = p + ggrepel::geom_text_repel(data=labeled, aes(label=label), size=2.5, color="black", max.overlaps=Inf)
        }
        else {
            p = p + geom_text(data=labeled, aes(label=label), size=2.5, color="black", vjust=-0.6, check_overlap=TRUE)
        }
    }

    print(p)
}


plot_MA = function(featureNames, logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20, top_gene_labels_show=20) {

    df = .de_plot_df(featureNames, logCounts, logFoldChange, FDR, top_gene_labels_show)
    .de_plot_draw(df, xlab, ylab, title, pch, hline=0)
}


plot_Volcano = function(featureNames, logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20, top_gene_labels_show=20) {

    df = .de_plot_df(featureNames, logFoldChange, FDR, FDR, top_gene_labels_show)

    # FDR of 0 would give an infinite -log10; cap those at the smallest nonzero FDR
    fdr = df$FDR
    min_nonzero = suppressWarnings(min(fdr[fdr > 0]))
    if (! is.finite(min_nonzero)) { min_nonzero = 1e-300 }
    fdr[fdr <= 0] = min_nonzero
    df$y = -1 * log10(fdr)

    .de_plot_draw(df, xlab, ylab, title, pch, hline=-1*log10(0.05))
}


plot_MA_and_Volcano = function(featureNames, logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {

    plot_MA(featureNames, logCounts, logFoldChange, FDR, xlab=xlab, ylab=ylab, title=title)
    plot_Volcano(featureNames, logFoldChange, FDR)

}
