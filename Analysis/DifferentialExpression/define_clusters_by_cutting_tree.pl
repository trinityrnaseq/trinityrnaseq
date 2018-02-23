#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use File::Basename;
use FindBin;

my $usage = <<__EOUSAGE__;

###################################################################################
#
# -K <int>          define K clusters via k-means algorithm
#
#  or, cut the hierarchical tree:
#
# --Ktree <int>     cut tree into K clusters
#
# --Ptree <float>   cut tree based on this percent of max(height) of tree 
#
# -R <string>  the filename for the store RData (file.all.RData)
#
#  misc:
#
#   --lexical_column_ordering       reorder column names according to lexical ordering
#   --no_column_reordering  
#
###################################################################################


__EOUSAGE__

    ;


my $Kmeans;
my $Ktree;
my $help_flag = 0;
my $R_data_file;
my $pct_height = 0;

my $lexically_order_columns;
my $no_column_reordering;

&GetOptions ( 'h' => \$help_flag,
              'K=i' => \$Kmeans,
              'Ktree=i' => \$Ktree,
              'Ptree=f' => \$pct_height,
              'R=s' => \$R_data_file,
              
              'lexical_column_ordering' => \$lexically_order_columns,
              
              'no_column_reordering' => \$no_column_reordering,
              
              );




if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, don't understand args: @ARGV";
}

unless (($Kmeans || $Ktree || $pct_height) && $R_data_file) {
    die $usage;
}

if ($pct_height && $pct_height < 1) {
    die "Error, specify --Ptree as percent value > 1\n\n";
}

main: {
    
    unless (-s $R_data_file) {
        die "Error, cannot find pre-existing R-session data as file: $R_data_file";
    }
    
    
    my $R_script = "__tmp_define_clusters.R";
    
    open (my $ofh, ">$R_script") or die "Error, cannot write to file $R_script";

    print $ofh "library(cluster)\n";
    #print $ofh "library(gplots)\n";
    print $ofh "library(Biobase)\n";
    print $ofh "library(fastcluster)\n";
    print $ofh "source(\"$FindBin::RealBin/R/heatmap.3.R\")\n";
    
    print $ofh "load(\"$R_data_file\")\n";
    
    print $ofh "data = heatmap_data\n";
    
    my $core_filename;
    my $outdir;
    
    if ($Kmeans) {
        print $ofh "kmeans_clustering <- kmeans(data, centers=$Kmeans, iter.max=100, nstart=5)\n";
        $core_filename = "clusters_fixed_Kmeans_${Kmeans}.heatmap";
        $outdir = basename($R_data_file) . ".clusters_fixed_Kmeans_" . $Kmeans;
        print $ofh "gene_partition_assignments = kmeans_clustering\$cluster\n";
        
    }
    elsif ($Ktree) {
        print $ofh "gene_partition_assignments <- cutree(as.hclust(hc_genes), k=$Ktree)\n";
        $core_filename = "clusters_fixed_Ktree_${Ktree}.heatmap";
        $outdir = basename($R_data_file) . ".clusters_fixed_Ktree_" . $Ktree;
        
    } 
    else {
        print $ofh "gene_partition_assignments <- cutree(as.hclust(hc_genes), h=$pct_height/100*max(hc_genes\$height))\n";
        $core_filename = "clusters_fixed_P_${pct_height}.heatmap";
        $outdir = basename($R_data_file) . ".clusters_fixed_P_" . $pct_height;
    }
    
    # write gene order in heatmap clustering
    print $ofh "write.table(gene_partition_assignments[hc_genes\$order], file=\"$core_filename.heatmap_gene_order.txt\", quote=F, sep='\t')\n";
    
    print $ofh "max_cluster_count = max(gene_partition_assignments)\n";
    
    print $ofh "outdir = \"" . $outdir . "\"\n";
    print $ofh "dir.create(outdir)\n";
    
    # make another heatmap:
    print $ofh "partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)\n";
    print $ofh "gene_colors_dframe = data.frame(clusters=gene_partition_assignments, colors=partition_colors[gene_partition_assignments])\n";
    print $ofh "write.table(gene_colors_dframe, file=\"$core_filename.gene_cluster_colors.dat\", quote=F, sep='\t')\n";
    
    print $ofh "gene_colors = as.matrix(partition_colors[gene_partition_assignments])\n";
    print $ofh "pdf(\"$core_filename.heatmap.pdf\")\n";


    if ($lexically_order_columns) {
        print $ofh "data = data[,order(colnames(data))]\n";
    }

    if ($lexically_order_columns || $no_column_reordering) {
        print $ofh "heatmap.3(data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=F, col=myheatcol, RowSideColors=gene_colors, scale=\"none\", density.info=\"none\", trace=\"none\", key=TRUE, cexCol=1, margins=c(10,10))\n";
    }
    else {
        print $ofh "heatmap.3(data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale=\"none\", density.info=\"none\", trace=\"none\", key=TRUE, cexCol=1, margins=c(10,10))\n";
    }
    print $ofh "dev.off()\n";
    
    
    print $ofh "gene_names = rownames(data)\n";
    print $ofh "num_cols = length(data[1,])\n";
    
    
    print $ofh "for (i in 1:max_cluster_count) {\n";
    print $ofh "    partition_i = (gene_partition_assignments == i)\n";
    
    print $ofh "    partition_data = data[partition_i,,drop=F]\n";
    
#    print $ofh "    # if the partition involves only one row, then it returns a vector instead of a table\n";
#        ;
#    print $ofh "    if (sum(partition_i) == 1) {\n";
#    print $ofh "          dim(partition_data) = c(1,num_cols)\n";
#    print $ofh "          colnames(partition_data) = colnames(data)\n";
#    print $ofh "          rownames(partition_data) = gene_names[partition_i]\n";
#    print $ofh "    }\n";
    
    if ($lexically_order_columns) {
        print $ofh "partition_data = partition_data[,order(colnames(partition_data)), drop=F]\n";
    }
    
    elsif (! $no_column_reordering) {
        ## order based on sample clustering
        print $ofh "partition_data = partition_data[,hc_samples\$order, drop=F]\n";
    }
    
    print $ofh "    outfile = paste(outdir, \"/subcluster_\", i, \"_log2_medianCentered_fpkm.matrix\", sep='')\n";
    print $ofh "    write.table(partition_data, file=outfile, quote=F, sep=\"\\t\")\n";
    print $ofh "}\n";
    

    close $ofh;


    &process_cmd("R --vanilla -q < $R_script");

    
    ###################################################
    ## Generate the expression plots for each cluster
    ###################################################

    chdir $outdir or die "Error, cannot cd into $outdir";
    
    my $cmd = "$FindBin::RealBin/plot_expression_patterns.pl subcluster\*fpkm.matrix";
    &process_cmd($cmd);
    
    


    exit(0);
    

}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd $cmd died with ret $ret";
    }

    return;
}
