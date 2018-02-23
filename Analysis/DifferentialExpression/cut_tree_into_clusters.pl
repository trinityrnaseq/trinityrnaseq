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
# --target <string>  'genes' or 'samples'
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
###################################################################################


__EOUSAGE__

    ;


my $Kmeans;
my $Ktree;
my $help_flag = 0;
my $R_data_file;
my $pct_height = 0;
my $target = undef;

my $lexically_order_columns;
my $no_column_reordering;

&GetOptions ( 'h' => \$help_flag,
              
              'target=s' => \$target,
              
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

unless ($target && $target =~ /^(genes|samples)$/ && ($Kmeans || $Ktree || $pct_height) && $R_data_file) {
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
    
    
    my $core_filename;
    my $outdir;
    
    if ($target eq 'samples') { 
        print $ofh "hc_features = hc_samples\n";
    }
    else {
        print $ofh "hc_features = hc_genes\n";
    }
    
    if ($Kmeans) {
        print $ofh "kmeans_clustering <- kmeans(data, centers=$Kmeans, iter.max=100, nstart=5)\n";
        $core_filename = "clusters_fixed_Kmeans_${Kmeans}";
        print $ofh "partition_assignments = kmeans_clustering\$cluster\n";
        
    }
    elsif ($Ktree) {
        print $ofh "partition_assignments <- cutree(as.hclust(hc_features), k=$Ktree)\n";
        $core_filename = "clusters_fixed_Ktree_${Ktree}";
    } 
    else {
        print $ofh "partition_assignments <- cutree(as.hclust(hc_features), h=$pct_height/100*max(hc_features\$height))\n";
        $core_filename = "clusters_fixed_P_${pct_height}";
        
    }
    
    print $ofh "partition_assignments = as.data.frame(partition_assignments)\n";
    print $ofh "cluster_names = paste('cl_', partition_assignments[,1], sep='')\n";
    print $ofh "partition_assignments = cbind(cluster_names, rownames(partition_assignments))\n";
    print $ofh "write.table(partition_assignments[order(cluster_names),], col.names=F, row.names=F, file=\'$core_filename.${target}_clusters\', quote=F, sep='\t')\n";
    
    close $ofh;
    

    &process_cmd("R --vanilla -q < $R_script");

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
