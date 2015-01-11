#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use POSIX qw (floor ceil);

my $usage = <<__EOUSAGE__;


########################################################################
#
# --matrix <string>              matrix file
#
#
########################################################################
#
# Dimensions:
#
#  --width_per_plot <float>      default: 2.5
#  --height_per_plot <float>     default: 2.5
#
# Layout:
#
#  --plots_per_row <int>         default: 2
#  --plots_per_col <int>         default: 3
#
# Misc:
#
#  --barplot
#  --log2
#
########################################################################

__EOUSAGE__

    ;



my $width_per_plot = 2.5;
my $height_per_plot = 2.5;

my $matrix;

my $plots_per_row = 2;
my $plots_per_col = 3;
my $help_flag;

my $barplot_flag = 0;
my $log2_flag = 0;

&GetOptions( 'h' => \$help_flag,
             
             'matrix=s' => \$matrix,
             'width_per_plot=f' => \$width_per_plot,
             'height_per_plot=f' => \$height_per_plot,
             'plots_per_row=i' => \$plots_per_row,
             'plots_per_col=i' => \$plots_per_col,
             'barplot' => \$barplot_flag,
             'log2_flag' => \$log2_flag,
    );





if ($help_flag) {
    die $usage;
}

unless ($matrix) {
    die $usage;
}

main: {

    
    my $R_script = "__tmp_plot_clusters.R";

    open (my $ofh, ">$R_script") or die "Error, cannot write to $R_script";
   
    #print $ofh "postscript(file=\"my_cluster_plots.eps\", horizontal=FALSE, width=$width, height=$height, paper=\"special\")\n";
    print $ofh "pdf(file=\"$matrix.per_gene_plots.pdf\")\n";
    print $ofh "par(mfrow=c($plots_per_col, $plots_per_row))\n";
    #print $ofh "# png(file=\"my_cluster_plots.png\");\n"; 
    print $ofh "all_data = read.table(file=\"$matrix\", header=T, com=\'\', row.names=1)\n";
    print $ofh "gene_names = rownames(all_data)\n";

    if ($log2_flag) {
        print $ofh "all_data = log2(all_data+1);\n";
    }
    
    print $ofh "for (i in 1:length(all_data[,1])) {\n";
    print $ofh "    data = all_data[i,]\n";
    print $ofh "    ymin = min(data); ymax = max(data);\n";
    if ($barplot_flag) {
        print $ofh "    barplot(as.numeric(data), cex.names=0.75, names.arg=colnames(data), las=2, main=gene_names[i])\n";
    }
    else {
        print $ofh "    plot(as.numeric(data), type='l', ylim=c(ymin,ymax), main=gene_names[i], col='blue', xaxt='n', xlab='', ylab='')\n";
        print $ofh "    axis(side=1, at=1:length(data), labels=colnames(all_data), las=2)\n";
    }
    print $ofh "}\n";
    print $ofh "dev.off()\n";
    
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
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}
