#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use POSIX qw (floor ceil);

my $usage = <<__EOUSAGE__;


########################################################################

 usage: $0 [opts] cluster_1 cluster_2 ...

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
#  --plots_per_col <int>         default: 2
#
########################################################################

__EOUSAGE__

    ;



my $width_per_plot = 2.5;
my $height_per_plot = 2.5;

my $plots_per_row = 2;
my $plots_per_col = 2;
my $help_flag;

&GetOptions( 'h' => \$help_flag,
             'width_per_plot=f' => \$width_per_plot,
             'height_per_plot=f' => \$height_per_plot,
             'plots_per_row=i' => \$plots_per_row,
             'plots_per_col=i' => \$plots_per_col,
             );



my @cluster_files = @ARGV;

if ($help_flag) {
    die $usage;
}

unless (@cluster_files) {
    die $usage;
}

main: {

    ## ensure each cluster file can be found as a file
    foreach my $file (@cluster_files) {
        unless (-s $file) {
            die "Error, cannot find file \"$file\" ";
        }
                
    }


    my $R_script = "__tmp_plot_clusters.R";

    open (my $ofh, ">$R_script") or die "Error, cannot write to $R_script";
    print $ofh "files = c(\"" . join("\",\"", @cluster_files) . "\")\n";

    #print $ofh "postscript(file=\"my_cluster_plots.eps\", horizontal=FALSE, width=$width, height=$height, paper=\"special\")\n";
    print $ofh "pdf(file=\"my_cluster_plots.pdf\")\n";

    print $ofh "par(mfrow=c($plots_per_col, $plots_per_row))\n";
    print $ofh "par(cex=0.6)\n";
    print $ofh "par(mar=c(7,4,4,2))\n";
    print $ofh "# png(file=\"my_cluster_plots.png\");\n"; 
    print $ofh "for (i in 1:length(files)) {\n";
    print $ofh "    data = read.table(files[i], header=T, row.names=1)\n";
    print $ofh "    ymin = min(data); ymax = max(data);\n";
    print $ofh "    plot_label = paste(files[i], ', ', length(data[,1]), \" trans\", sep='')\n";
    print $ofh "    plot(as.numeric(data[1,]), type='l', ylim=c(ymin,ymax), main=plot_label, col='lightgray', xaxt='n', xlab='', ylab='centered log2(fpkm+1)')\n";
    print $ofh "    axis(side=1, at=1:length(data[1,]), labels=colnames(data), las=2)\n";
    print $ofh "    for(r in 2:length(data[,1])) {\n";
    print $ofh "        points(as.numeric(data[r,]), type='l', col='lightgray')\n";
    print $ofh "    }\n";
    print $ofh "    points(as.numeric(colMeans(data)), type='o', col='blue')\n";
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
