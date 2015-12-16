#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 RSEM.isoforms.fpkm\n\n";

my $fpkm_file = $ARGV[0] or die $usage;

my $Rscript = "$fpkm_file.R";
open (my $ofh, ">$Rscript");
print $ofh "source(\"$FindBin::RealBin/R/expression_analysis_lib.R\")\n";
print $ofh "pdf(\"$fpkm_file.genes_vs_minFPKM.pdf\")\n";
print $ofh "plot_expressed_gene_counts(\"$fpkm_file\", title=\"expressed transcripts vs. min FPKM\", fpkm_range=seq(0,5,0.01), outfile=\"$fpkm_file.genes_vs_minFPKM.dat\")\n";
print $ofh "dev.off()\n";
close $ofh;

&process_cmd("R --vanilla -q < $Rscript");


exit(0);


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
