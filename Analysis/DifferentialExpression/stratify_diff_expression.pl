#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

my $usage = "\n\nusage: $0 edgeR_dir\n\n";

my $edgeR_dir = $ARGV[0] or die $usage;

unless (-d $edgeR_dir) {
    die "Error, require a directory name for the parameter.\n\n$usage\n";
}

my ($fpkm_matrix_file) = <$edgeR_dir/*FPKM>;
unless ($fpkm_matrix_file) {
    die "Error, cannot locate the FPKM matrix file at $edgeR_dir";
}

for my $fold_change (1..8) {

    for my $pvalue (2..10) {

        my $cmd = "$FindBin::RealBin/analyze_diff_expr.pl --matrix $fpkm_matrix_file -C $fold_change -P 1e-$pvalue";
        &process_cmd($cmd);
        
        


    }
    
}

exit(0);

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }
    
    return;
}
