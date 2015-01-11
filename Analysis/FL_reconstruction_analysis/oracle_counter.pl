#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 results.oracle\n\n";

my $oracle_results = $ARGV[0] or die $usage;

main: {

    my %genes;
    my %isoforms;

    open (my $fh, $oracle_results) or die "Error, cannot open file $oracle_results";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $TF = $x[4];

        if ($TF eq "T") {
            my ($trans, $gene) = split(/;/, $acc);
            $genes{$gene}++;
            
            $isoforms{$acc}++;
        }
    }

    close $fh;

    my $num_genes = scalar(keys %genes);
    my $num_trans = scalar(keys %isoforms);

    print "\n\nOracle_counts:\t$num_genes\t$num_trans\n\n";

    exit(0);
    
}
