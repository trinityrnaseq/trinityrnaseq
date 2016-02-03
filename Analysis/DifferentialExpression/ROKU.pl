#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use BHStats;

my $usage = "usage: $0 fpkm.matrix\n\n";

my $matrix = $ARGV[0] or die $usage;


my $pseudocount = 0.1;

main: {

    open (my $fh, $matrix) or die $!;
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        
        my @x = split(/\t/);
        
        my $gene_id = shift @x;

        my $bi = &tukey_biweight(@x);
        
        my $sum = 0;

        foreach my $val (@x) {
            $val = abs($val - $bi) + $pseudocount;
            $sum += $val + $pseudocount;
        }
        
        my $entropy = 0;
        foreach my $val (@x) {
            my $ratio = $val/$sum;
            $entropy += $ratio * log($ratio)/log(2);
        }
        
        $entropy *= -1;

        print "$gene_id\t$entropy\n";
    }
    
    exit(0);
}

