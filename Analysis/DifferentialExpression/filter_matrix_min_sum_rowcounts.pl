#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 matrix.fpkm min_sum_counts\n\n";

my $fpkm_matrix = $ARGV[0] or die $usage;
my $min_sum = $ARGV[1] or die $usage;

open (my $fh, $fpkm_matrix) or die $!;
my $header = <$fh>;
print $header;

while (<$fh>) {
    my $line = $_;
    chomp;
    my @x = split(/\t/);
    my $acc = shift @x;
    my $val = shift @x;

    foreach my $other (@x) {
        $val += $other;
    }
    
    if ($val >= $min_sum) {
        print $line;
    }
    
}

exit(0);

