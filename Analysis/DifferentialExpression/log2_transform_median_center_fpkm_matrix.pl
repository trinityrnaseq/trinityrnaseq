#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 matrix.fpkm\n\n";

my $matrix = $ARGV[0] or die $usage;


main: {
    
    open (my $fh, $matrix) or die "Error, cannot open file $matrix";
    my $header = <$fh>;
    print $header;
    while (<$fh>) {
        chomp;
        my ($acc, @vals) = split(/\t/);
        
        foreach my $val (@vals) {
            $val = log($val+1)/log(2);
        }

        my $median = &median(@vals);

        foreach my $val (@vals) {
            $val -= $median;
            $val = sprintf("%.2f", $val);
        }

        print join("\t", $acc, @vals) . "\n";
    }
    close $fh;

    exit(0);
    
}

####
sub median {
    my @nums = @_;
    
    @nums = sort {$a<=>$b} @nums;
        
    my $count = scalar (@nums);
    if ($count %2 == 0) {
        ## even number:
        my $half = $count / 2;
        return ( ($nums[$half-1] + $nums[$half]) / 2);
    }
    else {
        ## odd number. Return middle value
        my $middle_index = int($count/2);
        return ($nums[$middle_index]);
    }
}
