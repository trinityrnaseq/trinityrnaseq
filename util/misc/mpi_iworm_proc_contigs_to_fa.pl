#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 tmp.MPIiworm.rank.txt\n\n";

my $file = $ARGV[0] or die $usage;

main: {
    open (my $fh, $file) or die $!;
    
    my $seq_counter = 0;
    while (<$fh>) {
        chomp;
        $seq_counter++;

        print ">s$seq_counter\n$_\n";
        
    }
    
    
    exit(0);
}

        
