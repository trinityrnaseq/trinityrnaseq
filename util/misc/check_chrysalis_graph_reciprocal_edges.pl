#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 graphFromIwormFasta.out\n\n";

my $file = $ARGV[0] or die $usage;

main: {

    my %graph;
    
    open(my $fh, $file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);

        my $node_id = shift @x;
        my $spacer = shift @x;

        foreach my $other_node (@x) {
            $graph{$node_id}->{$other_node}++;
        }
    }
    close $fh;


    my $found_missing_recip = 0;
    foreach my $node (keys %graph) {

        foreach my $other_node (keys %{$graph{$node}}) {

            if (! exists $graph{$other_node}->{$node}) {
                $found_missing_recip = 1;

                print STDERR "Error, have $node\->$other_node, but missing $other_node\->$node\n";
            
            }
        }
    }

    if ($found_missing_recip) {
        die "Error, missing recips found.\n";
    }
    else {
        print STDERR "\n\nAll good. :-)\n\n";
        
        exit(0);
    }

    
}
