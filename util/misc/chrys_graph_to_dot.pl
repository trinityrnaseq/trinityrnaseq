#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 comp.raw.graph\n\n";

my $chrys_graph = $ARGV[0] or die $usage;

main: {

    open (my $fh, $chrys_graph) or die "Error, cannot open file $chrys_graph";
    
    print "digraph G {\n";

    while (<$fh>) {
        chomp;
        my ($id, $prev_id, @rest) = split(/\t/);
        if (defined ($id) && defined($prev_id)) {
            print "    $prev_id->$id\n";
        }
    }
    print "}\n";

    exit(0);
}


        
