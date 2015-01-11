#!/usr/bin/env perl

use strict;
use warnings;

my @entries;

while (<>) {

    if (/>(\w+).*path=\[(.*)\]/) {
        
        my $acc = $1;
        my $path = $2;

        my @node_descr = split(/\s+/, $path);

        my @nodes;
        
        foreach my $node (@node_descr) {
            $node =~ s/:.*$//;
            push (@nodes, $node);
        }

        push (@entries, {
            acc => $acc,
            nodes => join(" ", @nodes),
        });
    }
}


@entries = sort {$a->{nodes} cmp $b->{nodes}}  @entries;


foreach my $entry (@entries) {
    print join("\t", $entry->{acc}, $entry->{nodes}) . "\n";
}



exit(0);


        
