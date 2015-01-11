#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 venn.class\n\n";

my $venn_class = $ARGV[0] or die $usage;

main: {

    print join("\t", "#acc", "up", "down", "up+down", "both", "UP_list", "Down_list") . "\n";

    open (my $fh, $venn_class) or die $!;
    while (<$fh>) {
        chomp;
        my ($acc, $class_list_up, $class_list_down) = split(/\t/);
        
        my @class_up = split(/,/, $class_list_up);
        my @class_down = split(/,/, $class_list_down);

        my %both;
        for (@class_up, @class_down) {
            $both{$_}++;
        }
        

        print join("\t", $acc, scalar(@class_up), scalar(@class_down), 
                   scalar(@class_up) + scalar(@class_down),
                   scalar(keys %both),
                   
                   $class_list_up, $class_list_down) . "\n";
    }

    close $fh;

    exit(0);
}

