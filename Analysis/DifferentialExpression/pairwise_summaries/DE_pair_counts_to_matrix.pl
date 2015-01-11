#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 DE_counts.pairs\n\n";

my $pairs_info = $ARGV[0] or die $usage;

main: {

    my %sample_pair_to_counts;

    open (my $fh, $pairs_info) or die $!;
    while (<$fh>) {
        chomp;
        my ($sampleA, $sampleB, $count) = split(/\t/);
        $sample_pair_to_counts{$sampleA}->{$sampleB} = $count;
        $sample_pair_to_counts{$sampleB}->{$sampleA} = $count;
    }
    close $fh;


    my @samples = sort (keys %sample_pair_to_counts);

    print "\t" . join("\t", @samples) . "\n";
    foreach my $sampleA (@samples) {
        print "$sampleA";
        foreach my $sampleB (@samples) {
            my $count = $sample_pair_to_counts{$sampleA}->{$sampleB} || 0;
            print "\t$count";
        }
        print "\n";
    }


    exit(0);
}


