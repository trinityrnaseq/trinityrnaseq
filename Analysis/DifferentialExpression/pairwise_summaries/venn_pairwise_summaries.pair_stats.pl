#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "venn_pairwise_out.dat\n\n";

my $venn_txt = $ARGV[0] or die $usage;

main: {

    my %combos;
    
    my %sample_pair_to_combo_count;

    open (my $fh, $venn_txt) or die $!;

    while (<$fh>) {
        chomp;
        my ($feature, $sampleA, $sampleB, $combo) = split(/\t/);
        
        my @combo_renamed;
        foreach my $c (split(/,/, $combo)) {
            $c =~ /^([^\.]+)/ or die "Error, cannot parse combo $c";
            push (@combo_renamed, $1);
        }
        $combo = join(",", @combo_renamed);
        
        my $sample_pair = join(",", $sampleA, $sampleB);

        $sample_pair_to_combo_count{$sample_pair}->{$combo}++;
        $combos{$combo}++;
        
    }
    close $fh;

    my @combos = sort keys %combos;

    print join("\t", "#sample", @combos) . "\n";

    foreach my $sample_pair (sort keys %sample_pair_to_combo_count) {

        print $sample_pair;
        foreach my $combo (@combos) {
            my $val = $sample_pair_to_combo_count{$sample_pair}->{$combo} || 0;
            print "\t$val";
        }
        print "\n";
        
    }
    
    exit(0);
}


