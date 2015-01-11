#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 summary_file num_top_ranks\n\n";

my $summary_file = $ARGV[0] or die $usage;
my $num_top_ranks = $ARGV[1] or die $usage;

main: {


    my %method_pairs;
    
    open (my $fh, $summary_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $methods = $x[3];
        my $rankings = $x[4];

        my @meth = split(/,/, $methods);
        my @ranks = split(/,/, $rankings);

        for (my $i = 0; $i < $#meth; $i++) {
            
            my $method_i = $meth[$i];
            my $rank_i = $ranks[$i];

            for (my $j = $i + 1; $j <= $#meth; $j++) {

                my $method_j = $meth[$j];
                my $rank_j = $ranks[$j];

                if ($rank_i <= $num_top_ranks || $rank_j <= $num_top_ranks) {

                    my $method_pair = join(",", sort ($method_i, $method_j));
                    
                    push (@{$method_pairs{$method_pair}}, {  $method_i => $rank_i,
                                                             $method_j => $rank_j,
                                                         });

                }
            }
        }
    }
    close $fh;

    
    foreach my $method_pair (keys %method_pairs) {
        my ($method_i, $method_j) = split(/,/, $method_pair);

        open (my $ofh, ">meth_compare.$method_i-$method_j.ranks") or die $!;
        print $ofh join("\t", $method_i, $method_j) . "\n";
        
        my @structs = @{$method_pairs{$method_pair}};
        
        foreach my $struct (@structs) {
            my $rank_i = $struct->{$method_i};
            my $rank_j = $struct->{$method_j};

            print $ofh join("\t", $rank_i, $rank_j) . "\n";
        }
    }


    exit(0);
}
