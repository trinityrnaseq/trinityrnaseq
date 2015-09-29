#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 blastx.outfmt6 trinity.matrix\n\n";

my $blastx_file = $ARGV[0] or die $usage;
my $trinity_matrix = $ARGV[1] or die $usage;

main: {

    my %blastx_top_hit;  # store for transcript and component
    
    {
        open (my $fh, $blastx_file) or die $!;
        while (<$fh>) {
            if (/^(.*|c\d+)_/) {
                
                my $comp_id = $1;
                
                my @x = split(/\t/);
                my $trans_id = $x[0];
                my $hit = $x[1];
                $blastx_top_hit{$comp_id} = $hit;
                $blastx_top_hit{$trans_id} = $hit;
            }
        }
        close $fh;
    }

    
    open (my $fh, $trinity_matrix) or die "Error, cannot open file $trinity_matrix";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        if (my $hit = $blastx_top_hit{$acc}) {
            $x[0] = "$acc|$hit";
        }

        print join("\t", @x) . "\n";
    }
    close $fh;


    exit(0);
}


    
            
