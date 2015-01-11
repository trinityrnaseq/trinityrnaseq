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
            if (/^comp/) {
                /^(comp\d+_c\d+)(_seq\d+)/ or die "Error, cannot parse trinity accession";
                my $comp_id = $1;
                my $trans_id = $2;
                $trans_id = $comp_id . $trans_id;

                my @x = split(/\t/);
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


    
            
