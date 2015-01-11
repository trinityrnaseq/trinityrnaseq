#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 chrysalis_component_listing.txt\n\n";

my $comp_list_file = $ARGV[0] or die $usage;


main: {

    my $ret_val = 0;

    open (my $fh, $comp_list_file) or die "Error, cannot open file $comp_list_file";
    while (<$fh>) {
        chomp;
        my ($comp_id, $comp_base) = split(/\t/);
        my $butterfly_fasta_file = "$comp_base.graph.allProbPaths.fasta";
        
        if (-e $butterfly_fasta_file) {
            open (my $bfly_fh, $butterfly_fasta_file) or die "Error, cannot open file $butterfly_fasta_file";
            while (<$bfly_fh>) {
                print $_;
            }
            close $bfly_fh;
        }
        else {
            print STDERR "Error, no fasta file reported as: $butterfly_fasta_file\n";
            $ret_val = 1;
        }
    }
    
    close $fh;
    
    exit($ret_val);
}

