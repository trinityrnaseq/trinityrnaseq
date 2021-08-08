#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 chrysalis_component_listing.txt min_seq_length\n\n";

my $comp_list_file = $ARGV[0] or die $usage;
my $min_seq_length = $ARGV[1] or die $usage;


main: {

    my $ret_val = 0;

    my %seen;
    
    open (my $fh, $comp_list_file) or die "Error, cannot open file $comp_list_file";
    while (<$fh>) {
        chomp;
        my ($comp_id, $comp_base) = split(/\t/);
        my $butterfly_fasta_file = "$comp_base.graph.allProbPaths.fasta";
        
        if (-e $butterfly_fasta_file) {

            my $fasta_reader = new Fasta_reader($butterfly_fasta_file);
            while (my $seq_obj = $fasta_reader->next()) {
                my $sequence = $seq_obj->get_sequence();
                if (length($sequence) >= $min_seq_length) {

                    if ($seen{$sequence}) {
                        print STDERR "-duplicate sequence detected, excluding it.\n";
                        next;
                    }
                    else {
                        $seen{$sequence} = 1;
                    }
                    
                    print $seq_obj->get_FASTA_format();
                }
            }
        }
        else {
            print STDERR "Error, no fasta file reported as: $butterfly_fasta_file\n";
            $ret_val = 1;
        }
    }
    
    close $fh;
    
    exit($ret_val);
}

