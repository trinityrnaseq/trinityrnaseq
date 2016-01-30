#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Nuc_translator;

my $usage = "usage: $0 file.sam\n\n";
my $sam_file = $ARGV[0] or die $usage;

main: {


    my $sam_reader = new SAM_reader($sam_file);

    while (my $sam_entry = $sam_reader->get_next()) {

        my $read_name = $sam_entry->reconstruct_full_read_name();
        my $sequence = $sam_entry->get_sequence();
        
        if ((! $sam_entry->is_query_unmapped()) && $sam_entry->get_query_strand() eq '-') {
            $sequence = &reverse_complement($sequence);
        }

        print ">$read_name\n$sequence\n";
    }


    exit(0);
}

