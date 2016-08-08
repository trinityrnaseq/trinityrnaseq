#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use SAM_reader;
use SAM_entry;

use Carp;
use Data::Dumper;

my $usage = "usage: $0 file.sam number_records_to_check=10\n\n";


my $sam_file = $ARGV[0] or die $usage;
my $num_records_to_check = $ARGV[1] or die $usage;

main: {

    my $sam_reader = new SAM_reader($sam_file);

    my $num_records_checked = 0;

    while (my $sam_entry = $sam_reader->get_next()) {
        
        my $core_read_name = $sam_entry->get_read_name();

        if (! $sam_entry->is_paired()) {
            confess "ERROR, only paired reads should exist in bam file.  Encountered unpaired read: " . Dumper($sam_entry);
        }
        
        $num_records_checked++;

        if ($num_records_checked >= $num_records_to_check) {
            last;
        }
    }
    
    exit(0); # all good!

}

