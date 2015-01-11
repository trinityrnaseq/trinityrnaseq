#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;


my $usage = "usage: $0 coordSorted.sam max_per_position\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $max_per_pos = $ARGV[1] or die $usage;

main: {

    my $sam_reader = new SAM_reader($sam_file);


    my $curr_scaffold = "";
    my $curr_pos = -1;

    my $counter = 0;

    while (my $sam_entry = $sam_reader->get_next()) {


        if ($sam_entry->is_query_unmapped()) {
            next;
        }
        
        my $scaff = $sam_entry->get_scaffold_name();
        my $pos = $sam_entry->get_scaffold_position();

        if ($curr_scaffold ne $curr_scaffold || $curr_pos != $pos) {
            
            $counter = 0; # re init
            $curr_scaffold = $scaff;
            $curr_pos = $pos;
            
            
        }

        $counter++;

        if ($counter <= $max_per_pos) {
            
            print $sam_entry->toString() . "\n";
        }

    }


    exit(0);
}




