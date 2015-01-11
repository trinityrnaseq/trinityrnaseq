#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;


my $usage = "usage: $0 file.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;

main: {



    my $sam_reader = new SAM_reader($sam_file);

    print "chr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\n";
    
    my @entries;
    my $curr_read_name = "";

    while (my $sam_entry = $sam_reader->get_next()) {
        my $core_read_name = $sam_entry->get_core_read_name();
        if ($core_read_name ne $curr_read_name) {
            if (scalar(@entries) > 1) {
                &process_entries(@entries);
            }
            @entries = ();
            
        }
        $curr_read_name = $core_read_name;
        push (@entries, $sam_entry);
    }
    
    # get last ones
    if (scalar(@entries) > 1) {
        &process_entries(@entries);
        
    }
    
    exit(0);
}

####
sub process_entries {
    my @entries = @_;

    my @left_entries;
    my @right_entries;

    foreach my $entry (@entries) {
        
        if ($entry->get_scaffold_name() eq "*") { next; } # unaligned read
        
        if ($entry->is_first_in_pair()) {
            push (@left_entries, $entry);
        }
        elsif ($entry->is_second_in_pair()) {
            push (@right_entries, $entry);
        }
        else {
            die "Error, not left or right seq of pair: " . Dumper($entry);
        }
    }

    if (scalar(@left_entries) > 1 || scalar(@right_entries) > 1) {
        print STDERR "Error, should only have single read mappings, but have more than one:"
            . Dumper(\@left_entries) . Dumper(\@right_entries);
        
        return;
    }

    if (scalar(@left_entries) == 1 && scalar(@right_entries) == 1) {
        ## all good.
        
        my $left_entry = $left_entries[0];
        my $right_entry = $right_entries[0];

        # "chr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\n";

        print join("\t", 
                   $left_entry->get_scaffold_name(), $left_entry->get_scaffold_position(), $left_entry->get_query_strand(),
                   $right_entry->get_scaffold_name(), $right_entry->get_scaffold_position(), $right_entry->get_query_strand())
            . "\n";
        
    }
    
    
    return;
}
