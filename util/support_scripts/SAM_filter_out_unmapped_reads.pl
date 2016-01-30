#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 file.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;



main: {

	my $sam_reader = new SAM_reader($sam_file);
    
    
    my $filtered_count = 0;
    my $total_count = 0;
    

	while ($sam_reader->has_next()) {
		
		my $sam_entry = $sam_reader->get_next();
        $total_count++;
        
        if ($sam_entry->is_query_unmapped()) {
            $filtered_count++;
        }
        else {
            print $sam_entry->toString() . "\n";
        }
    }

    print STDERR "-filtered $filtered_count of $total_count SAM entries as unaligned = " . sprintf("%.2f", $filtered_count / $total_count * 100) . "\% of SAM file corresponding to unmapped reads.  Note, aligned reads may be counted multiple times due to multiple-mappings, so not the same as percent of unaligned reads.\n";
    
	exit(0);

}
