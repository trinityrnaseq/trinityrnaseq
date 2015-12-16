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
        
        if (! $sam_entry->is_proper_pair()) {
            $filtered_count++;
        }
        else {
            print $sam_entry->toString() . "\n";
        }
    }
    
    print STDERR "-filtered $filtered_count of $total_count indiv read alignments as not proper pairs = " . sprintf("%.2f", $filtered_count / $total_count * 100) . "\% of alignments\n";
    
	exit(0);

}
