#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 nameSorted.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;



main: {

	my $sam_reader = new SAM_reader($sam_file);
    
    
    my $filtered_count = 0;
    my $total_count = 0;
    
    my $prev_core_read_name = "";
    my @entries;
    
	while ($sam_reader->has_next()) {
		
		my $sam_entry = $sam_reader->get_next();
        
        my $core_read_name = $sam_entry->get_core_read_name();
        if ($core_read_name ne $prev_core_read_name) {
            &process_entries(@entries);
            @entries = ();
        }
        
        push (@entries, $sam_entry);
    
        $prev_core_read_name = $core_read_name;
    }
    

    if (@entries) {
        &process_entries(@entries);
    }
    

	exit(0);

}


####
sub process_entries {
    my @entries = @_;

    my %scaffolds;
    foreach my $entry (@entries) {
        my $scaff = $entry->get_scaffold_name();
        $scaffolds{$scaff}++;
    }

    my $num_scaff = scalar(keys %scaffolds);
    if ($num_scaff == 1) {
        foreach my $entry (@entries) {
            print $entry->toString() . "\n";
        }
    }

    return;
}
