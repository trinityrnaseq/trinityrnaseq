#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Overlap_piler;

my $usage = "usage: $0 file.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;

main: {

	my $sam_reader = new SAM_reader($sam_file);

    my %core_to_coords;
    
	while ($sam_reader->has_next()) {

		my $sam_entry = $sam_reader->get_next();

		if ($sam_entry->is_query_unmapped()) {
			next;
		}

		my $read_name = $sam_entry->get_read_name();
		my $scaff_name = $sam_entry->get_scaffold_name();
		
		my $strand = $sam_entry->get_query_strand();

		my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();

        
        my $core_scaff = join("$;", $sam_entry->get_core_read_name(), $scaff_name);

		my @coords;
		foreach my $segment (@$genome_coords_aref) {
			my ($lend, $rend) = @$segment;
            
            push (@{$core_to_coords{$core_scaff}}, [$lend, $rend]);
		}

    }
    
    foreach my $core_scaff (keys %core_to_coords) {

        my @coords = @{$core_to_coords{$core_scaff}};
        my ($read_core_name, $scaff_name) = split(/$;/, $core_scaff);
        
		@coords = sort {$a<=>$b} @coords;
		@coords = &Overlap_piler::simple_coordsets_collapser(@coords);
        

        my $span_lend = $coords[0]->[0];
		my $span_rend = $coords[$#coords]->[1];
        
		my @lengths;
		my @starts;
		my $num_segments = 0;
		foreach my $segment (@coords) {
			my ($lend, $rend) = @$segment;
			
			my $length = $rend - $lend + 1;
			push (@lengths, $length);
			push (@starts, $lend - $span_lend);
			$num_segments++;
		}
        
		$span_lend--; # coordinate is zero-based, and rend is exclusive

		print join("\t", 
				   $scaff_name,
				   $span_lend,
				   $span_rend,
				   $read_core_name,
				   0,
				   "+",
				   $span_lend,
				   $span_rend,
				   ".",
				   $num_segments,
				   join(",", @lengths),
				   join(",", @starts),
				   ) . "\n";
	}


	exit(0);
}
