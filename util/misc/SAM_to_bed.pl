#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 file.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;

main: {

	my $sam_reader = new SAM_reader($sam_file);

	while ($sam_reader->has_next()) {

		my $sam_entry = $sam_reader->get_next();

		if ($sam_entry->is_query_unmapped()) {
			next;
		}

		my $read_name = $sam_entry->get_read_name();
		my $scaff_name = $sam_entry->get_scaffold_name();
		
		my $strand = $sam_entry->get_query_strand();

		my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();

		my @coords;
		foreach my $segment (@$genome_coords_aref) {
			my ($lend, $rend) = @$segment;

			push (@coords, $lend, $rend);
		}

		@coords = sort {$a<=>$b} @coords;
		my $span_lend = shift @coords;
		my $span_rend = pop @coords;

		my @lengths;
		my @starts;
		my $num_segments = 0;
		foreach my $segment (@$genome_coords_aref) {
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
				   $read_name,
				   0,
				   $strand,
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
