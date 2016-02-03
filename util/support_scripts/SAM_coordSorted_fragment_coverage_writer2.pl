#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\nusage: $0 coordSorted.file.sam [max_pair_distance=10000]\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $MAX_PAIR_DIST = $ARGV[1] || 10000;


main: {
	
	my @coverage;
	
	my $sam_reader = new SAM_reader($sam_file);
	
	my $current_scaff = "";
	
	my $counter = 0;
	while ($sam_reader->has_next()) {

		$counter++;
		if ($counter % 10000 == 0) {
			print STDERR "\r[$counter]  ";
		}

		my $sam_entry = $sam_reader->get_next();
		
		if ($sam_entry->get_scaffold_name() ne $current_scaff) {
			&report_coverage(\@coverage) if @coverage;
			@coverage = ([1,0]); # reinit
			$current_scaff = $sam_entry->get_scaffold_name();
			print "variableStep chrom=$current_scaff\n";
		}
		
		&add_coverage($sam_entry, \@coverage);
		

	}

	if (@coverage) {
		&report_coverage(\@coverage);
	}
		

	exit(0);

}


####
sub report_coverage {
	my ($coverage_aref) = @_;
	
	foreach my $cov_info (@$coverage_aref) {
		my ($pos, $cov) = @$cov_info;

		print join("\t", $pos, $cov) . "\n";
	}

	return;
}





####
sub add_coverage {
	my ($sam_entry, $coverage_aref) = @_;

	my $scaffold = $sam_entry->get_scaffold_name();
	my $scaff_pos = $sam_entry->get_scaffold_position();
	
	if (@$coverage_aref) {
		my $last_pos;
		while (@$coverage_aref && $coverage_aref->[0]->[0] < $scaff_pos) {
			my $pos_info = shift @$coverage_aref;
			print join("\t", @$pos_info) . "\n";
			$last_pos = $pos_info->[0];
		}
		
		## add any intervening positions to wig.
		if ($last_pos) {
			$last_pos++;
			while ($last_pos < $scaff_pos) {
				print "$last_pos\t0\n";
				$last_pos++;
			}
		}
	}
	
	unless (@$coverage_aref) {
		# prime it
		push (@$coverage_aref, [$scaff_pos, 0]);
	}
	
	
	my @genome_coords;
	
	my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();
	
	foreach my $coordset (@$genome_coords_aref) {
		my ($lend, $rend) = @$coordset;
		push (@genome_coords, $lend, $rend);
	}
	
	if ($sam_entry->is_paired() && $sam_entry->is_proper_pair()) {
		
		my $mate_scaffold = $sam_entry->get_mate_scaffold_name();
		my $mate_position = $sam_entry->get_mate_scaffold_position();
		
		
		if ( ($mate_scaffold eq $scaffold  || $mate_scaffold eq '=')
			 && 
			 $mate_position > $scaff_pos
			 &&
			 $mate_position - $scaff_pos <= $MAX_PAIR_DIST) {
			
			push (@genome_coords, $mate_position-1);
		}
	}
	
	@genome_coords = sort {$a<=>$b} @genome_coords;

	my $lend = shift @genome_coords;
	my $rend = pop @genome_coords;
	
	
	if ($rend - $lend + 1 < $MAX_PAIR_DIST) {  ## same deal for long introns.

		my $current_first_pos = $coverage_aref->[0]->[0];
		my $current_last_pos = $coverage_aref->[$#$coverage_aref]->[0];

		#print "Current first pos: $current_first_pos\ncurrent last: $current_last_pos\n";
		
		## add positions for currently missing entries.
		for (my $i = $current_last_pos + 1; $i <= $rend; $i++) {
			push (@$coverage_aref, [$i, 0]);
		}
		
		my $first_index = $lend - $current_first_pos;
		my $last_index = $rend - $current_first_pos;

		## add coverage:
		for (my $i = $first_index; $i <= $last_index; $i++) {
			
			$coverage_aref->[$i]->[1]++;
			
		}
	}
	
	return;
}

