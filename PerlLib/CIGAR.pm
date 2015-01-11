package CIGAR;

use strict;
use warnings;

use Nuc_translator;

####
sub construct_cigar {
	my ($genome_coords_aref, $query_coords_aref, $read_length, # required
		$genome_seq_sref, $strand # optional
		) = @_;

	my $cigar = "";
		

	for (my $i = 0; $i <= $#$genome_coords_aref; $i++) {
		
		my ($curr_genome_lend, $curr_genome_rend) = @{$genome_coords_aref->[$i]};
		my ($curr_query_lend, $curr_query_rend) = @{$query_coords_aref->[$i]};

		if ($i == 0) {

			if ($curr_query_lend > 1) {
				$cigar .= ($curr_query_lend - 1) . "S";
			}
		}
		else {
			my ($prev_genome_lend, $prev_genome_rend) = @{$genome_coords_aref->[$i-1]};
			my ($prev_query_lend, $prev_query_rend) = @{$query_coords_aref->[$i-1]};
			
			if ( (my $delta_genome = $curr_genome_lend - $prev_genome_rend) > 1) {
				my $deletion_intron_char = ($genome_seq_sref && $strand) ? &_check_intron_consensus($prev_genome_rend, $curr_genome_lend, $genome_seq_sref, $strand) : 'D'; # intron or deletion?
				$cigar .= ($delta_genome-1) . $deletion_intron_char;
			}
			if ( (my $delta_query = $curr_query_lend - $prev_query_rend) > 1) {
				$cigar .= ($delta_query-1) . "I";
			}
		}
		

		my $len = $curr_genome_rend - $curr_genome_lend + 1;
		
		$cigar .= "$len" . "M";
		
		if ($i == $#$genome_coords_aref) {
			
			if ($curr_query_rend < $read_length) {
				$cigar .= ($read_length - $curr_query_rend) . "S";
			}
		}
		
	}

	return($cigar);
}


####
sub _check_intron_consensus {
	my ($left_segment_bound, $right_segment_bound, $genome_sref, $strand) = @_;

	my $left_dinuc = uc substr($$genome_sref, $left_segment_bound, 2);

	my $right_dinuc = uc substr($$genome_sref, $right_segment_bound-3, 2);

	if ($strand eq '-') {
		my ($left_dinuc_copy, $right_dinuc_copy) = ($left_dinuc, $right_dinuc);
		$left_dinuc = &reverse_complement($right_dinuc_copy);
		$right_dinuc = &reverse_complement($left_dinuc_copy);
	}
	
	if (  
		( ($left_dinuc eq "GT" || $left_dinuc eq "GC") && $right_dinuc eq "AG")
		||
		($left_dinuc eq "CT" && $right_dinuc eq "AC")
		) {

		return("N"); # got splice pair
	}
	else {
		return("D"); # deletion
	}
	

}








1; #EOM
