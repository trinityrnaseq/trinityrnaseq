#!/usr/bin/env perl

package Overlap_info;

use strict;
use warnings;
use List::Util qw (min max);
use Carp;
use Data::Dumper;


## works on coordinate pairs:  [a1, a2], [b1, b2]
sub overlap {
	my ($coordsA_aref, $coordsB_aref) = @_;

	my ($lendA, $rendA) = sort {$a<=>$b} @$coordsA_aref;
	
	my ($lendB, $rendB) = sort {$a<=>$b} @$coordsB_aref;

	if ($lendA <= $rendB && $rendA >= $lendB) {
		return(1);
	}
	else {
		return(0);
	}
}


sub contains {
	my ($larger, $smaller) = @_;

	my ($smaller_lend, $smaller_rend) = sort {$a<=>$b} @$smaller;
	my ($larger_lend, $larger_rend) = sort {$a<=>$b} @$larger;

	if ($smaller_lend >= $larger_lend && $smaller_rend <= $larger_rend) {
		return(1);
	}
	
	else {
		#print "no containment: [$larger_lend-$larger_rend] no containment of [$smaller_lend-$smaller_rend]\n";
		return(0);
	}
}


## works on coordinate pairs:  [a1, a2], [b1, b2]
sub overlap_length {
	my ($coordsA_aref, $coordsB_aref) = @_;

	if (&overlap($coordsA_aref, $coordsB_aref)) {
		
		my ($lendA, $rendA) = sort {$a<=>$b} @$coordsA_aref;
		
		my ($lendB, $rendB) = sort {$a<=>$b} @$coordsB_aref;
		
		if ($lendA > $lendB) {
			# swap em
			($lendA, $rendA, $lendB, $rendB) = ($lendB, $rendB, $lendA, $rendA);
		}
		
		my $overlap_lend = max($lendA, $lendB);
		my $overlap_rend = min($rendA, $rendB);

		my $overlap_len = $overlap_rend - $overlap_lend + 1;
		return($overlap_len);
				
	}
	else {
		return(0);
	}
}



## works on sets of coordinates: ([ [a1,a2], [a3,a4], ...]) , ([ [b1,b2], [b3,b4], ... ])
sub sum_overlaps {
	my ($coordsets_A_aref, $coordsets_B_aref) = @_;


	## Be sure that no intra-set coordinates overlap each other.  TODO: force this sanity check
	
	my $sum = 0;
	
	foreach my $coordset_A_aref (@$coordsets_A_aref) {

		foreach my $coordset_B_aref (@$coordsets_B_aref) {

			$sum += &overlap_length($coordset_A_aref, $coordset_B_aref);
		}
	}

	return($sum);
}


## overlap, compatible, and A contains B
sub compatible_overlap_A_contains_B {
    my ($coordsets_A_aref, $coordsets_B_aref) = @_;

    my ($A_lend, $A_rend) = &get_coordset_span(@$coordsets_A_aref);
    my ($B_lend, $B_rend) = &get_coordset_span(@$coordsets_B_aref);
    
    
    if (&contains([$A_lend, $A_rend], [$B_lend, $B_rend])
        &&
        &compatible_overlap($coordsets_A_aref, $coordsets_B_aref) ) {

        return(1);
    }
    else {
        return(0);
    }
}



## overlap and encode identical internal boundaries.
sub compatible_overlap {
	my ($coordsets_A_aref, $coordsets_B_aref) = @_;

	my $verbose_check = 0;
    
	print "* Compatibility check between: " . Dumper($coordsets_A_aref) . " and " . Dumper($coordsets_B_aref) if $verbose_check;
	

	my @A_coordsets = &_order_coordsets(@$coordsets_A_aref);
	my @B_coordsets = &_order_coordsets(@$coordsets_B_aref);

	my @A_span = &get_coordset_span(@A_coordsets);
	my @B_span = &get_coordset_span(@B_coordsets);
	
	unless (&overlap(\@A_span, \@B_span)) {
		print "-no overlap of spans: @A_span, @B_span\n" if $verbose_check;
		return(0);
	}
	

	my ($i, $j);
	my $found_overlap_flag = 0;
	
  overlap_search:
	for ($i = 0; $i <= $#A_coordsets; $i++) {
		
		for ($j = 0; $j <= $#B_coordsets; $j++) {
			
			if (&overlap($A_coordsets[$i], $B_coordsets[$j])) {
				
				
				$found_overlap_flag = 1;
				


				last overlap_search;
			}
		}
	}
	
	unless ($found_overlap_flag) {
		print "-no overlap of segments.\n" if $verbose_check;
		return(0);
	}
	
	## if its not the first segment of A or B, then they're incompatible.
	if (! ($i == 0 || $j == 0) ) {
		print "-overlap doesn't anchor at first segment of either entry. $i, $j\n" if $verbose_check;
		return(0);
	}
	
	
	while ($i <= $#A_coordsets && $j <= $#B_coordsets) {

		my @a_coords = @{$A_coordsets[$i]};
		my @b_coords = @{$B_coordsets[$j]};

	   
		## left junction check.
		if ($i != 0 && $j != 0) {
			## left bounds should have matching edges.
			if ($a_coords[0] != $b_coords[0]) {
				print "-(internal seg) left bounds fail to match: " . Dumper(\@A_coordsets) . Dumper(\@B_coordsets) if $verbose_check;
				return(0);
			}
		}
		## right junction check:
		if ($i != $#A_coordsets && $j != $#B_coordsets) {
			if ($a_coords[1] != $b_coords[1]) {
				print "-(internal seg) right bounds fail to match: " . Dumper(\@A_coordsets) . Dumper(\@B_coordsets) if $verbose_check;
				return(0);
			}
		}
		
		## left bound check
		if ($i == 0 && $j != 0) {
			if ($a_coords[0] < $b_coords[0]) {
				print "-(i first seg), fail left: " . Dumper(\@A_coordsets) . Dumper(\@B_coordsets) if $verbose_check;
				return(0);
			}
		}
		if ($i != 0 && $j == 0) {
			if ($b_coords[0] < $a_coords[0]) {
				print "-(j first seg), fail left: " . Dumper(\@A_coordsets) . Dumper(\@B_coordsets) if $verbose_check;
				return(0);
			}
		}

		## right bound check.
		if ($i == $#A_coordsets && $j != $#B_coordsets) {
			if ($a_coords[1] > $b_coords[1]) {
				print "-(i last seg), fail right: " . Dumper(\@A_coordsets) . Dumper(\@B_coordsets) if $verbose_check;
				return(0);
			}
		}

		if ($i != $#A_coordsets && $j == $#B_coordsets) {
			if ($b_coords[1] > $a_coords[1]) {
				print "-(j laset seg), fail right: " . Dumper(\@A_coordsets) . Dumper(\@B_coordsets) if $verbose_check;
				return(0);
			}
		}
		
		$i++;
		$j++;

	}
	
	## must be compatible
	print "-made it. Must be compatible.\n" if $verbose_check;
	return(1);
	
	

}

####
sub _order_coordsets {
	my (@coordsets) = @_;

	my @ret_coords;

	foreach my $coordpair (@coordsets) {
		my ($lend, $rend) = sort {$a<=>$b} @$coordpair;
		push (@ret_coords, [$lend, $rend]);
	}


	@ret_coords = sort {$a->[0]<=>$b->[0]} @ret_coords;

	return(@ret_coords);
}

sub get_gap_coords {
	my ($coordsets_aref) = @_;

	my @coordsets = &_order_coordsets(@$coordsets_aref);
	
	my @gaps;
	
	for (my $i = 1; $i <= $#coordsets; $i++) {
		
		my $prev_rend = $coordsets[$i-1]->[1];
		my $curr_lend = $coordsets[$i]->[0];
		
		push (@gaps, [$prev_rend + 1, $curr_lend - 1]);
	}

	return(@gaps);
}



sub get_coordset_span {
	my @coordsets = @_;

	my @coords;
	foreach my $coordset (@coordsets) {
		push (@coords, @$coordset);
	}

	my $min_coord = min(@coords);
	my $max_coord = max(@coords);
	
	return($min_coord, $max_coord);
}


sub coordsets_to_string {
    my @coordsets = @_;

    my @text;
    foreach my $coordset (@coordsets) {
        my ($lend, $rend) = @$coordset;
        push (@text, "($lend,$rend)");
    }
    
    my $text_line = join("--", @text);
    return($text_line);
}

    


sub _coords_are_identical {
	my ($coordpair_A_aref, $coordpair_B_aref)  = @_;
	
	if ($coordpair_A_aref->[0] == $coordpair_B_aref->[0]
		&&
		$coordpair_A_aref->[1] == $coordpair_B_aref->[1]) {

		return(1);
	}

	else {
		return(0);
	}
}
	


1; #EOM

