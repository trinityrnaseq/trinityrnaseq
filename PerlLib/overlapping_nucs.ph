#!/usr/local/bin/perl

####
sub nucs_in_common {
    my ($e5, $e3, $g5, $g3) = @_;
	

    ($e5, $e3) = sort {$a<=>$b} ($e5, $e3);
    ($g5, $g3) = sort {$a<=>$b} ($g5, $g3);

	unless (&coordsets_overlap([$e5,$e3], [$g5, $g3])) {
		return(0);
	}
	
    my $length = abs ($e3 - $e5) + 1;
    my $diff1 = ($e3 - $g3);
    $diff1 = ($diff1 > 0) ? $diff1 : 0;
    my $diff2 = ($g5 - $e5);
    $diff2 = ($diff2 > 0) ? $diff2 : 0;
    my $overlap_length = $length - $diff1 - $diff2;
    return ($overlap_length);
}

####
sub coordsets_overlap {
	my ($coordset_A_aref, $coordset_B_aref) = @_;

	my ($lend_A, $rend_A) = sort {$a<=>$b} @$coordset_A_aref;

	my ($lend_B, $rend_B) = sort {$a<=>$b} @$coordset_B_aref;

	if ($lend_A <= $rend_B && $rend_A >= $lend_B) {
		## yes, overlap
		return (1);
	}
	else {
		return (0);
	}
}


####
sub coordset_A_encapsulates_B {
	my ($coordset_A_aref, $coordset_B_aref) = @_;

	my ($lend_A, $rend_A) = sort {$a<=>$b} @$coordset_A_aref;
	
	my ($lend_B, $rend_B) = sort {$a<=>$b} @$coordset_B_aref;
	
	if ($lend_A <= $lend_B && $rend_B <= $rend_A) {
		return(1); # true
	}
	else {
		return(0); # false;
	}
}



1;

