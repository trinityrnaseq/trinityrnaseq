#!/usr/local/bin/perl

package main;
our $SEE;


package CDNA::Genome_based_cDNA_assembler;

=head1 NAME

CDNA::Genome_based_cdna_assembler

=cut

=head1 DESCRIPTION

This module is used to assemble compatible cDNA alignments.  The algorithm is as follows:
must describe this here.

=cut

use strict;
use CDNA::CDNA_alignment;
use Data::Dumper;

my $DELIMETER = "$;,";
my $FUZZLENGTH = 20;


=item new()

=over 4

B<Description:> instantiates a new cDNA assembler obj.

B<Parameters:> $sequence_sref

$sequence_sref is a reference to a scalar containing the genomic sequence string.

B<Returns:> $obj_href

$obj_href is the object reference newly instantiated by this new method.

=back

=cut

sub new () {
    my $package_name = shift;
    my $self = {};
    bless ($self, $package_name);
    $self->_init(@_);
    return ($self);
}


sub _init {
    my $self = shift;
    my ($sequence_ref) = @_;
    $self->{incoming_alignments} = []; #these are the alignments to be assembled.
    $self->{assemblies} = []; #contains list of all singletons and assemblies.
    $self->{sequence_ref} = $sequence_ref;
    $self->{fuzzlength} = $FUZZLENGTH;  #default setting.
}


=item assemble_alignments()

=over 4

B<DESCRIPTION:> assembles a series of cDNA aligmnments into one or more cDNA assemblies

B<Parameters:> @alignments

@alignments is an array of CDNA::CDNA_alignment objects

B<Returns:> none.

=back

=cut

sub assemble_alignments {
    my $self = shift;
    my @alignments = @_;
    @alignments = reverse sort {$a->{length}<=>$b->{length}} @alignments; #keep in order of decreasing alignment length.
    $self->{incoming_alignments} = [@alignments];
    #return;
 
    ## Algorithm: given cdna, merge with rest of cDNAs iteratively until merging complete.
    ##            find unmerged cDNA, merge it if possible with all individual cdna entries.
    ##            continue until all cDNAs have merged products, if possible.
    
    my $num_alignments = $#alignments + 1;
    for (my $i = 0; $i < $num_alignments; $i++) {
	my $seed_alignment = $alignments[$i];
	if ($seed_alignment->{merged}) {next;}
	print "Checking seed $i\n" if $SEE;
	my $merged_alignment = $seed_alignment; #initialize to seed alignment.
	my $merged_flag = 1;
	my $round = 0;
	while ($merged_flag) {
	    $merged_flag = 0;
	    $round++;
	    print "Merging round $round\n" if $SEE;
	    for (my $j=0; $j < $num_alignments; $j++) {
		if ($i == $j) {next;} #no self comparisons.
		my $other_alignment = $alignments[$j];
		if ($self->already_contains($merged_alignment, $other_alignment)) { next;}
		if ($self->can_merge($merged_alignment, $other_alignment)) {
		    my $initial_merged_fli_status = $merged_alignment->is_fli();
		    my $initial_merged_orient = $merged_alignment->get_orientation();
		    $merged_alignment = $self->merge_alignments($merged_alignment, $other_alignment);
		    $alignments[$i]->{merged} = 1; #set seed alignment merge flag.
		    $alignments[$j]->{merged} = 1; #set other alignment merge flag.
		    $merged_flag = 1; #indicates something actually merged this round.
		    
		    $merged_alignment->remap_cdna_segment_coords();
		    print "merged: " . $merged_alignment->toToken() . "\n" if $SEE;
		    
		}
	    }
	}
	push (@{$self->{assemblies}}, $merged_alignment); #either an assembled product, or something that won't ever merge.
    }
}

=item get_assemblies()

=over 4

B<Description:> returns all the alignment assemblies resulting from the assembly procedure.

B<Parameters:> none.

B<Returns:> @assemblies

@assemblies is an array of CDNA::CDNA_alignment objects.

use the get_acc() method of the alignment object to retrieve all the accessions of the cDNAs that were merged into the assembly.

=back

=cut


sub get_assemblies {
    my $self = shift;
    return (@{$self->{assemblies}});
}


# private method.  Determines if two alignments are compatible with one another.
sub can_merge () {
    my $self = shift;
    my ($a1, $a2) = @_;
    print "Checking to see if can merge: " . $a1->get_acc() . ", " . $a2->get_acc() . "\n" if $::SEE;
    ## See if the coord spans overlap
    my ($a1_lend, $a1_rend) = $a1->get_coords();
    my ($a2_lend, $a2_rend) = $a2->get_coords();
    unless (&overlap($a1_lend, $a1_rend, $a2_lend, $a2_rend)) {
	print "failed merge: No overlap between alignment spans. ($a1_lend, $a1_rend) vs. ($a2_lend, $a2_rend)\n" if $::SEE;
	return(0);
    }
    
    ## Make sure the spliced orientation is equivalent if appropriate
    my $a1_num_segs = $a1->get_num_segments();
    my $a2_num_segs = $a2->get_num_segments();
    my $a1_spliced_orientation = $a1->get_spliced_orientation();
    my $a2_spliced_orientation = $a2->get_spliced_orientation();
    my $a1_is_fli = $a1->is_fli();
    my $a2_is_fli = $a2->is_fli();

    my $fuzzlength = $self->{fuzzlength};

    if ($a1_num_segs > 1 && $a2_num_segs > 1) { #if more than one segment, then spliced orientation is relevant.
	if ($a1_spliced_orientation ne $a2_spliced_orientation) {
	    print "failed merge: $a1_num_segs segments vs. $a2_num_segs and opposite spliced orientations.\n" if $::SEE;
	    return (0);
	}
    }
    
    if ($a1_is_fli && $a2_is_fli && ($a1_spliced_orientation ne $a2_spliced_orientation)) { #fli's must have same orient.
	print "failed merge: (a1-fli: $a1_is_fli, a2-fli: $a2_is_fli) and opposite orientations.\n" if $SEE;
	return (0);
    }

    ## Check all overlapping segments to ensure non-conflicting segments.
    my @a1_segments = $a1->get_alignment_segments();
    my @a2_segments = $a2->get_alignment_segments();
    
    
    ## align segment orders between a1 and a2
    my ($starting_a1, $starting_a2);
    for (my $i = 0; $i <= $#a1_segments; $i++) {
	my $a1_seg = $a1_segments[$i];
	my ($a1_lend, $a1_rend) = $a1_seg->get_coords();
	for (my $j = 0; $j <= $#a2_segments; $j++) {
	    my $a2_seg = $a2_segments[$j];
	    my ($a2_lend, $a2_rend) = $a2_seg->get_coords();
	    if (&overlap($a1_lend, $a1_rend, $a2_lend, $a2_rend)) {
		$starting_a1 = $i;
		$starting_a2 = $j;
		last;
	    }
	} 
	if (defined ($starting_a1) && defined ($starting_a2)) {
	    last;
	}
    }
    
    unless (defined ($starting_a1) && defined ($starting_a2)) {
	print "failed merge:  can't align two segments between overlapping alignments.\n" if $::SEE;
	return (0);
    }
    
    unless ($starting_a1 == 0 || $starting_a2 == 0) {
	print "failed merge: segment alignment doesn't begin at either cDNA terminus.\n" if $::SEE;
	return (0);
    }
    
    while ($starting_a1 <= $#a1_segments && $starting_a2 <= $#a2_segments) {
	my $a1_segment = $a1_segments[$starting_a1];
	my $a2_segment = $a2_segments[$starting_a2];
	
	my ($a1_lend, $a1_rend) = $a1_segment->get_coords();
	my ($a2_lend, $a2_rend) = $a2_segment->get_coords();

	if (&overlap($a1_lend, $a1_rend, $a2_lend, $a2_rend)) {
	    ## See if have splice sites, do they exist and are they identical
	    if ($a1_segment->has_left_splice_junction() || $a2_segment->has_left_splice_junction()) {
		if ($a1_segment->has_left_splice_junction() && $a2_segment->has_left_splice_junction() && $a1_lend != $a2_lend) {
		    print "failed merge:\tboth left splice, but unequal coords: L1 ($a1_lend), L2 ($a2_lend)\n" if $::SEE;
		    return (0);
		} elsif ($a1_segment->has_left_splice_junction() && ($a2_lend + $fuzzlength < $a1_lend)) { #alignment extends beyond a splice junction.
		    print "failed merge:\tL1 left splice, L2 ($a2_lend) < L1 ($a1_lend)\n" if $::SEE;
		    return (0);
		} elsif ($a2_segment->has_left_splice_junction() && ($a1_lend + $fuzzlength < $a2_lend)) {
		    print "failed merge:\tL2 left splice, L1 ($a1_lend) < L2 ($a2_lend)\n" if $::SEE; 
		    return (0);
		}
		
	    }
	    if ($a1_segment->has_right_splice_junction() || $a2_segment->has_right_splice_junction()) {
		if ($a1_segment->has_right_splice_junction() && $a2_segment->has_right_splice_junction() && $a1_rend != $a2_rend) {
		    print "failed merge:\tboth right splice, but unequal coords: R1($a1_rend), R2($a2_rend)\n" if $::SEE;
		    return (0);
		} elsif ($a1_segment->has_right_splice_junction() && ($a2_rend - $fuzzlength > $a1_rend)) {
		    print "failed merge:\tR1 right splice, R2($a2_rend) > R1 ($a1_rend)\n" if $::SEE;
		    return (0);
		} elsif ($a2_segment->has_right_splice_junction() && ($a1_rend - $fuzzlength > $a2_rend)) {
		    print "failed merge:\tR2 right splice, R1 ($a1_rend) > R2 ($a2_rend)\n" if $::SEE;
		    return (0);
		}
	    }
	} else {
	    print "failed merge: Two ordered segments don't overlap. ($a1_lend, $a1_rend) , ($a2_lend, $a2_rend)\n" if $::SEE;
	    return (0);
	}
	$starting_a1++;
	$starting_a2++;
	
    }

    ## Passed all tests
    print "Merge tests PASSED.\n" if $::SEE;
    return (1);

}



# private method
# merges two alignment objects together into an assembly.
sub merge_alignments () {
    my ($self, $a1, $a2) = @_;
    print "Merging <" . $a1->get_acc() . ">, <" . $a2->get_acc() . ">\n" if $::SEE;
    my $a1_fli = $a1->is_fli();
    my $a2_fli = $a2->is_fli();
    print "a1_fli: $a1_fli, a2_fli: $a2_fli\n" if $SEE;
    ## Determine fli status for merged product.
    my $merged_fli_status = ($a1_fli || $a2_fli);
    
    my $merged_orientation = $self->determine_merged_orientation($a1, $a2);
    
    #get a1 segment cooridnates;
    my %leftsplicecoords; #preferrentially use splice coords over Fuzzlength extensions.
    my %rightsplicecoords;
    my %a1_coords;
    my @a1_segments = $a1->get_alignment_segments();
    foreach my $seg (@a1_segments) {
    	my ($lend, $rend) = $seg->get_coords();
	$a1_coords{$lend} = $rend;
	if ($seg->has_left_splice_junction()) {
	    $leftsplicecoords{$lend} = 1;
	}
	if ($seg->has_right_splice_junction()) {
	    $rightsplicecoords{$rend} = 1;
	}
    }
    
    # get a2 segment coordinates:
    my %a2_coords;
    my @a2_segments = $a2->get_alignment_segments();
    foreach my $seg (@a2_segments) {
	my ($lend, $rend) = $seg->get_coords();
	$a2_coords{$lend} = $rend;
	if ($seg->has_left_splice_junction()) {
	    $leftsplicecoords{$lend} = 1;
	}
	if ($seg->has_right_splice_junction()) {
	    $rightsplicecoords{$rend} = 1;
	}
    }
    
    my %merged_coords;
    
    #print "Coord dumps:\n" . Dumper (\%a1_coords) . Dumper (\%a2_coords) . "\n";
    
    #print "\n\nBEGIN\n";
    ## merge the overlapping coordinate sets between a1 and a2 alignments.
    foreach my $a1_lend (keys %a1_coords) {
	my $a1_rend = $a1_coords{$a1_lend};
	my ($merged_lend, $merged_rend);
	foreach my $a2_lend (keys %a2_coords) {
	    my $a2_rend = $a2_coords{$a2_lend};
	    if (&overlap($a1_lend, $a1_rend, $a2_lend, $a2_rend)) { #overlap
		## Determine merged lend;
		if ($leftsplicecoords{$a1_lend}) {
		    $merged_lend = $a1_lend;
		} elsif ($leftsplicecoords{$a2_lend}) {
		    $merged_lend = $a2_lend;
		} else {
		    $merged_lend = min ($a1_lend, $a2_lend);
		}
		# Determine merged rend
		if ($rightsplicecoords{$a1_rend}) {
		    $merged_rend = $a1_rend;
		} elsif ($rightsplicecoords{$a2_rend}) {
		    $merged_rend = $a2_rend;
		} else {
		    $merged_rend = max ($a1_rend, $a2_rend);
		}
		last;
	    }
	}
	if ($merged_lend && $merged_rend) {
	    #print "Adding  overlapped \$merged_coords{$merged_lend} = $merged_rend\n";
	    $merged_coords{$merged_lend} = $merged_rend;
	} else { #must not have been any overlap; keep a1 coordset
	    #print "Keeping a1 coords:  \$merged_coords{$a1_lend} = $a1_rend\n";
	    $merged_coords{$a1_lend} = $a1_rend;
	}
    }
    #print "adding unconsumed a2 coords:\n";
    ## add non-overlapping a2-segments
    foreach my $a2_lend (keys %a2_coords) {
	my $overlap = 0;
	my $a2_rend = $a2_coords{$a2_lend};
	foreach my $m_lend (keys %merged_coords) {
	    my $m_rend = $merged_coords{$m_lend};
	    if (&overlap($a2_lend, $a2_rend, $m_lend, $m_rend)) {
		$overlap = 1;
		last;
	    }
	}
	if (!$overlap) {
	    #print "Consuming a2 non-overlapping coords.\n";
	    $merged_coords{$a2_lend} = $a2_rend; #added a2 non-overlapping coordset
	} else {
	    #print "coords overlapped, not consuming.\n";
	}
    }
    print Dumper (\%merged_coords) if $::SEE;
    ## Create a new alignment based on a1 and a2
    my @alignment_segments;
    my $merged_length = 0;
    foreach my $end5 (keys %merged_coords) {
	my $end3 = $merged_coords{$end5};
	my $alignment_seg = new CDNA::Alignment_segment($end5, $end3);
	push (@alignment_segments, $alignment_seg);
	$merged_length += abs ($end3 - $end5) + 1;
    }
    my $new_alignment = new CDNA::CDNA_alignment($merged_length, \@alignment_segments, $self->{sequence_ref});
    my $new_acc = $self->merge_accs($a1->get_acc(), $a2->get_acc());
    $new_alignment->set_acc($new_acc);
    $new_alignment->set_fli_status($merged_fli_status);
    $new_alignment->force_spliced_validation($merged_orientation);
    #print "END.\n\n";
    return ($new_alignment);
}


#private method
# returns the minimum of an array of numerical values.
sub min {
    my @x = @_;
    @x = sort {$a<=>$b} @x;
    my $y = shift @x;
    return ($y);
}


#private method
# returns the maximum of an array of numerical values.
sub max {
    my @x = @_;
    @x = sort {$a<=>$b} @x;
    my $y = pop @x;
    return ($y);
}



#private method
# returns true/false, determines whether two coordinate sets overlap each other.
sub overlap {
    my ($a1_lend, $a1_rend, $a2_lend, $a2_rend) = @_;
    #print "Checking overlap @_\t";
    if ($a2_rend >= $a1_lend && $a2_lend <= $a1_rend) { #overlap
	#print "YES\n";
	return (1);
    } else {
	#print "NO\n";
	return (0);
    }
}





=item toAlignIllustration()

=over 4

B<Description:> illustrates the individual cDNAs to be assembled along with the final products.

B<Parameters:> $max_line_chars(optional)

$max_line_chars is an integer representing the maximum number of characters in a single line of output to the terminal.  The default is  100.

B<Returns:> $alignment_illustration_text

$alignment_illustration_text is a string containing a paragraph of text which illustrates the alignments and assemblies. An example is below:

 --->    <-->  <----->     <--->    <----------------	(+)gi|1199466

 --->    <-->  <----->     <--->    <------------   (+)gi|1209702
                        
---->    <-->  <----	(+)AV827070

---->    <-->  <---	(+)AV828861

---->    <-->  <---      (+)AV830936

 --->    <-->  <-	(+)H36350

ASSEMBLIES: (1)

---->    <-->  <----->     <--->    <---------------- (+) gi|1199466, gi|1209702, AV827070, AV828861, AV830936, H36350




=back

=cut


sub toAlignIllustration () {
    my $self = shift;
    my $max_line_chars = shift;
    $max_line_chars = ($max_line_chars) ? $max_line_chars : 100; #if not specified, 100  chars / line is default.
    
    ## Get minimum coord for relative positioning.
    my @coords;
    my @alignments = @{$self->{incoming_alignments}};
    foreach my $alignment (@alignments) {
	my @c = $alignment->get_coords();
	push (@coords, @c);
    }
    @coords = sort {$a<=>$b} @coords;
    print "coords: @coords\n" if $::SEE;
    my $min_coord = shift @coords;
    my $max_coord = pop @coords;
    my $rel_max = $max_coord - $min_coord;
    my $alignment_text = "";
    ## print each alignment followed by assemblies:
    my $num_alignments = $#alignments + 1;
    $alignment_text .= "Individual Alignments: ($num_alignments)\n";
    my $i = 0;
    foreach my $alignment (@alignments) {
	$alignment_text .= (sprintf ("%3d ", $i)) . $alignment->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
	$i++;
    }
    
    my @assemblies = @{$self->{assemblies}};
    my $num_assemblies = $#assemblies + 1;
    $alignment_text .= "\n\nASSEMBLIES: ($num_assemblies)\n";
    foreach my $assembly (@assemblies) {
	$alignment_text .= "    " . $assembly->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
    }

    return ($alignment_text);
}



=over 4

=item set_fuzzlength()

B<Description:> Sets the fuzzlength parameter.

B<Parameters:> int 

B<Returns:> none.

The fuzzlength is the length allowed to be fuzzy at the terminus of all alignments when compared to overlapping exons containining nearby splice sites.

=back

=cut

sub set_fuzzlength {
    my $self = shift;
    my $fuzzlength = shift;

    $self->{fuzzlength} = $fuzzlength;
}



#private method
# determines if two assemblies have the same accessions.
sub composition_same () {
    my ($name1, $name2) = @_;
    if ($name1 eq $name2) {
	return (1);
    } else {
	return (0);
    }
}

#private method
# creates a new name based on the accessions composing two assemblies.
sub merge_accs () {
    my $self = shift;
    my @names = @_;
    my @nameaccs;
    foreach my $name (@names) {
	my @nameacclist = split (/$DELIMETER/, $name);
	push (@nameaccs, @nameacclist);
    }
    my %unique;
    foreach my $acc (@nameaccs) {
	$unique{$acc} = 1;
    }
    my @accs = sort {$a cmp $b} keys %unique;
    my $combo = join ($DELIMETER, @accs);
    return ($combo);
}

#private method
# checks to see if an assembly already contains all the accessions built into a second assembly.
sub already_contains () {
    ## Checks to see if a1 contains all accessions of a2
    my $self = shift;
    my ($a1, $a2) = @_;
    my $a1_acc = $a1->get_acc();
    my $a2_acc = $a2->get_acc();
    my @a1_accs = split (/$DELIMETER/, $a1_acc);
    my @a2_accs = split (/$DELIMETER/, $a2_acc);
    my %a1_acc_hash;
    foreach my $acc (@a1_accs) {
	$a1_acc_hash{$acc} = 1;
    }
    my %a2_acc_hash;
    foreach my $acc (@a2_accs) {
	$a2_acc_hash{$acc} = 1;
    }
    
    foreach my $acc (keys %a2_acc_hash) {
	unless ($a1_acc_hash{$acc}) {
	    return (0);
	}
    }
    ## if still here, then must have all a2 accs in a1.
    return (1);
}


#private
sub determine_merged_orientation {
    my ($self, $a1, $a2) = @_;
    ## Determine the orientation of the merged product.
    my $num_a1_segments = $a1->get_num_segments();
    my $num_a2_segments = $a2->get_num_segments();
    my $a1_spliced_orientation = $a1->get_spliced_orientation();
    my $a2_spliced_orientation = $a2->get_spliced_orientation();
    print "a1_spliced_orient: $a1_spliced_orientation, a2_spliced_orient: $a2_spliced_orientation\n" if $SEE;
    my $merged_orientation = '+'; #initialize to a default.
    if ($a1_spliced_orientation eq $a2_spliced_orientation) {
	$merged_orientation = $a1_spliced_orientation; 
	print "Same spliced orientation: $merged_orientation\n" if $SEE;
    } elsif ($a1->is_fli() || $a2->is_fli()) {
	$merged_orientation = ($a1->is_fli()) ? $a1_spliced_orientation : $a2_spliced_orientation;
	print "a1 is fli, using a1 orient: $merged_orientation\n" if ($a1->is_fli() && $SEE);
	print "a2 is fli, using a2 orient: $merged_orientation\n" if ($a2->is_fli() && $SEE);
   
    } elsif ($num_a1_segments > 1 || $num_a2_segments > 1) {
	$merged_orientation = ($num_a1_segments > 1) ? $a1_spliced_orientation : $a2_spliced_orientation;
	print "a1 has multiple segments, using a1 orient: $merged_orientation\n" if ($num_a1_segments > 1 && $SEE);
	print "a2 has multiple segments, using a2 orient: $merged_orientation\n" if ($num_a2_segments > 1 && $SEE);
	
    } else {
	print "Can't predict the orientation of the merged alignments ($a1_spliced_orientation vs. $a2_spliced_orientation).  Using default '+'\n" if $SEE;
    }
    return ($merged_orientation);
}



1; #EOM





