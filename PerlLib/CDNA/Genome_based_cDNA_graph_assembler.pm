#!/usr/local/bin/perl

package main;
our $SEE;


package CDNA::Genome_based_cDNA_graph_assembler;

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
use base ("CDNA::Genome_based_cDNA_assembler");

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
    $self->CDNA::Genome_based_cDNA_assembler::_init(@_);
    $self->{show_matrix} = 0;
    $self->{compatibilities} = []; #retains results of compatibility comparisons:
    $self->{incapsulations} = [];
    $self->{LobjsForient} = [];  #non-fli single-segment alignments forced in forward orientation.
    $self->{LobjsRorient} = []; #non-fli single segment alignments forced in reverse orientation.
    $self->{indices_included} = []; #remember what indices are built into assemblies.
}


=item assemble_alignments()

=over 4

B<DESCRIPTION:> assembles a series of cDNA aligmnments into one or more cDNA assemblies using a directed acyclic graph.

B<Parameters:> @alignments

@alignments is an array of CDNA::CDNA_alignment objects

B<Returns:> none.

=back

=cut


sub assemble_alignments {
    my $self = shift;
    my @alignments = @_;
    @alignments = sort {$a->{lend}<=>$b->{lend}} @alignments; #keep in order of lend across genomic sequence to provide a layout.
    $self->{incoming_alignments} = [@alignments];
    
    my $num_alignments = $#alignments + 1;

    ## Initialize globals.
    $self->{incapsulations} = []; #initialize.
    $self->{compatibilities} = {}; #initialize.
    $self->{indices_included} = []; #initialize
    
    $self->{assemblies} = []; #initialize
    
    #precompute compatible alignment pairs and full encapsulations.
    $self->determine_compatibilities_and_encapsulations();
    
    ######################################
    # Do Forward Scan Dynamic Programming.
    my $LobjsForient = $self->{LobjsForient} = [];
    my $LobjsRorient = $self->{LobjsRorient} = [];
    $self->force_flexorient('+');
    $self->do_full_Fscan($LobjsForient);
    $self->force_flexorient('-');
    $self->do_full_Fscan($LobjsRorient);

    &Describe_containment($LobjsForient, $LobjsRorient) if $SEE;
    
    ## Get the highest scoring Lobj
    
    my $struct = $self->get_top_alignment_indices();
    
    ## Create the highest scoring assembly.
    my $assembly = $self->create_assembly($struct);
    my $num_cdnas_included = $struct->{num_cdnas};
    
    if ( $num_cdnas_included == $num_alignments) { return();}
    
   
    #######################################
    ## Do Reverse Scan Dynamic Programming (if cDNAs exist and not included in assembly above).
    $self->force_flexorient('+');
    $self->do_full_Rscan($LobjsForient);
    $self->force_flexorient('-');
    $self->do_full_Rscan($LobjsRorient);
    ## Create assemblies for each cDNA not included yet.
    my @asmbl_sets;
    my %seen;
    for (my $i = 0; $i < $num_alignments; $i++) {
	unless ($self->{indices_included}->[$i]) {
	    print "Not Included: $i, Doing Back trace starting at $i.\n" if $SEE;
	    my $struct = $self->get_top_alignment_indices_starting_index($i);
	    push (@asmbl_sets, $struct);
	}
    }
    @asmbl_sets = reverse sort {$a->{num_cdnas}<=>$b->{num_cdnas}} @asmbl_sets;
    my $all_included_flag = 0; 
    while ((!$all_included_flag) && @asmbl_sets) {
	my $asmbl_set = shift @asmbl_sets;
	$all_included_flag = $self->test_and_add_asmbl($asmbl_set);
    }
    unless ($all_included_flag) {
	die "FATAL: All assemblies not included!!\n";
    }
    
}

# private.
sub do_full_Fscan {
    my $self = shift;
    my $Lobjs = shift;

    print "\n# Doing full Fscan\n" if $SEE;
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
   
    for (my $i = 0; $i < $num_alignments; $i++) {
	my $Lobj = Lobject->new($i, $self);
	if ($i != 0) { #first cDNA, base case.
	    ## Must compare to previous alignments:
	    my $top_score = 0;
	    my $top_scoring_index = -1;
	    my $i_alignment = $alignments_aref->[$i];
	    for (my $j = $i - 1; $j >= 0; $j--) {
		my $j_alignment = $alignments_aref->[$j];
		my $prevLobj = $Lobjs->[$j];
		print "Comparing $i to $j\n" if $SEE;
		my $compatibility = $self->{compatibilities}->{$i}->{$j};
		my $j_contains_i = $prevLobj->{contained_cdna_indices}->{$i};
		my $i_contains_j = $Lobj->{contained_cdna_indices}->{$j};
		print "\tcompat: $compatibility\tj_contains_i: $j_contains_i\n" if $SEE;
		if ($compatibility && (!($j_contains_i||$i_contains_j))) {
		    unless (&runtime_compatible($i_alignment, $j_alignment)) {next;}
		    my $curr_Lscore = $Lobjs->[$j]->{LscoreF};
		    my $Cscore = &get_Cscore($Lobj, $prevLobj);
		    my $curr_total_score = $curr_Lscore + $Cscore;
		    
		    print "\tFSCAN score: ($i,$j) = $curr_total_score\n" if $SEE;
		    if ($curr_total_score > $top_score) {
			$top_scoring_index = $j;
			$top_score = $curr_total_score;
			print "\t*making topscore.\n" if $SEE;
			
		    }
		}
	    }
	    if ($top_scoring_index > -1) {
		my $topLobj = $Lobjs->[$top_scoring_index];
		$Lobj->{fromLptr} = $topLobj;
		$Lobj->{LscoreF} = $top_score;
	
	    }
	    	    
	}
	$Lobjs->[$i] = $Lobj;
    }
}

# private.
sub do_full_Rscan {
    my $self = shift;
    my $Lobjs = shift;


    print "\n# Doing full Rscan\n" if $SEE;
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    
    for (my $i = $num_alignments - 2; $i >= 0; $i--) {
	my $Lobj = $Lobjs->[$i];
	my $i_alignment = $alignments_aref->[$i];
	## Must compare to previous alignments:
	my $top_score = 0;
	my $top_scoring_index = -1;
	for (my $j = $i + 1; $j < $num_alignments; $j++) {
	    my $j_alignment = $alignments_aref->[$j];
	    my $nextLobj = $Lobjs->[$j];
	    print "Comparing $i to $j\n" if $SEE;
	    my $compatibility = $self->{compatibilities}->{$i}->{$j};
	    my $j_contains_i = $nextLobj->{contained_cdna_indices}->{$i};
	    my $i_contains_j = $Lobj->{contained_cdna_indices}->{$j};
	    print "\tcompat: $compatibility\tj_contains_i: $j_contains_i\n" if $SEE;
	    if ($compatibility && (!($j_contains_i||$i_contains_j))) {
		unless (&runtime_compatible($i_alignment, $j_alignment)) { next;}
		my $curr_Lscore = $Lobjs->[$j]->{LscoreR};
		my $Cscore = &get_Cscore($Lobj, $nextLobj);
		my $curr_total_score = $curr_Lscore + $Cscore;
		
		print "\tRSCAN score: ($i,$j) =  $curr_total_score\n" if $SEE;
		if ($curr_total_score > $top_score) {
		    $top_scoring_index = $j;
		    $top_score = $curr_total_score;
		    print "\t*making topscore.\n" if $SEE;
		    
		}
	    }
	}
	if ($top_scoring_index > -1) {
	    my $topLobj = $Lobjs->[$top_scoring_index];
	    $Lobj->{toLptr} = $topLobj;
	    $Lobj->{LscoreR} = $top_score;
	}
	
    }
}


# private.
sub get_Cscore {
    ## Cscore is the number of cdnas contained in a but not in b including a itself
    my ($a, $b) = @_;
    my $score = 0;
    my $a_container_ref = $a->{contained_cdna_indices};
    my $b_container_ref = $b->{contained_cdna_indices};
    foreach my $a_index (keys %$a_container_ref) {
	unless ($b_container_ref->{$a_index}) {
	    $score++;
	}
    }
    return ($score);
}

    
=over 4

=item encapsulates()

B<Description:> Returns true (1) if alignment A encapsulates the span of alignment B

B<Parameters:> ($alignment_A, $alignment_B)

Alignments A and B are of type CDNA::CDNA_alignment

B<Returns:> [1|0]

=back

=cut


sub encapsulates {
    my $self = shift;
    my ($alignmentA, $alignmentB) = @_;
    my ($alend, $arend) = $alignmentA->get_coords();
    my ($blend, $brend) = $alignmentB->get_coords();
    
    if ($blend >= $alend && $brend <= $arend) {
	return (1);
    } else {
	return (0);
    }
}


sub determine_compatibilities_and_encapsulations {
    my $self = shift;
    
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    
   
    #initialize incapsulations list
    for (my $i = 0; $i < $num_alignments; $i++) {
        $self->{incapsulations}->[$i] = [];
    }
    
    #compare each downstream alignment to the upstream alignment for compatiblity and encapsulation:
    for (my $i = 0; $i < $num_alignments; $i++) {
        for (my $j = $i + 1; $j < $num_alignments; $j++) {
            
            my $alignment_i = $alignments_aref->[$i];
            my $alignment_j = $alignments_aref->[$j];
            
            
            my $can_merge = 0;
            if ($self->can_merge($alignment_i, $alignment_j)) {
                $can_merge = 1;
                # set compatibility flags.
                $self->{compatibilities}->{$j}->{$i} = 1;
                $self->{compatibilities}->{$i}->{$j} = 1;
            }
            
            my $align_i_single_seg = ($alignment_i->get_num_segments() == 1) ? 1:0;
            my $align_j_single_seg = ($alignment_j->get_num_segments() == 1) ? 1:0;
            
            
            ## Analyze compatible alignments for containment.  Single segment alignments may conflict if in opposite orientations and fli-status, so should check these incompatible alignments for containment as well to avoid multiple identical assemblies from being generated.
            
            if ($can_merge || ($align_i_single_seg && $align_j_single_seg)) {
                
                ## Check for encapsulation:
                if ($self->encapsulates($alignment_j, $alignment_i)) {
                    push (@{$self->{incapsulations}->[$j]}, $i);
                    print "$j contains $i\n" if $SEE;
                } 
                if ($self->encapsulates($alignment_i, $alignment_j)) {
                    push (@{$self->{incapsulations}->[$i]}, $j);
                    print "$i contains $j\n" if $SEE;
                }
                
            }
        }
    }
}

sub force_flexorient {
    my $self = shift;
    my $orient = shift;
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    ## Fix orientations of fli and multi-segment alignments:
    for (my $i = 0; $i < $num_alignments; $i++) {
        my $alignment = $alignments_aref->[$i];
        my $num_segments = $alignment->get_num_segments();
        my $spliced_orientation = $alignment->get_spliced_orientation();
        
        if ($spliced_orientation =~ /^[+-]$/) { # set specifically
            $alignment->{fixed_orient} = $spliced_orientation;  ## adding tag, using only in this module.
        } else {
            print "Setting $i to flex orient: $orient.\n" if $SEE;
            $alignment->{fixed_orient} = $orient;
        }
    }
    
}


sub back_trace {
    my $self = shift;
    my $top_scoring_index = shift;
    my $Lobjs = shift;

    print "Back_trace: " if $SEE;
    my $Lobj = $Lobjs->[$top_scoring_index];
    ## Traverse the fromLptrs
    
    my %unique_indices;
    my @trace_indices;    
    while ($Lobj != 0) {
	my $index = $Lobj->{myIndex};
	print "$index " if $SEE;
	push (@trace_indices, $index);
	# get composed accs:
	foreach my $contained_index ($Lobj->get_contained_indices()) {
	    $unique_indices{$contained_index} = 1;
	}
	$Lobj = $Lobj->{fromLptr};
    }
    print "\n" if $SEE;
    my @unique = keys %unique_indices;
    my $num_cdnas_included = $#unique + 1; 
    my $struct = {trace_indices=>\@trace_indices,
		  num_cdnas=>$num_cdnas_included,
		  all_indices=>\@unique};
    return ($struct);
}


sub forward_trace {
    my $self = shift;
    my $index = shift;
    my $Lobjs = shift;

    print "Forward_trace: " if $SEE;
    my $Lobj = $Lobjs->[$index];
    ## Traverse the fromLptrs
    my @trace_indices;
    my %unique_indices;
    while ($Lobj != 0) {
	my $index = $Lobj->{myIndex};
	print " $index" if $SEE;
	push (@trace_indices, $index);
	# get composed accs:
	foreach my $contained_index ($Lobj->get_contained_indices()) {
	    $unique_indices{$contained_index} = 1;
	}
	$Lobj = $Lobj->{toLptr};
    }
    

    print "\n" if $SEE;
    my @unique = keys %unique_indices;
    my $num_cdnas_included = $#unique + 1; 
    my $struct = {trace_indices=>\@trace_indices,
		  num_cdnas=>$num_cdnas_included,
		  all_indices=>\@unique};
    
    return ($struct);
}


sub create_assembly {
    my $self = shift;
    my $struct = shift;
    my $indices_aref = $struct->{trace_indices};
    my $alignments_aref = $self->{incoming_alignments};
    my @indices = sort {$a<=>$b} @$indices_aref;
    my $all_indices_ref = $struct->{all_indices};
    my @all_indices = sort {$a<=>$b} @$all_indices_ref;
    print "Creating assembly: @indices\n" if $SEE;
   
    my $first_index = shift @indices;
    my $alignment = $alignments_aref->[$first_index];
    $self->{indices_included}->[$first_index] = 1;
    print "Index: $first_index\t" . $alignment->toToken() . "\n" if $SEE;
    my $assembly = $alignment->clone();
    
    my @alignments_to_assemble;
    # first index already included in assembly; ...now, include the others.
    foreach my $index (@indices) {
	my $alignment = $alignments_aref->[$index];
	print "Index: $index\t" . $alignment->toToken() . "\n" if $SEE;
	$assembly = $self->merge_alignments($assembly, $alignment);
    }
    print "assembly contains: @all_indices\n" if $SEE;
    my @all_accs;
    
    foreach my $index (@all_indices) {
	$self->{indices_included}->[$index] = 1;
	my $acc = $alignments_aref->[$index]->get_acc();
	push (@all_accs, $acc);
    }
    my $acc = $self->merge_accs(@all_accs);
    $assembly->set_acc($acc);
        
    push (@{$self->{assemblies}}, $assembly);;
    
    return ($assembly);
}

sub runtime_compatible {
    my ($a, $b) = @_;
    my $a_fixed_orient = $a->{fixed_orient};
    my $b_fixed_orient = $b->{fixed_orient};
    print "a($a_fixed_orient) vs. b($b_fixed_orient)\n" if $SEE;
    if ($a_fixed_orient eq $b_fixed_orient) {
	print "runtime compatible.\n" if $SEE;
	return (1);
    } else {
	print "NOT runtime compatible\n" if $SEE;
	return (0);
    }
}

	
sub get_top_alignment_indices {
    my $self = shift;
    my $LobjsForient = $self->{LobjsForient};
    my $LobjsRorient = $self->{LobjsRorient};
    
    ## Process Forient (singletons forced Forward orient)
    my $top_Forient_index = $self->top_scoring_index_from_Lobjs($LobjsForient, 'LscoreF');
    my $ForientStruct = $self->back_trace($top_Forient_index, $LobjsForient);
    
    ## Process Rorient (singletons forced Reverse orient)
    my $top_Rorient_index = $self->top_scoring_index_from_Lobjs($LobjsRorient, 'LscoreF');
    my $RorientStruct = $self->back_trace($top_Rorient_index, $LobjsRorient);
    if ($ForientStruct->{num_cdnas} >= $RorientStruct->{num_cdnas}) {
	return ($ForientStruct);
    } else {
	return ($RorientStruct);
    }
}

sub get_top_alignment_indices_starting_index {
    my $self = shift;
    my $index = shift;
    my $LobjsForient = $self->{LobjsForient};
    my $LobjsRorient = $self->{LobjsRorient};
    my $LobjsForient_num_contained = $LobjsForient->[$index]->get_contained_indices();
    my $LobjsRorient_num_contained = $LobjsRorient->[$index]->get_contained_indices();
    
    ## Process Forient (singletons forced forward orient)
    my $ForientBtraceStruct = $self->back_trace($index, $LobjsForient);
    my $ForientFtraceStruct = $self->forward_trace($index, $LobjsForient);
    my $ForientNumCdnas = $ForientBtraceStruct->{num_cdnas} + $ForientFtraceStruct->{num_cdnas} -$LobjsForient_num_contained;
        
    ## Process Rorient (singletons forced reverse orient)
    my $RorientBtraceStruct = $self->back_trace($index, $LobjsRorient);
    my $RorientFtraceStruct = $self->forward_trace($index, $LobjsRorient);
    my $RorientNumCdnas = $RorientBtraceStruct->{num_cdnas} + $RorientFtraceStruct->{num_cdnas} - $LobjsRorient_num_contained;

    my ($Ftrace, $Btrace,$total_cdnas);
    if ($ForientNumCdnas > $RorientNumCdnas) {
	($Ftrace, $Btrace) = ($ForientFtraceStruct, $ForientBtraceStruct);
	$total_cdnas = $ForientNumCdnas;
    } else {
	($Ftrace, $Btrace) = ($RorientFtraceStruct, $RorientBtraceStruct);
	$total_cdnas = $RorientNumCdnas;
    }
    
    ## Create new struct
    my @unique_trace = &unique_entries(@{$Ftrace->{trace_indices}}, @{$Btrace->{trace_indices}});
    my @all_indices = &unique_entries (@{$Ftrace->{all_indices}}, @{$Btrace->{all_indices}});
    if ($#all_indices +1 != $total_cdnas) {
	die "Error: total number of alignments not calculated correctly.\n";
    }
    print "FnB trace from index[$index] yields assembly containing indices [@all_indices]\n" if $SEE;
    my $struct = {trace_indices=>\@unique_trace,
		  num_cdnas=>$total_cdnas,
		  all_indices=>\@all_indices};
    return ($struct);
}
 



sub unique_entries {
    my @x = @_;
    my %z;
    foreach my $y (@x) {
	$z{$y}=1;
    }
    return (keys %z);
}


sub top_scoring_index_from_Lobjs {
    my $self = shift;
    my $Lobjs = shift;
    my $score_type = shift;
    
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    my $top_scoring_index = -1;
    my $top_score = -1;

    for (my $i = 0; $i < $num_alignments; $i++) {
	my $Lscore = $Lobjs->[$i]->{$score_type};
	if ($Lscore > $top_score) {
	    $top_scoring_index = $i;
	    $top_score = $Lscore;
	}
    }
    return ($top_scoring_index);
}


sub test_and_add_asmbl {
    my $self = shift;
    my $asmbl_set = shift;
    my $all_indices_ref = $asmbl_set->{all_indices};
    my $indices_included_ref = $self->{indices_included};
## If the asmbl_set contains a cDNA missing from an assembly inclusion, then add it.
    
    
    my $num_alignments = $#{$self->{incoming_alignments}} + 1;
    my $contains_unincorporated = 0;
    foreach my $index (@$all_indices_ref) { 
	if (! $indices_included_ref->[$index]) {
	    $contains_unincorporated = 1;
	    last;
	}
    }
    if ($contains_unincorporated) {
	$self->create_assembly($asmbl_set);
    }
    
    ## Test to see if all cDNAs are accounted for in assemblies.
    my $all_included_flag = 1;
    
    for (my $i =0; $i < $num_alignments; $i++) {
	if (! $indices_included_ref->[$i]) {
	    print "$i not included!\n" if $SEE;
	    $all_included_flag = 0;
	}
    }
    return ($all_included_flag);
}

sub Describe_containment {
    my ($LobjsForient, $LobjsRorient) = @_;
    for (my $i = 0; $i <= $#$LobjsForient; $i++) {
	print "Fixed Forient: index $i contains: " . join (" ", $LobjsForient->[$i]->get_contained_indices()) . "\n";
	print "Fixed Rorient: index $i contains: " . join (" ", $LobjsRorient->[$i]->get_contained_indices()) . "\n";
    }
    print "Lscores:\n";
    for (my $i = 0; $i <= $#$LobjsForient; $i++) {
	print "$i: F-forced singleton orient Lscores: F: " .   $LobjsForient->[$i]->{LscoreF} . " R: " . $LobjsForient->[$i]->{LscoreR} . "\n";
	print "$i: R-forced singleton orient Lscores: F: " .   $LobjsRorient->[$i]->{LscoreF} . " R: " . $LobjsRorient->[$i]->{LscoreR} . "\n";
    }
}





#######################
package Lobject;

sub new {
    my $packagename = shift;
    my $index = shift;
    my $assembler_ref = shift;
    my $self = {
	contained_cdna_indices =>{ $index => 1  #include as containing itself.
				   }, 
	    myIndex => $index,  #remember Lobj position.
	    LscoreF => 1, #Forward scan Lscore.
	    LscoreR => 1, #Reverse scan Lscore.
	    toLPtr => 0,  # used in Rscan for forward tracking.
	    fromLPtr => 0 # used in Fscan for backtracking
	    };
    
    ## populate list of contained entries.
    my @contained_indices = @{$assembler_ref->{incapsulations}->[$index]};
    foreach my $contained_index (@contained_indices) {
	$self->{contained_cdna_indices}->{$contained_index} = 1;
	$self->{LscoreF}++;
	$self->{LscoreR}++;
    }
    
    bless ($self, $packagename);
    return ($self);
}

sub get_contained_indices {
    my $self = shift;
    return (keys %{$self->{contained_cdna_indices}});
}


1; #EOM





