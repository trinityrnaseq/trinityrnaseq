#!/usr/local/bin/perl

##########################
### Class CDNA_alignment
##########################
package main;
our $SEE;

=head1 NAME

CDNA::CDNA_alignment

=cut

=head1 DESCRIPTION

This module provides an object specification for storing and manipulating cDNA alignments.  The alignment coordinates along the genomic sequence are given along with a genomic sequence.  Splice sites are validated and the spliced orientation is determined.

=cut



package CDNA::CDNA_alignment;
use strict;
use Exons_to_geneobj;
use Gene_obj;
use CDNA::Alignment_segment;
use Carp qw (cluck confess croak);
use GFF_maker;

## Global Vars
our $ALLOW_ATAC_splice_pairs = 1; # ON by default; turn it off if you prefer.

##


=over 4

=item new()

B<Description:> Instantiates a new CDNA::CDNA_alignment object. 

B<Parameters:> $cdna_length, $Alignment_obj_aref, $sequence_sref

B<$cdna_length> is the length of the complete cDNA sequence.

B<$Alignment_obj_aref> is a reference to a list of CDNA::Alignment_segment objects like so:

$Alignment_obj_aref = \@Alignment_obj_list

B<$sequence_sref> should be a reference to the string containing the genomic sequence like so:
    
$sequence = "gatc.....";

$sequence_sref = \$sequence;

B<Returns:> $obj_ref

returns a reference to a CDNA_alignment object.  This object is validated requiring consensus splice sites, and the spliced orientation of the cDNA is determined based on the validating orientation.

=back

=cut


sub new {
    my $packagename = shift;
    my ($cdna_length, $Alignment_obj_aref, $sequence_ref) = @_; 
    unless (@$Alignment_obj_aref) {
        die "No alignment segments available to create alignment object.\n";
    }
    
    my $self = {
        acc => undef(), # cDNA accession.
        title =>undef(), # com_name, title, header, whatever you want to call the cdna.
        genomic_seq => $sequence_ref,
        genome_acc => undef,
        alignment_segs=>$Alignment_obj_aref,
        orientation => undef(),
        lend=>undef(),
        rend=>undef(),
        length=>0, # alignment span length.
        num_aligned_nts => 0, #number of cDNA nucleotides found in alignment
        cdna_length => $cdna_length, # the length of the cDNA sequence.
        avg_per_id => 0, # the average percent identity of this alignment.
        error_flag=>0, #default w/o errors. Used to indicate non-consensus splice sites.
        spliced_orientation => '?', # [+-?] depending on validating consenus splice sites; relative to genomic sequence orientation.
        num_segments=>0,
        percent_cdna_aligned => 0,   # provides the percentage of the cDNA length that is in the alignment.
        is_fli => 0 #indicates whether the cDNA is a full-length insert (ie. expected to be complete).
        };
    bless ($self, $packagename);
    
    ## Convert alignment to object form:
    $self->process_alignment();
    
    if ($self->{num_segments} > 1 && (ref $sequence_ref)) {
        $self->identify_splice_junctions($sequence_ref);
    }
    return ($self);
}

sub process_alignment {
    my $self = shift;
    $self->determine_alignment_attributes();
    my @alignment_segments = $self->get_alignment_segments();
    for (my $i = 0; $i <= $#alignment_segments; $i++) {
        my $segment = $alignment_segments[$i];
        if ($#alignment_segments == 0) { #single segment in alignment
            $segment->set_type("single");
        } elsif ($i == 0) {
            $segment->set_type("first");
        } elsif ($i == $#alignment_segments) {
            $segment->set_type("last");
        } else { #must be internal
            $segment->set_type("internal");
        }
    }
    $self->set_num_segments($#alignment_segments + 1);
    $self->verify_contiguity();
    
}

sub identify_splice_junctions {
    my $self = shift;
    my $sequence_ref = shift;
    my $orientation = $self->get_orientation();
    
    my $opposite_orientation = ($orientation eq '+') ? '-' : '+';
    
    ## Since some cDNAs are supplied in the opposite orientation, the splice sites may be on the reverse strand.
    ## In this case, we must change the orientation
    my $error_flag = 0;
    my $error_text = "";
    my $validating_orientation = $orientation; #initialize.
    my %error_lengths; #used to find best orientation (least errors)
    foreach my $orient ($orientation, $opposite_orientation) {
        $error_text = $self->validate_splice_junctions($orient, $sequence_ref);
        $error_lengths{$orient} = length $error_text;
        if ($error_text) {
            print "ERRORS in Splice sites given orientation $orient:\n$error_text\n" if $::SEE;
            unless ($error_flag) { print "Trying other strand might help?\n" if $::SEE;}
            $error_text = ""; #rest
            $error_flag = 1;
        } else {
            $error_flag = 0;
            $validating_orientation = $orient;
            last;
        }
    }
    if ($error_flag) {
        print "Sorry, still contains problematic splice sites:\n$error_text\n" if $::SEE;
        print "Setting error flag for this alignment\n" if $::SEE;
        $self->set_error_flag("Splice site validations failed");
        print "ERROR: SETTING ERROR_FLAG\n" if $::SEE;
        ## revalidate splice sites using orientation that generated the least errors:
        my @orients = sort {$error_lengths{$a}<=>$error_lengths{$b}} ('+', '-');
        my $best_orient = shift @orients;
        $self->validate_splice_junctions($best_orient, $sequence_ref); #required for appropriate token printing.
    } else {
        $self->set_spliced_orientation($validating_orientation);
    }
}

sub validate_splice_junctions {
    my $self = shift;
    my ($orient) = shift;
    my $sequence_ref = shift;
    my $errors = "";
    my (@splice_boundary_pairs) = &get_consensus_splice_sites($orient);
    my @segments = $self->get_alignment_segments();
    
    my $num_segments = scalar (@segments);
    if ($num_segments == 1) {
        ## no introns
        die "Error, trying to validate splice junctions for single segment alignment! ";
    }
    
    ## analyze introns
    for (my $i = 1; $i <= $#segments; $i++) {
        my ($prev_segment, $curr_segment) = ($segments[$i-1], $segments[$i]);
        
        my ($prev_lend, $prev_rend) = $prev_segment->get_coords();
        my ($curr_lend, $curr_rend) = $curr_segment->get_coords();

        my $splice_chars_left = uc substr($$sequence_ref, $prev_rend, 2);
        my $splice_chars_right = uc substr($$sequence_ref, $curr_lend -2 -1, 2);
        
        $prev_segment->set_right_splice_site_chars($splice_chars_left);
        $curr_segment->set_left_splice_site_chars($splice_chars_right);
        
        ## check boundaries:
        my $splice_pair_OK = 0;
        
      CONSENSUS_PAIR:
        foreach my $consensus_pair (@splice_boundary_pairs) {
            my ($left_chars, $right_chars) = @$consensus_pair;
            
            if ($left_chars eq $splice_chars_left && $right_chars eq $splice_chars_right) {
                # found consensus
                ## further validate any AT-AC for donor site extended consensus
                my ($intron_lend, $intron_rend) = ($prev_rend + 1, $curr_lend - 1);
                if ($ALLOW_ATAC_splice_pairs) {
                    if ( ( $orient eq '+' && $left_chars eq 'AT')
                         ||
                         ($orient eq '-' && $right_chars eq 'AT') ) {
                        unless (&_validates_AT_AC_donor_extended_consensus($orient, $intron_lend, $intron_rend, $sequence_ref)) {
                            next CONSENSUS_PAIR;
                        }
                    }
                }
                
                ## got consensus splice pair
                $splice_pair_OK = 1;
                ## left and right local are relative to intron
                # in methods, they reference segment
                $prev_segment->set_right_splice_junction(1);
                $curr_segment->set_left_splice_junction(1);
                last;
            }
        }
        unless ($splice_pair_OK) {
            
            $errors .= "nonconsensus splice pair [$splice_chars_left-$splice_chars_right]\n";
            
            ## validate boundaries separately.
            ## this is useful if only one site is nonconsensus
            foreach my $consensus_pair (@splice_boundary_pairs) {
                my ($left_chars, $right_chars) = @$consensus_pair;
                if ($left_chars eq $splice_chars_left) {
                    $prev_segment->set_right_splice_junction(1);
                }
                if ($right_chars eq $splice_chars_right) {
                    $curr_segment->set_left_splice_junction(1);
                }
            }
            
            ## make nonconsensus lower case 
            unless ($prev_segment->has_right_splice_junction()) {
                $prev_segment->set_right_splice_site_chars( lc $splice_chars_left);
            }
            unless ($curr_segment->has_left_splice_junction()) {
                $curr_segment->set_left_splice_site_chars( lc $splice_chars_right);
            }
            
        }
    }
    
    return ($errors);
}


sub verify_contiguity {
    my $self = shift;
    my $num_segments = $self->get_num_segments();
    my $orient = $self->get_orientation();
    if ($num_segments > 1) {
        my @alignment_segments = $self->get_alignment_segments();
        for (my $i = 1; $i <= $#alignment_segments; $i++) {
            my $prev_seg = $alignment_segments[$i-1];
            my $curr_seg = $alignment_segments[$i];
            my $error_flag = 0;
            my $diff = $curr_seg->{mlend} - $prev_seg->{mrend};
            if ($orient eq '+' && $diff != 1) {
                $error_flag = 1;
            }elsif ($orient eq '-' && $diff != -1) {
                $error_flag = 1;
            }
            if ($error_flag) {
                $self->set_error_flag("Incontiguous alignment");
                return();
            }
        }
    }
}



sub get_consensus_splice_sites () {
    
    my $orientation = shift;
    
    my @pairs;

    if ($orientation eq '+') {
        
        ## Forward pairs
        # GT-AG
        # GC-AG
        # AT-AC
        
        push (@pairs, ['GT', 'AG'], ['GC', 'AG']);
        
        if ($ALLOW_ATAC_splice_pairs) {
            push (@pairs, ['AT', 'AC']);
        }
              
    }

    else {
        ## Rev Comp of above:
        # CT-AC
        # CT-GC
        # GT-AT
        
        push (@pairs, ['CT', 'AC'], ['CT', 'GC']);
        
        if ($ALLOW_ATAC_splice_pairs) {
            push (@pairs, ['GT', 'AT']);
        }
    }


    return (@pairs);
}


####
sub _validates_AT_AC_donor_extended_consensus {
    my ($orient, $intron_lend, $intron_rend, $sequence_ref) = @_;

    my $forward_consensus = 'ATATCC';
    my $reverse_consensus = 'GGATAT';
        
    if ($orient eq '+') {
        my $long_donor_seq = uc substr($$sequence_ref, $intron_lend - 1, 6);
        if ($long_donor_seq eq $forward_consensus) { 
            return (1);
        }
    }
    elsif ($orient eq '-') {
        my $long_donor_seq = uc substr($$sequence_ref, $intron_rend -6, 6);
        if ($long_donor_seq eq $reverse_consensus) {
            return (1);
        }
    }

    ## got here, didn't fit long consensus
    return (0);
}





sub set_orientation {
    my $self = shift;
    my $orientation = shift;
    $self->{orientation} = $orientation;
}


=over 4

=item get_aligned_orientation()

B<Description:> Provides the orientation of the incoming cDNA alignment.

B<Parameters:> none.

B<Returns:> [+-] 

=back

=cut


sub get_aligned_orientation {
    my $self = shift;
    return ($self->{orientation});
}

sub get_orientation { # deprecated in favor of get_aligned_orientation()
    my $self = shift;
    return ($self->get_aligned_orientation());
}


sub set_title { 
    my $self = shift;
    my $title = shift;
    $self->{title} = $title;
}




=over 4

=item get_title()

B<Description:> Provides the title for the cDNA, generally the header from a fasta file

B<Parameters:> none.

B<Returns:> string or undef 

=back

=cut




sub get_title {
    my $self = shift;
    return ($self->{title});
}



sub set_coords {
    my $self = shift;
    my ($lend, $rend) = @_;
    $self->{lend} = $lend;
    $self->{rend} = $rend;
}


=over 4

=item get_coords()

B<Description:> Provides the coordinate span for the alignment along the genomic sequence.  Orientation is not implied, coordinates always provided relevant to the forward orientation.  Use the get_orientation method to determine the strand.

B<Parameters:> none.

B<Returns:> $lend, $rend

$lend, $rend are integer values indicating the beginning and end coordinates of the alignment on the genomic sequence.  The genomic sequence begins at position 1.

=back

=cut



sub get_coords {
    my $self = shift;
    return ($self->{lend}, $self->{rend});
}





=over 4

=item get_mcoords()

B<Description:> returns the cDNA sequence coordinates corresponding to the lend and rend returned by get_coords()

    ie. my ($lend, $rend) = $alignment->get_coords();  // always in forward orientation (lend < rend)

    my ($mlend, $mrend) = $alignment->get_mcoords(); 

    $mlend corresponds to $lend

    $mrend corresponds to $rend


B<Parameters:>none 

B<Returns:> ($mlend, $mrend)

=back

=cut


sub get_mcoords {
    my $self = shift;

    my @alignment_segments = $self->get_alignment_segments();
    
    my $leftmost_seg = shift @alignment_segments;
    my $rightmost_seg = pop @alignment_segments;
    
    unless ($rightmost_seg) {
        $rightmost_seg = $leftmost_seg; #must only be one in which case left = right
    }
    
    my ($mlend, $whatever1) = $leftmost_seg->get_mcoords();
    my ($whatever2, $mrend) = $rightmost_seg->get_mcoords();
    
    return ($mlend, $mrend);
}




=over 4

=item get_intron_coords()

B<Description:> returns list of intron coordinates

B<Parameters:> none

B<Returns:> ( [intron_lend, intron_rend], [intron_lend, intron_rend], ...)


    lend always less than rend

=back

=cut


sub get_intron_coords {
    my $self = shift;

    my @segments = $self->get_alignment_segments();
    
    my @seg_coords;
    foreach my $segment (@segments) {
        my ($exon_lend, $exon_rend) = sort {$a<=>$b} $segment->get_coords();
        
        push (@seg_coords, [$exon_lend, $exon_rend]);
    }

    my @intron_coords;
    
    if (scalar (@seg_coords) >= 2) {
        ## actually have introns
        @seg_coords = sort {$a->[0]<=>$b->[0]} @seg_coords;
        
        my $curr_seg = shift @seg_coords;
        while (@seg_coords) {
            my $next_seg = shift @seg_coords;

            my ($curr_lend, $curr_rend) = @$curr_seg;

            my ($next_lend, $next_rend) = @$next_seg;
            
            my ($intron_lend, $intron_rend) = ($curr_rend + 1, $next_lend - 1);
            
            push (@intron_coords, [$intron_lend, $intron_rend]);

            $curr_seg = $next_seg;
        }
    }


    return (@intron_coords);
}


sub determine_alignment_attributes {
    my $self = shift;
    my @alignment_segments = $self->get_alignment_segments();
    my @coords;
    my $orientation;
    my $num_nts_matched = 0;
    my $per_id_x_length = 0;
    foreach my $segment (@alignment_segments) {
        my ($lend, $rend) = $segment->get_coords();
        my ($mlend, $mrend) = $segment->get_mcoords();
        my $orient = $segment->get_orientation();
        if (!$orientation && $orient =~ /[+-]/) {
            $orientation = $orient;
        }
        my $seg_length;
        if ($mrend && $mlend) {
            $seg_length = abs($mrend - $mlend) + 1;
        } else {
            $seg_length = abs ($rend - $lend) + 1;
        }
        $num_nts_matched += $seg_length;
        my $per_id = $segment->get_per_id();
        $per_id_x_length += $per_id * $seg_length;
        push (@coords, $lend, $rend);
    }
    $self->set_orientation($orientation);
    foreach my $segment (@alignment_segments) {
        $segment->set_orientation($orientation);
    }
    @coords = sort {$a<=>$b} @coords;
    my $lend = shift @coords;
    my $rend = pop @coords;
    $self->set_coords($lend, $rend);
    $self->{length} = abs ($rend - $lend) + 1;
    $self->{num_aligned_nts} = $num_nts_matched;
    $self->{avg_per_id} = $per_id_x_length / $num_nts_matched;
    if ($self->{avg_per_id} > 100) { die "Error, can't have average per_id > 100%\n";}
    if ($self->{cdna_length} > 0) {
        $self->{percent_cdna_aligned} = $num_nts_matched / $self->{cdna_length} * 100;
    }     
}


=over 4
    
=item get_alignment_segments()
    
B<Description:> Returns the alignment segments which comprise an alignment object, ordered by left coordinate position.

B<Parameters:> none

B<Returns:> @alignment_segments

@alignment_segments is a list of CDNA::Alignment_segment objects (see below).

=back

=cut

sub get_alignment_segments() {
    my $self = shift;
    return (sort {$a->{lend}<=>$b->{lend}} @{$self->{alignment_segs}});
}



=over 4

=item add_alignment_segment()

B<Description:> Adds a CDNA::Alignment_segment object to the list of segments of this CDNA_alignment object.

B<Parameters:> CDNA::Alignment_segment object.

B<Returns:> none.

=back

=cut


sub add_alignment_segment() {
    my $self = shift;
    my $segment = shift;
    push (@{$self->{alignment_segs}}, $segment);
}


=over 4

=item delete_all_segments()

B<Description:> Empties the current list of Alignment_segment objects.

B<Parameters:> None.

B<Returns:> None.

=back

=cut




sub delete_all_segments () {
    my $self = shift;
    @{$self->{alignment_segs}} = (); #empty the array
}


sub set_error_flag () {
    my $self = shift;
    my $error = shift;
    if ($self->{error_flag}) {
        $self->{error_flag} .= $error;
    } else {
        $self->{error_flag} = $error;
    }
}


=over 4

=item get_error_flag()

B<Description:> Provides the status of the aligments validation.

B<Parameters:> none.

B<Returns:> [error_text|0]

A text string containing the error is returned if an error exists.  Otherwise, zero is returned indicating the lack of errors.

=back

=cut


sub get_error_flag () {
    my $self = shift;
    return ($self->{error_flag});
}


=over 4

=item toString()

B<Description:> Returns an alignment as lines of text providing coordinate information. 

B<Parameters:> none.

B<Returns:> alignment_string

=back

=cut


sub toString() {
    my $self = shift;
    my $text = "\n\nAlignment: orientation: " . $self->{orientation} . "\n"
        . "coords: " . $self->{lend} . "-" . $self->{rend} . "\n";
    
    my @segments = $self->get_alignment_segments();
    foreach my $segment (@segments) {
        $text .= $segment->toString();
    }
    return ($text);
}


####
sub to_GFF3_format {
    my $self = shift;
    my %preferences = @_;
    

    my $seq_id = $preferences{seq_id} || $self->{genome_acc} || confess "Need seq_id in preferences, or set genome_acc attribute of obj";
    my $match_id = $preferences{match_id} or confess "Error, require match_id attribute";
    my $source = $preferences{source} || "PASA";
    
    my $orientation = $self->get_orientation();
    


    my @alignment_segments = $self->get_alignment_segments();
    my $acc = $self->get_acc() or confess "Error, accession for cdna_alignment is not available";
    

    my $GFF3_text = "";

    foreach my $alignment_segment (@alignment_segments) {
        my ($genome_lend, $genome_rend) = $alignment_segment->get_coords();
        my ($cdna_lend, $cdna_rend) = sort {$a<=>$b} $alignment_segment->get_mcoords();
                

        my $gff_struct = { seq_id => $seq_id,
                           source => $source,
                           type => "cDNA_match",
                           lend => $genome_lend,
                           rend => $genome_rend,
                           strand => $orientation,
                           attributes => "ID=$match_id; Target=$acc $cdna_lend $cdna_rend +",
                       };

        if (my $per_id = $alignment_segment->get_per_id()) {
            $gff_struct->{score} = $per_id;
        }
        
        $GFF3_text .= &GFF_maker::get_GFF_line($gff_struct);
        
    }
    
    return ($GFF3_text);
    
}




####
sub to_GTF_format {
    my $self = shift;
    my %preferences = @_;
    

    my $seq_id = $preferences{seq_id} || $self->{genome_acc} || confess "Need seq_id in preferences, or set genome_acc attribute of obj";
    
	my $gene_id = $preferences{gene_id} || confess "Need gene_id in preferences";
	my $transcript_id = $preferences{transcript_id} || $self->get_acc();
    	

	my $source = $preferences{source} || "PASA";
    
    my $orientation = $self->get_orientation();
    

	my ($trans_lend, $trans_rend) = sort {$a<=>$b} $self->get_coords();

    my @alignment_segments = $self->get_alignment_segments();
    my $acc = $self->get_acc() or confess "Error, accession for cdna_alignment is not available";
    

    my $gtf_text = join("\t", ( $seq_id,
								$source,
								"transcript",
								$trans_lend,
								$trans_rend,
								".",
								$orientation,
								".",
								"gene_id \"$gene_id\"; transcript_id \"$transcript_id\";")
						) . "\n";
	
	
    foreach my $alignment_segment (@alignment_segments) {
        my ($genome_lend, $genome_rend) = $alignment_segment->get_coords();
        my ($cdna_lend, $cdna_rend) = sort {$a<=>$b} $alignment_segment->get_mcoords();
                

		my $per_id = $alignment_segment->get_per_id() || ".";
        
		$gtf_text .= join("\t", ( $seq_id,
								  $source,
								  "exon",
								  $genome_lend,
								  $genome_rend,
								  $per_id,
								  $orientation,
								  ".",
								  "gene_id \"$gene_id\"; transcript_id \"$transcript_id\";")
						  ) . "\n";
		
		
        
    }
    
    return ($gtf_text);
    
}


sub provide_cdna_segment_coords () {
    my $self = shift;
    my @segments = $self->get_alignment_segments();
    my $coord_mapper = $self->{genome_cdna_coord_mapper};
    foreach my $segment (@segments) {
        my ($lend, $rend) = $segment->get_coords();
        my $mlend = $coord_mapper->{$lend};
        my $mrend = $coord_mapper->{$rend};
        $segment->set_mcoords($mlend, $mrend);
    }
}

=over 4

=item remap_cdna_segment_coords()

B<Description:> Method is used on an assembled alignment to renumber the coordinates of the assembled cDNA product, starting at 1 and ending at the assembly length.  Orientation is determined by the validated spliced orientation.

B<Parameters:> none.

B<Returns:> none.

=back

=cut

sub remap_cdna_segment_coords () { ## Used when cDNA mapped coordinates aren't known, or when an assembly of other alignments was generated.
    ## reassigns cdna coords to alignment coords starting at 1 and ending at determined length.
    ## spliced orientation determines how the coords will map
    ## if spliced orient is ambiguous, aligned orient is used.

    my $self = shift;
    my %coord_mapper;
    
    
    my @alignment_segments = $self->get_alignment_segments(); #remember, already in increasing order.
    
    my $spliced_orient = $self->get_spliced_orientation();
    my $aligned_orient = $self->get_aligned_orientation();

    my $transcript_orientation = ($spliced_orient =~ /[\+\-]/) ? $spliced_orient : $aligned_orient;
    
    if ($transcript_orientation eq '-') {
        @alignment_segments = reverse @alignment_segments;
    }
    
    my $curr_pos = 0;
    foreach my $segment (@alignment_segments) {
        my ($lend, $rend) = $segment->get_coords();
        my $seglength = abs ($rend - $lend) + 1;
        my $mlend = $curr_pos + 1;
        my $mrend = $curr_pos + $seglength;
        $curr_pos += $seglength;
        if ($transcript_orientation eq '-') {
            ($mlend, $mrend) = ($mrend, $mlend);
        }
        
        $segment->set_orientation($transcript_orientation);
        #print "setting $mlend, $mrend\n";
        $coord_mapper{$lend} = $mlend;
        $coord_mapper{$rend} = $mrend;
        $segment->set_mcoords($mlend, $mrend);
    }

    ## reset aligned_orient to transcript_orient if different
    ## this should only happen when aligned orient and spliced orient are opposite.
    ## in which case, the aligned orient is reset to the spliced orient.
    
    if ($aligned_orient ne $transcript_orientation) {
        $self->set_orientation($transcript_orientation);
    }

}

=over 4
    
=item toToken()

B<Description:> similar to the toString() method, but returns a pretty line of text summarizing the alignment and splice site data. 

B<Parameters:> none.

B<Returns:> text_line

Here is an example:

orient(-/-) align: 103762(2613)-104359(2016)E<gt>CT....ACE<lt>104452(2015)-105229(1238)E<gt>CT....ACE<lt>105315(1237)-105482(1070)E<gt>CT....ACE<lt>105582(1069)-105842(809)E<gt>CT....ACE<lt>105935(808)-105985(758)E<gt>CT....ACE<lt>106071(757)-106316(512)E<gt>CT....ACE<lt>106394(511)-106619(286)E<gt>CT....ACE<lt>106712(285)-106879(118)E<gt>CT....ACE<lt>107427(117)-107543(1)

or

orient(+/+) align: 83074(1)-83318(245)E<gt>GT....AGE<lt>83637(246)-83702(311)E<gt>GT....AGE<lt>83796(312)-83846(362)E<gt>GT....AGE<lt>83938(363)-84017(442)E<gt>GT....AGE<lt>84308(443)-84352(487)E<gt>GT....AGE<lt>84467(488)-84507(528)


The orient specification includes (cDNA sequence alignment orientation/ spliced orientation).   These are different when the reverse-complement of the sequence is provided, ascertained by the aligned orientation with consensus splice sites.

=back

=cut

sub toToken () {
    my $self = shift;
    my $orientation = $self->get_orientation();
    my $spliced_orientation = $self->get_spliced_orientation();
    my @alignment_segments = $self->get_alignment_segments();
    my $assembled_token = "orient(a$orientation/s$spliced_orientation) align: ";
    for (my $i = 0; $i <= $#alignment_segments; $i++) {
        my $segment = $alignment_segments[$i];
        $assembled_token .= $segment->toToken();
        unless ($i == $#alignment_segments) {
            $assembled_token .= "....";
        }
    }
    return ($assembled_token);
}


sub set_num_segments {
    my $self = shift;
    my $num_segments = shift;
    $self->{num_segments} = $num_segments;
}

=over 4

=item get_num_segments()

B<Description:> method provides the number of segments composing an alignment.

B<Parameters:> none.

B<Returns:> int

=back

=cut

sub get_num_segments {
    my $self = shift;
    return ($self->{num_segments});
}

sub set_spliced_orientation {
    my $self = shift;
    my $orientation = shift;
    $self->{spliced_orientation} = $orientation;
}


=over 4

=item get_spliced_orientation ()

B<Description:> provides the validating spliced orientation for an alignment.  In some cases this will be different from the alignment orientation; for example, in cases where the cDNA sequence is provided in the reverse orientation.

B<Parameters:> none

B<Returns:> [+|-|undef()]

undef is returned if the alignment did not validate properly.  See get_error_flag()

=back

=cut


sub get_spliced_orientation {
    my $self = shift;
    return ($self->{spliced_orientation});
}


=over 4

=item set_acc()

B<Description:> Method sets the accession field of the cDNA sequence.   

B<Parameters:> string

Provide the accession for the cDNA corresponding to this alignment.

B<Returns:> none.

=back

=cut

sub set_acc () {
    my $self = shift;
    my $acc = shift;
    $self->{acc} = $acc;
}

=over 4

=item get_acc()

B<Description:> Method provides the accession for the cDNA in the alignment. 

B<Parameters:> none

B<Returns:> string

=back

=cut


sub get_acc () {
    my $self = shift;
    return ($self->{acc});
}


=over 4

=item set_fli_status()

B<Description:> sets the is_fli attribute of the cDNA alignment, indicative of a full-length insert clone (or complete cDNA sequence).

B<Parameters:> [1|0]

1 = true, 0 = false.

B<Returns:> [1|0]

=back

=cut


sub set_fli_status {
    my $self = shift;
    my $is_fli_status = shift;
    $self->{is_fli} = $is_fli_status;
}


=over 4

=item is_fli()

B<Description:>Provides the full-length insert status of the cDNA. 

B<Parameters:> none.

B<Returns:> [0|1]

=back

=cut

sub is_fli {
    my $self = shift;
    return ($self->{is_fli});
}






=over 4

=item toAlignIllustration()

B<Description:> Provides a single line of text which illustrates the gapped alignment.  See the example below.

B<Parameters:> none.

B<Returns:> string

Here is an example of an illustrated alignment:

------>     <--->   <-----  (-)asmbl_6711

=back

=cut


sub toAlignIllustration () {
    my ($self, $subtract, $rel_max, $max_line_chars) = @_;
    my $spliced_orient = $self->get_spliced_orientation();
    my $orient = $self->get_orientation();
    my @segments = $self->get_alignment_segments();
    my @chars = ();
    my $converter = sub {my $coord = shift; 
                         return ( int ( ($coord - $subtract)/$rel_max * $max_line_chars + 0.5));
                     };
    foreach my $segment (@segments) {
        my ($lend, $rend) = $segment->get_coords();
        my $l_rel = &$converter($lend);
        #print "lend: $lend -> l_rel: $l_rel\n" if $::SEE;
        my $r_rel = &$converter($rend);
        #print "rend: $rend -> r_rel: $r_rel\n" if $::SEE;
        
        for (my $i = $l_rel; $i <= $r_rel; $i++) {
            $chars[$i] = '-';
        }
        if ($segment->has_left_splice_junction()) {
            $chars[$l_rel] = '<';
        } elsif ( (! $segment->is_first()) && (! $segment->is_single_segment())) {
            $chars[$l_rel] = '|';
        }
        
        if ($segment->has_right_splice_junction()) {
            $chars[$r_rel] = '>';
        } elsif ( (! $segment->is_last()) && (! $segment->is_single_segment())) {
            $chars[$r_rel] = '|';
        }
    }
    
    #fill rest of line with spaces.
    for (my $i = 0; $i <= $#chars; $i++) {
        unless ($chars[$i]) {
            $chars[$i] = ' ';
        }
    }
    my $outline = join ("", @chars);
    my $acc = $self->get_acc();
    $acc =~ tr/\t\n\000-\037\177-\377/\t\n/d; #remove any control characters from accession.
    my $fli_status = ($self->is_fli()) ? " FL" : "";;
    return ($outline . "\t(a$orient/s$spliced_orient)" . $acc . $fli_status);
}


=over 4

=item get_gene_obj_via_alignment()

B<Description:> Creates a Gene_obj object based on an alignment using the Exons_to_geneobj.pm module.

B<Parameters:> 

B<Returns:> Gene_obj

The object returned is of the type Gene_obj defined in Gene_obj.pm

=back

=cut


sub get_gene_obj_via_alignment {
    my $self = shift;
    my $partial_info_href = shift; # { 5prime => 0|1, 3prime => 0|1 }, optional
    
    my $orient = $self->get_spliced_orientation();
    
	if ($orient eq '+' || '-') {
		return($self->_get_gene_obj_via_alignment_by_orient($orient, $partial_info_href));
	}
	elsif ($orient eq '?') {
		## find the orientation that provides the longest ORF
		
		my $plus_orient_gene = $self->_get_gene_obj_via_alignment_by_orient('+', $partial_info_href);
		my $minus_orient_gene = $self->_get_gene_obj_via_alignment_by_orient('-', $partial_info_href);
		
		if ($plus_orient_gene->get_CDS_length() >= $minus_orient_gene->get_CDS_length()) {
			return($plus_orient_gene);
		}
		else {
			return($minus_orient_gene);
		}
		
	}
	else {
		confess "cannot process spliced orientation of $orient ";
	}
}
	   





sub _get_gene_obj_via_alignment_by_orient {
	my ($self, $orient, $partial_info_href) = @_;
	
	my @alignment_segments = $self->get_alignment_segments();
    my %coords;
    foreach my $segment (@alignment_segments) {
        my ($end5, $end3) = $segment->get_coords();
        
        if ($orient eq '-') {
            ($end5, $end3) = ($end3, $end5); #force coordinates to contain orientation info.
        }
        $coords{$end5} = $end3;
    }
    my $genomic_seq_ref = $self->{genomic_seq};
    my $gene_obj;
    if ($genomic_seq_ref) {
	## find ORF
        $gene_obj = Exons_to_geneobj::create_gene_obj(\%coords, $genomic_seq_ref, $partial_info_href);
    } else {
        ## No ORF
        $gene_obj = new Gene_obj;
        $gene_obj->populate_gene_obj({}, \%coords);
    }
    
    return ($gene_obj);
}





=over 4

=item force_spliced_validation()

B<Description:> Routine used for testing purposes.  Any fake alignment can be created and set to validate to the corresponding orientation using this routine.

B<Parameters:> [+|-]

B<Returns:> none.

=back

=cut


sub force_spliced_validation {
    my $self = shift;
    my $orientation = shift;

    unless ($orientation eq '+' || $orientation eq '-') {
        croak ("cannot force spliced validation to $orientation.\n");
        return;
    }
    
    my ($right_splice, $left_splice) = ($orientation eq '+') ? ('XX','YY') : ('YY','XX');
    
    my @alignment_segments = $self->get_alignment_segments();
    foreach my $alignment_segment (@alignment_segments) {
        if ($alignment_segment->is_internal() || $alignment_segment->is_last()) {
            $alignment_segment->set_left_splice_junction(1);
            $alignment_segment->set_left_splice_site_chars($left_splice);
        }
        if ($alignment_segment->is_internal() || $alignment_segment->is_first()) {
            $alignment_segment->set_right_splice_junction(1);
            $alignment_segment->set_right_splice_site_chars($right_splice);
        }
        $alignment_segment->set_orientation($orientation);
    }
    
    $self->set_spliced_orientation($orientation);
}



=over 4

=item clone()

B<Description:>Clones a CDNA_alignment object into a new CDNA_alignment object with same attributes.  Performs a deep copy, so all alignment segments contained within the cloned CDNA_alignment object are also clones. 

B<Parameters:> none.

B<Returns:> new CDNA_alignment  

=back

=cut

sub clone {
    my $self = shift;
    my $packagename = ref $self;
    
    my $clone = {};
    bless ($clone, $packagename);
    foreach my $key (keys %$self) {
        $clone->{$key} = $self->{$key};
    }
    $clone->{alignment_segs} = [];
    foreach my $alignment_segment ($self->get_alignment_segments()) {
        $clone->add_alignment_segment($alignment_segment->clone());
    }
    
    return ($clone);
}




=over 4

=item get_genomic_seq_ref()

B<Description:> Returns a scalar refernence to the genomic sequence.

B<Parameters:> none.

B<Returns:> string_ref

=back

=cut



sub get_genomic_seq_ref {
    my $self = shift;
    return ($self->{genomic_seq});
}


=over 4

=item extractSplicedSequence()

B<Description:> Returns a string corresponding to the spliced cDNA sequence.

B<Parameters:> none.

B<Returns:> scalar

=back

=cut



sub extractSplicedSequence {
    my $self = shift;
    my ($genomic_seq_ref) = @_;
    my $genomic_seq = $genomic_seq_ref;

    unless ($genomic_seq) {
        $genomic_seq = $self->{genomic_seq};
        unless (ref $genomic_seq) {
            confess "Can't extract the spliced sequence when no genomic sequence reference is available.\n";
        }
    }
    
    my @segments = $self->get_alignment_segments();
    my $splicedSequence = "";
    my $toggle = 0;
    foreach my $segment (@segments) {
        my ($lend, $rend) = $segment->get_coords();
        my $length = abs ($rend - $lend) + 1;
        my $exonseq = substr ($$genomic_seq, $lend - 1, $length);

        # alternate case among segments to facilitate manual identification of junctions.
        if ($toggle) {
            $exonseq = lc $exonseq;
            $toggle = 0;
        }
        else {
            $exonseq = uc $exonseq;
            $toggle = 1;
        }
        
        $splicedSequence .= $exonseq;
    }
    my $orient = $self->get_orientation();
    if ($orient eq "-") {
        #reverse complement the sequence:
        $splicedSequence = reverse ($splicedSequence);
        $splicedSequence =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
    }
    return ($splicedSequence);
}





=over 4

=item get_cDNA_to_genomic_coordinates()

B<Description:> converts a cDNA sequence -relative coordinates to the corresponding genomic sequence coordinates.

B<Parameters:> @coordinates

list of integers

B<Returns:> @converted_coordinates

list of integers.   

=back

=cut

sub get_cDNA_to_genomic_coordinates {
    my $self = shift;
    my @coordinates = @_;
    
    my @ret_coordinates = ();
    
    my @segments = $self->get_alignment_segments();
    foreach my $cdna_coord (@coordinates) {
        my $corresponding_segment;
        foreach my $seg (@segments) {
            my ($mlend, $mrend) = sort {$a<=>$b} $seg->get_mcoords();
            if ($cdna_coord >= $mlend && $cdna_coord <= $mrend) {
                $corresponding_segment = $seg;
                last;
            }
        }
        unless (ref $corresponding_segment) {
            confess "Error, cDNA coordinate ($cdna_coord) not found in segment list: " . $self->toToken();
        }
        my ($lend, $rend) = $corresponding_segment->get_coords();
        my ($mlend, $mrend) = $corresponding_segment->get_mcoords();
        my $orient = $corresponding_segment->get_orientation();
        
        my $genomic_coord = undef;
        if ($orient eq '+') {
            my $diff = $cdna_coord - $mlend;
            $genomic_coord = $lend + $diff;
        } elsif ($orient eq '-') {
            ## mlend > mrend
            my $diff = $cdna_coord - $mrend;
            $genomic_coord = $rend - $diff;
        } 
        
        push (@ret_coordinates, $genomic_coord);
    }
    
    return (@ret_coordinates);
}




=item get_genomic_to_cDNA_coordinates()

B<Description:> converts coordinates in the genome to coordinates in the cDNA:

B<Parameters:> @coordinates

list of integers

B<Returns:> @converted_coordinates

list of integers.   Undef is returned for each position that could not be converted to the genomic coordinate system because it was not found within the aligned region.

=back

=cut


####
sub get_genomic_to_cDNA_coordinates {

    my $self = shift;
    my @coordinates = @_;
    
    my @ret_coordinates = ();
    
    my @segments = $self->get_alignment_segments();
    foreach my $genomic_coord (@coordinates) {
        my $corresponding_segment;
        foreach my $seg (@segments) {
            my ($lend, $rend) = $seg->get_coords();
            if ($genomic_coord >= $lend && $genomic_coord <= $rend) {
                $corresponding_segment = $seg;
                last;
            }
        }
        unless (ref $corresponding_segment) {
            confess "Error, genomic coordinate ($genomic_coord) not found in segment list:" . $self->toToken();
        }
        my ($lend, $rend) = $corresponding_segment->get_coords();
        my ($mlend, $mrend) = $corresponding_segment->get_mcoords();
        my $orient = $corresponding_segment->get_orientation();
        my $diff = $rend - $genomic_coord;
        my $cdna_coord = undef;
        if ($orient eq '+') {
            $cdna_coord = $mrend - $diff;
        } elsif ($orient eq '-') {
            $cdna_coord = $mrend + $diff;
        } 
        
        push (@ret_coordinates, $cdna_coord);
    }
    
    return (@ret_coordinates);
}



####
sub overlaps_genome_span {
	my ($self, $other_alignment) = @_;

	my ($lend_A, $rend_A) = sort {$a<=>$b} $self->get_coords();

	my ($lend_B, $rend_B) = sort {$a<=>$b} $other_alignment->get_coords();

	if (&_overlap($lend_A, $rend_A, $lend_B, $rend_B)) {
		return(1);
	}
	else {
		return(0);
	}
}


sub has_overlapping_segment {
	my ($self, $other_alignment) = @_;

	unless ($self->overlaps_genome_span($other_alignment)) {
		return(0);
	}
	

	my @self_segments = $self->get_alignment_segments();

	my @other_segments = $self->get_alignment_segments();

	foreach my $segment_A (@self_segments) {
		
		my ($lend_A, $rend_A) = sort {$a<=>$b} $segment_A->get_coords();


		foreach my $segment_B (@other_segments) {

			my ($lend_B, $rend_B) = sort {$a<=>$b} $segment_B->get_coords();
						
			if (&_overlap($lend_A, $rend_A, $lend_B, $rend_B) ) {
				
				return(1);
			}
			
		}
		
	}
	
	return(0); # no overlapping segment
}


=over 4

=item is_compatible()

B<Description:> Returns true (1) if this alignment is found compatible with the other_alignment_obj

B<Parameters:> ($other_alignment_obj, $fuzz_dist)

Alignments A and B are of type CDNA::CDNA_alignment

B<Returns:> [1|0]

Compatibility between alignment objects requires:
-within their region of overlap, introns are identical
-if both have spliced orientations, they must be on the same strand


=back

=cut


####
sub is_compatible {
    my $self = shift;
    my ($other_alignment, $fuzzlength) = @_;
    if (!defined $fuzzlength) {
        $fuzzlength = 0; ## no fuzzy termini allowed
    }
    
    ## The compatibility test requires:
    #  -alignments must have the same spliced orientation if not ?
    #  -alignments must overlap
    #  -alignments must have identical introns in their region of overlap (taking into account the fuzz distance for terminal exons)
    
    my ($a_lend, $a_rend) = $self->get_coords();
    my $a_spliced_orient = $self->get_spliced_orientation();
    my $a_num_segments = $self->get_num_segments();

    my ($b_lend, $b_rend) = $other_alignment->get_coords();
    my $b_spliced_orient = $other_alignment->get_spliced_orientation();
    my $b_num_segments = $other_alignment->get_num_segments();

    ## overlap test:
    unless (&_overlap($a_lend, $a_rend, $b_lend, $b_rend)) {
        return (0); # not compatible
    }

    ## transcribed orientation test:
    if ($a_spliced_orient ne $b_spliced_orient && $a_spliced_orient ne '?' && $b_spliced_orient ne '?') {
        # neither alignment is ambiguously oriented and they're transcribed on opposite strands:
        return (0); # not compatible
    }
    
    ## same introns test:
    if ($a_num_segments > 1 || $b_num_segments > 1) {
        my @a_introns = $self->get_intron_coords();
        my @b_introns = $self->get_intron_coords();
        
        my ($overlapping_lend, $overlapping_rend) = &_get_coords_of_overlap($a_lend, $a_rend, $b_lend, $b_rend);
        print "Overlapping coords between alignments: $overlapping_lend to $overlapping_rend\n" if $SEE;
        ## make adjustments to required overlap coordinates considering fuzzlength:
        if ($fuzzlength) {
            ## adjust left overlap boundary requirement
            $overlapping_lend = &_adjust_left_overlap_boundary_via_fuzzlength($overlapping_lend, $self, $other_alignment, $fuzzlength);
            $overlapping_rend = &_adjust_right_overlap_boundary_via_fuzzlength($overlapping_rend, $self, $other_alignment, $fuzzlength);
            print "\tadjusted overlapping coords to: $overlapping_lend to $overlapping_rend\n" if $SEE;
        }
        
        ## find introns within adjusted overlap range and ensure identity
        my @a_intron_coords = $self->get_intron_coords();
        my @b_intron_coords = $other_alignment->get_intron_coords();
        
        my @a_introns_in_range = &_get_overlapping_capped_introns($overlapping_lend, $overlapping_rend, \@a_intron_coords);
        my @b_introns_in_range = &_get_overlapping_capped_introns($overlapping_lend, $overlapping_rend, \@b_intron_coords);
        
        if (@a_introns_in_range || @b_introns_in_range) {
            ## ensure identity:
            my %all_introns;
            my %a_introns;
            foreach my $coordset (@a_introns_in_range) {
                my $key = join (",", @$coordset);
                $a_introns{$key} = 1;
                $all_introns{$key} = 1;
            }
            my %b_introns;
            foreach my $coordset (@b_introns_in_range) {
                my $key = join (",", @$coordset);
                $b_introns{$key} = 1;
                $all_introns{$key} = 1;
            }
            foreach my $intron_key (keys %all_introns) {
                unless ($a_introns{$intron_key} && $b_introns{$intron_key}) {
                    return (0); # not compatible, an intron difference in the overlapping region exists.
                }
            }
        }
        
    }

    ## if got this far, passed all compatibility tests.

    return (1); # yes, compatible.
}

####
sub _get_contained_coords {
    my ($lend, $rend, $coordsets_aref) = @_;
    my @contained_coords;
    foreach my $coordset (@$coordsets_aref) {
        my ($coord_lend, $coord_rend) = sort {$a<=>$b} @$coordset;
        if ($lend <= $coord_lend && $coord_rend <= $rend) {
            ## coordset contained
            push (@contained_coords, $coordset);
        }
    }
    return (@contained_coords);
}

####
# get introns that overlap lend and rend, and set termini of intron coords to these values if they extend beyond them.
sub _get_overlapping_capped_introns {
    my ($lend, $rend, $coordsets_aref) = @_;
    my @contained_coords;
    foreach my $coordset (@$coordsets_aref) {
        my ($coord_lend, $coord_rend) = sort {$a<=>$b} @$coordset;
        if ($lend <= $coord_rend && $rend >= $coord_lend) {
            ## coordset contained
            if ($coord_lend < $lend) {
                $coord_lend = $lend;
            }
            if ($coord_rend > $rend) {
                $coord_rend = $rend;
            }
            push (@contained_coords, [$coord_lend, $coord_rend]);
        }
    }
    return (@contained_coords);
}


####
sub _adjust_left_overlap_boundary_via_fuzzlength {
    my ($overlapping_lend, $alignment_a, $alignment_b, $fuzzlength) = @_;
    
    ## make adjustments to take into account fuzzlength and existing intron coordinates and adjacent segments
    ## we trust short aligment segments that precede an intron, even if they're shorter than the fuzzlength
    
    my $a_overlapping_segment = $alignment_a->find_segment_containing_coord($overlapping_lend);
    my $b_overlapping_segment = $alignment_b->find_segment_containing_coord($overlapping_lend);
    my @bounds = ();
    if ($a_overlapping_segment) {
        my ($a_seg_lend, $a_seg_rend) = $a_overlapping_segment->get_coords();
        push (@bounds, $a_seg_rend);
    }
    if ($b_overlapping_segment) {
        my ($b_seg_lend, $b_seg_rend) = $b_overlapping_segment->get_coords();
        push (@bounds, $b_seg_rend);
    }
    if (@bounds) {
        my $max_bound = max_coord(@bounds);
        my $delta = $max_bound - $overlapping_lend;
        if ($delta < 0) {
            confess "Error, delta left bound is less than zero";
        }
        my $fuzz_employed = min_coord($fuzzlength, $delta);
        $overlapping_lend += $fuzz_employed;
    }
    return ($overlapping_lend);
}

####
sub _adjust_right_overlap_boundary_via_fuzzlength {
    my ($overlapping_rend, $alignment_a, $alignment_b, $fuzzlength) = @_;
    
    ## make adjustments to take into account fuzzlength and existing intron coordinates and adjacent segments
    ## we trust short aligment segments that precede an intron, even if they're shorter than the fuzzlength
    
    my $a_overlapping_segment = $alignment_a->find_segment_containing_coord($overlapping_rend);
    my $b_overlapping_segment = $alignment_b->find_segment_containing_coord($overlapping_rend);
    my @bounds = ();
    if ($a_overlapping_segment) {
        my ($a_seg_lend, $a_seg_rend) = $a_overlapping_segment->get_coords();
        push (@bounds, $a_seg_lend);
    }
    if ($b_overlapping_segment) {
        my ($b_seg_lend, $b_seg_rend) = $b_overlapping_segment->get_coords();
        push (@bounds, $b_seg_lend);
    }
    if (@bounds) {
        my $min_bound = min_coord(@bounds);
        my $delta = $overlapping_rend - $min_bound;
        if ($delta < 0) {
            confess "Error, delta left bound is less than zero";
        }
        my $fuzz_employed = min_coord($fuzzlength, $delta);
        $overlapping_rend -= $fuzz_employed;
    }
    return ($overlapping_rend);
}





####
sub find_segment_containing_coord {
    my $self = shift;
    my ($coord) = @_;
    my @segments = $self->get_alignment_segments();
    foreach my $segment (@segments) {
        my ($lend, $rend) = $segment->get_coords();
        if ($lend <= $coord && $coord <= $rend) {
            return ($segment);
        }
    }
    return (undef); # none found
}


sub _overlap {
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

=over 4

=item encapsulates()

B<Description:> Returns true (1) if this alignment encapsulates the span of the other alignment object

B<Parameters:> ($other_alignment_obj, $fuzz_dist)

Alignments A and B are of type CDNA::CDNA_alignment

B<Returns:> [1|0]

=back

=cut


sub encapsulates {
    my $self = shift;
    my ($alignmentB, $fuzz_dist) = @_;
    my $alignmentA = $self;
    
    if (! defined $fuzz_dist) {
        $fuzz_dist = 0;
    }
    
    my ($alend, $arend) = $alignmentA->get_coords();
    my ($blend, $brend) = $alignmentB->get_coords();
    
    if ($blend + $fuzz_dist >= $alend && 
        $brend - $fuzz_dist <= $arend) 
    {
        return (1);
    } else {
        return (0);
    }
}

####
sub _get_coords_of_overlap {
    my ($a_lend, $a_rend, $b_lend, $b_rend) = @_;
    unless (&_overlap($a_lend, $a_rend, $b_lend, $b_rend)) {
        confess "Error, trying to get coordinates of overlapping region for two features that do not overlap!";
    }
    my $overlapping_lend = &max_coord($a_lend, $b_lend);
    my $overlapping_rend = &min_coord($a_rend, $b_rend);

    return ($overlapping_lend, $overlapping_rend);
}

####
sub min_coord {
    my @coords = @_;
    @coords = sort {$a<=>$b} @coords;
    my $min_coord = shift @coords;
    return ($min_coord);
}

####
sub max_coord {
    my @coords = @_;
    @coords = sort {$a<=>$b} @coords;
    my $max_coord = pop @coords;
    return ($max_coord);
}


1; #EOM





