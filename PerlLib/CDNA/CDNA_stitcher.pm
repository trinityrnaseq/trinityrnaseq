package main;
our $SEE;

package CDNA::CDNA_stitcher;
use strict;
use CDNA::CDNA_alignment;
use CDNA::Gene_obj_alignment_assembler;  #used for converting gene_obj to alignment.
use Carp;

my $FUZZDIST = 20; #allow FUZZDIST nt extension beyond splice site.

sub new {
    my $packagename = shift;
    my $self = {
	fuzzlength=>$FUZZDIST  #default setting.
	};
    bless ($self, $packagename);
    return ($self);
}


sub set_fuzzlength {
    my $self = shift;
    my $fuzzdist = shift;
    $self->{fuzzlength} = $fuzzdist;
}

sub stitch_alignments {
    my $self = shift;
    my ($template_alignment, $thread_alignment) = @_;

    my $template_acc = $template_alignment->get_acc();
    my $thread_acc = $thread_alignment->get_acc();

    print "Stitching Template:\n$template_acc: " . $template_alignment->toToken() 
        . "\nby\n$thread_acc: " . $thread_alignment->toToken() . "\n\n" if $main::SEE;

    ## don't even think about stitching alignments of opposite spliced orienatations!!!

    my $template_spliced_orient = $template_alignment->get_spliced_orientation();
    my $thread_spliced_orient = $thread_alignment->get_spliced_orientation();

    if ($template_spliced_orient ne '?' && $thread_spliced_orient ne '?') {
        ## both have assigned transcribed orients
        if ($template_spliced_orient ne $thread_spliced_orient) {
            confess "Error, trying to stitch together oppositely transcribed alignments.  That's impossible. ";
        }
    }
    
    
    my $spliced_orient = ($thread_spliced_orient ne '?') ? $thread_spliced_orient : $template_spliced_orient;
    
    my $fuzzlength = $self->{fuzzlength};

    ## Algorithm
    ## -anchor the first thread segment to the template alignment
    ## -determine the boundaries of the matching anchor and first threaded segment
    ## -add all preceding template segments to the stitched alignment if splice compatible.
    ## -add internal cdna alignment segments
    ## -add terminal template alignment segments if splice compatible.

    my @template_segments = $template_alignment->get_alignment_segments();
    my @threaded_segments = $thread_alignment->get_alignment_segments();
    
    my $single_threaded_segment = ($#threaded_segments == 0) ? 1:0; #indicate only a single segment exists.

    my @stitched_segments;
    ## Anchor first threaded segment to the template alignment:
    my ($thread_lend, $thread_rend) = $threaded_segments[0]->get_coords();
    my $anchor_pos = -1;
    for (my $i = 0; $i <= $#template_segments; $i++) {
        my ($lend, $rend) = $template_segments[$i]->get_coords();
        if ($lend < $thread_rend && $rend > $thread_lend) { #overlap
            $anchor_pos = $i;
            last;
        }
    }
    print "Anchoring first cDNA segment. Anchor point: $anchor_pos\n" if $SEE;
    if ($anchor_pos > -1) { #found an anchor point in template alignment.
        
        
        ## add stitched segment
        my $anchor_segment = $template_segments[$anchor_pos];
        my ($anchor_lend, $anchor_rend) = $anchor_segment->get_coords();
        my $first_threaded_segment = $threaded_segments[0];
        my ($thread_lend, $thread_rend) = $first_threaded_segment->get_coords(); 
        # initialize to largest spread
        my ($new_lend) = ($anchor_lend < $thread_lend) ? $anchor_lend : $thread_lend;
        my ($new_rend) = ($anchor_rend > $thread_rend) ? $anchor_rend : $thread_rend;
        # adjust based on splice boundaries
        my ($has_left_splice_junction, $has_right_splice_junction);
        $has_left_splice_junction = $anchor_segment->has_left_splice_junction(); #initialize.
        if ($anchor_segment->has_left_splice_junction() && ($anchor_lend - $thread_lend <= $fuzzlength)) {
            $new_lend = $anchor_lend;
        } elsif ($thread_lend < $anchor_lend) {
            #using threaded segment as initial stitched segment
            $has_left_splice_junction = $first_threaded_segment->has_left_splice_junction();
        }
        
        $has_right_splice_junction = $first_threaded_segment->has_right_splice_junction();
        if ($first_threaded_segment->has_right_splice_junction()) {
            $new_rend = $thread_rend;
        } elsif ($anchor_segment->has_right_splice_junction() && ($thread_rend - $anchor_rend <= $fuzzlength)) {
            $new_rend = $anchor_rend;
            $has_right_splice_junction = $anchor_segment->has_right_splice_junction();
        }
        
        ## Add preceding template segments if current segment is left-splice compatible.
        if ($has_left_splice_junction) {
            print "Adding all template segments before position $anchor_pos.\n" if $SEE;
            ## add all template segments to the stitched alignment preceding the anchor point:
            for (my $i = 0; $i < $anchor_pos; $i++) {
                push (@stitched_segments, $template_segments[$i]->clone());
            }
        }
        
        ## add the current segment.
        my $newsegment = $anchor_segment->clone();
        $newsegment->set_coords($new_lend, $new_rend);
        $newsegment->set_left_splice_junction($has_left_splice_junction);
        $newsegment->set_right_splice_junction($has_right_splice_junction);
        print "adding current segment: lsplice: $has_left_splice_junction, rsplice: $has_right_splice_junction\n" if $SEE;
        push (@stitched_segments, $newsegment);
        
    } else { #no anchor pos, simply add the first threaded segment:
        print "No anchor position, simply adding the first threaded segment.\n" if $SEE;
        push (@stitched_segments, $threaded_segments[0]->clone());
    }
    
    ## Add all but last threaded segments to the stitched alignment:
    
    for (my $i = 1; $i < $#threaded_segments; $i++) {
        print "Adding internal cDNA alignment segment, index: $i\n" if $SEE;
        push (@stitched_segments, $threaded_segments[$i]->clone());
    }
    
    ## Anchor the last threaded segment to a segment within the template model:
    print "Anchoring the terminal segment.\n" if $SEE;
    $anchor_pos = -1;
    my $last_threaded_segment = $threaded_segments[$#threaded_segments];
    ($thread_lend, $thread_rend) = $last_threaded_segment->get_coords();
    # find last template segment overlapping last cdna segment
    for (my $i = $#template_segments; $i >= 0; $i--) {
        my ($lend, $rend) = $template_segments[$i]->get_coords();
        if ($lend < $thread_rend && $rend > $thread_lend) { #overlap
            $anchor_pos = $i;
            last;
        }
    }
    print "Terminal anchor position: $anchor_pos\n" if $SEE;
    if ($anchor_pos > -1) {
        #build composite segment
        my $anchor_segment = $template_segments[$anchor_pos];
        my $has_left_splice_junction = $last_threaded_segment->has_left_splice_junction(); #initialize.
        my $has_right_splice_junction = $anchor_segment->has_right_splice_junction();
        my ($anchor_lend, $anchor_rend) = $anchor_segment->get_coords();
        my $new_lend = ($anchor_lend < $thread_lend) ? $anchor_lend : $thread_lend;
        my $new_rend = ($anchor_rend > $thread_rend) ? $anchor_rend : $thread_rend;
        if ($last_threaded_segment->has_left_splice_junction()) {
            $new_lend = $thread_lend;
        } elsif ($anchor_segment->has_left_splice_junction() && ($anchor_lend - $thread_lend <= $fuzzlength)) {
            $new_lend = $anchor_lend;
            $has_left_splice_junction = $anchor_segment->has_left_splice_junction();
        }
        if ($anchor_segment->has_right_splice_junction() && ($thread_rend - $anchor_rend <= $fuzzlength)) {
            $new_rend = $anchor_rend;
        } else {
            $has_right_splice_junction = $last_threaded_segment->has_right_splice_junction();
        }
        if ($single_threaded_segment) { #only one segment (first == last)
            #just update existing status:
            print "Only one threaded segment, updating it's rend coord and status.\n" if $SEE;
            my $last_stitched_segment = $stitched_segments[$#stitched_segments];
            my ($curr_lend, $curr_rend) = $last_stitched_segment->get_coords();
            $last_stitched_segment->set_coords($curr_lend, $new_rend);
            $last_stitched_segment->set_right_splice_junction($has_right_splice_junction);
        } else {
            my $newsegment = $anchor_segment->clone();
            $newsegment->set_coords($new_lend, $new_rend);
            $newsegment->set_left_splice_junction($has_left_splice_junction);
            $newsegment->set_right_splice_junction($has_right_splice_junction);
            print "adding composite terminal segment. Lsplice: $has_left_splice_junction, Rsplice: $has_right_splice_junction.\n" if $SEE;
            push (@stitched_segments, $newsegment);
        }
    } else {
        print "No matching terminal exon in template.\n" if $SEE;
        if ($single_threaded_segment) {
            print "Single threaded segment remains unchanged.\n" if $SEE;
        } else {
            
            print "Adding last cDNA exon.\n" if $SEE;
            push (@stitched_segments, $last_threaded_segment->clone());
        }
    }
    
    
    ## Add terminal template segments if the last stitched segment contains a splice junction
    
    my $last_stitched_segment =  $stitched_segments[$#stitched_segments];
    if ($last_stitched_segment->has_right_splice_junction()) {
        print "Adding terminal template segments:\n" if $SEE;
        my ($last_lend, $last_rend) = $last_stitched_segment->get_coords();
        $anchor_pos = -1;
        for (my $i = $#template_segments; $i >= 0; $i--) {
            my ($anchor_lend, $anchor_rend) = $template_segments[$i]->get_coords();
            if ($anchor_rend > $last_lend && $anchor_lend < $last_rend) { #overlap
                $anchor_pos = $i;
                last;
            }
        }
        print "Last overlapping alignment segment in index found at: $anchor_pos\n" if $SEE;
        if ($anchor_pos > -1) {
            for (my $i = $anchor_pos + 1; $i <= $#template_segments; $i++) {
                print "\tadding terminal template segment index: $i\n" if $SEE;
                push (@stitched_segments, $template_segments[$i]->clone());
            }
        }
    }
    
    
    ## Done stitching segments together.  Create new alignment:
    my $cdna_length = 0;
    foreach my $segment (@stitched_segments) {
        my $length = $segment->get_length();
        $cdna_length += $length;
    }
    
    my $seq_ref = $template_alignment->get_genomic_seq_ref();
    
    my $stitched_alignment = new CDNA::CDNA_alignment ($cdna_length, \@stitched_segments, $seq_ref); #aligned orientation auto assigned.
    
    $stitched_alignment->set_spliced_orientation($spliced_orient);
    $stitched_alignment->set_acc($template_alignment->get_acc());
    
    print "Stitched alignment: " . $stitched_alignment->toToken() . "\n" if $SEE;
    
    return ($stitched_alignment);
}



sub stitch_alignment_into_gene {
    my $self = shift;
    my ($gene_obj, $cdna_obj, $sequence_ref) = @_;
    
    my $gene_obj_assembler = new CDNA::Gene_obj_alignment_assembler($sequence_ref);
    my $gene_based_cdna_obj = $gene_obj_assembler->gene_obj_to_cdna_alignment($gene_obj);
    
    my $stitched_alignment = $self->stitch_alignments($gene_based_cdna_obj, $cdna_obj);
    $stitched_alignment->set_fli(1);
    
    my $partials_href = $self->_analyze_partial_status($gene_obj);
    my $stitched_gene_obj = $stitched_alignment->get_gene_obj_via_alignment($partials_href);
    
    return ($stitched_gene_obj);
    
}


#### 
sub _analyze_partial_status {
    my $self = shift;
    my ($gene_obj) = shift;
    my %partials;
    if ($gene_obj->is_5prime_partial()) {
        $partials{"5prime"} = 1;
    }
    if ($gene_obj->is_3prime_partial()) {
        $partials{"3prime"} = 1;
    }
    return (\%partials);
}




1; #EOM





	     



    
    
    
