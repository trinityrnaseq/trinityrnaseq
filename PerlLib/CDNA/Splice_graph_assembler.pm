#!/usr/local/bin/perl

package main;
our $SEE;

package CDNA::Splice_graph_assembler;

use strict;
use warnings;
use Carp;
use Overlap_piler;
use CDNA::CDNA_alignment;
use CDNA::Alignment_segment;
use Data::Dumper;
no warnings "recursion";

our $FUZZ_DIST = 20; # bp's not to trust at non-splice termini; untrustworthy alignment extensions.

## allowable connections between gene structure components:
my %ACCEPTABLE_CONNECTIONS = ( "terminal_left_exon" => {"intron" => 1},
                               "internal_exon" => {"intron" => 1},
                               
                               "intron" => {"internal_exon" => 1,
                                            "terminal_right_exon" => 1,
                                        },
                               );


####
sub new {
    my $packagename = shift;

    my $self = {
        _graph_nodes => [], # each and every graph node (main repository)
        _incoming_alignments => [], # alignments to assemble
        _assemblies => [], # resulting alignment assemblies
        _unincorporated_alignments => [],

        _valid_splice_paths => [], # Splice_graph_path objects 
        _assembled_splice_paths => [], # compatible splice paths chained into splice path assemblies.
        
        ## other helpers
        _graph_node_hashkey_lookup => {},  # key is lend,rend,type,orient
        _graph_node_via_nodeID => {}, # node access by nodeID
        
        ## various node types: (all subsets of _graph_nodes and accessible via helpers above.
        _internal_exons => [],  
        _introns => [],
        _terminal_left_exons => [],
        _terminal_right_exons => [],
        _singleton_exons => [],
        
        ## misc
        _terminal_exon_nodeID_to_nonsplice_position_list => {}, # track coords of non-splice site ends in terminal exon supports
    };
    
    bless ($self, $packagename);

    return ($self);
}

####
sub assemble_alignments {
    my $self = shift;
    my @alignments = @_;
    
    ## build the splicing graph.
    $self->build_splicing_graph(@alignments);
    
    ## Chain together splice paths that are compatible:
    print "-Chaining compatible splice paths\n" if $SEE;
    $self->_chain_compatible_splice_paths();
    
    ## Extend splice paths into maximally scoring complete paths with terminal exons
    print "-Extending splice paths into maximal structures.\n" if $SEE;
    $self->_extend_splice_paths_to_termini();  ## products are included in the assemblies list.

    ## Add the singleton assemblies
    $self->_append_singletons_to_assembly_list(); ## convert the singletons to cdna alignments and add to assembly list

    ## assign the incoming alignments to the assemblies that contain them.
    print "-Assigning transcript alignments to assemblies\n" if $SEE;
    my $have_unincorporated_alignments_flag = $self->_correlate_assemblies_with_incoming_alignments();
    if ($have_unincorporated_alignments_flag) {
        ## find paths from unincorporated terminal segments.
        print "-Not all alignments were included.  Exploring assemblies from uninorporated alignments\n" if $SEE;
        $self->_explore_assemblies_from_unincorporated_alignments();
        print "-Assigning transcript alignments to assemblies, again...\n" if $SEE;
        $have_unincorporated_alignments_flag = $self->_correlate_assemblies_with_incoming_alignments();
        if ($have_unincorporated_alignments_flag) {
            ## something horribly has gone wrong.  
            confess "Error, not all incoming alignments are accounted for!\n" . $self->toString();
        }
    }
    
    print $self->toString() if $SEE;
}


####
sub build_splicing_graph {
    my $self = shift;
    my @alignments = @_;
 
    $self->set_incoming_alignments(@alignments);

    ## going to process the alignments in several phases:
    #  -decompose alignments into structural components
    #  -identify valid paths to connect components
    #  -extract assemblies using valid paths
    #  -associate transcripts with assemblies
 
    ## first, examine the alignments with the most segments first:
    @alignments = reverse sort {$a->get_num_segments() <=> $b->get_num_segments()} @alignments;
    
    
    ## build components of the splice graph
    ## add internal exons and all introns to splice graph: (unambiguous structures)
    ##    terminal structures are not as well defined.  Apply these later.

    print "-Decomposing alignments into gene structure components.\n" if $SEE;
    foreach my $alignment (@alignments) {
        if ($alignment->get_num_segments() > 1) {
            $self->_decompose_unambiguous_alignment_structures_add_nodes($alignment);
        }
    }
    
    ## Add terminal exons:
    ## Only add terminal exons where it's clear that it's not just part of an existing internal exon
    
    print "-Adding terminal exon segments\n" if $SEE;
    my @single_segments; # capture those alignments w/o introns
    
    foreach my $alignment (@alignments) {
        my $segment_orient = $alignment->get_spliced_orientation();
        if ($alignment->get_num_segments() > 1) {
            foreach my $segment ($alignment->get_alignment_segments()) {
                if ($segment->is_first() || $segment->is_last()) {
                    
                    my $exon_type = ($segment->is_first()) ? "terminal_left_exon" : "terminal_right_exon";
                    
                    $self->_try_terminal_exon_addition($exon_type, $segment, $segment_orient);
                    
                    ## at this point, introns and internal exons are scored only by perfect support (exact boundaries).
                    ## terminal exons, on the other hand, have evidence scored based on boundary match.
                    ## Now, need to augment scores of internal exons that encompass terminal exons
                    
                    
                    my ($segment_lend, $segment_rend) = $segment->get_coords();
                    
                    $self->_augment_internal_exon_scores_using_terminal_alignment_segment($exon_type, $segment_lend, 
                                                                                          $segment_rend, $segment_orient);
                }
            }
        }
        else {
            push (@single_segments, $alignment);
        }
    }
    
    ## reconstruct some internal exons based on overlapping right/left terminal exons
    print "-Merging overlapping right-to-left terminal exons\n" if $SEE;
    $self->_merge_overlapping_right_to_left_terminal_exons();
    
    ## apply intron-less segments as evidence to internal and terminal segments
    ## and instantiate single exons where they're not simply supporting existing structures.
    print "-Analyzing intronless alignments\n" if $SEE;
    $self->_analyze_intronless_alignments(@single_segments);
    
    print "-Building the splice graph\n" if $SEE;
    $self->_build_splice_graph(); # could/should have done this earlier during alignment parse for efficiency.
 
    return;
   
}

####
sub _find_internal_exons_encompassing_coords_and_share_boundary {
    my $self = shift;
    my ($exon_type, $segment_lend, $segment_rend, $segment_orient) = @_;
    
    unless ($exon_type =~ /left|right/) { 
        confess "invalide terminal exon type: $exon_type\n";
    }
    
    my @internal_exons_found;

     foreach my $internal_exon (@{$self->{_internal_exons}}) {
        my $internal_exon_orient = $internal_exon->get_orient();
        my ($internal_exon_lend, $internal_exon_rend) = $internal_exon->get_coords();
        
        ## if internal exon encompasses the terminal exon, add it to its evidence collection
        if ( ($segment_orient eq $internal_exon_orient)
             &&
             ($internal_exon_lend <= ($segment_lend + $FUZZ_DIST) && ($segment_rend - $FUZZ_DIST) <= $internal_exon_rend) # internal encompasses it
             &&
             ( ($exon_type eq "terminal_right_exon" && $internal_exon_lend == $segment_lend)
                 || 
               ($exon_type eq "terminal_left_exon" && $internal_exon_rend == $segment_rend) ) ## a splice boundary in common
           )
        {
            push (@internal_exons_found, $internal_exon);
        }
    }
    return (@internal_exons_found);
}

sub _find_terminal_exons_encompassing_segment {
    my $self = shift;
    my ($lend, $rend, $orient) = @_;

    my @nodes_encompassing_segment;
    foreach my $node_obj (@{$self->{_terminal_left_exons}}, @{$self->{_terminal_right_exons}}) {
        my ($node_lend, $node_rend) = $node_obj->get_coords();
        my $node_orient = $node_obj->get_orient();
        if ( ($node_orient eq $orient || $orient eq '?') && 
            ($lend + $FUZZ_DIST >= $node_lend) &&
            ($rend - $FUZZ_DIST <= $node_rend) ) {
            push (@nodes_encompassing_segment, $node_obj);
        }
    }

    return (@nodes_encompassing_segment);
}




####
sub get_assemblies {
    my $self = shift;
    return (@{$self->{_assemblies}});
}

####
sub get_incoming_alignments {
    my $self = shift;
    return (@{$self->{_incoming_alignments}});
}

####
sub get_unincorporated_alignments {
    my $self = shift;
    return (@{$self->{_unincorporated_alignments}});
}


####
sub _add_alignment_assembly {
    my $self = shift;
    my (@cdna_alignments) = @_;
    
    push (@{$self->{_assemblies}}, @cdna_alignments);
    return;
}
   


####
sub set_incoming_alignments {
    my $self = shift;
    my @alignments = @_;
    @{$self->{_incoming_alignments}} = @alignments;
    return;
}


####
sub _decompose_unambiguous_alignment_structures_add_nodes {
    my $self = shift;
    my ($alignment) = @_;
    
    my $spliced_orient = $alignment->get_spliced_orientation();
    
    my @path_nodes;

    ## Add internal exons
    my @segments = $alignment->get_alignment_segments();
    foreach my $segment (@segments) {
        #print "seg: " . $segment->toString() . "\n";
        if ($segment->is_internal()) {
            my ($exon_lend, $exon_rend) = $segment->get_coords();
            my $node = $self->_add_internal_exon($exon_lend, $exon_rend, $spliced_orient);
            push (@path_nodes, $node);
        }
    }
    
    ## Add introns
    my @intron_coords = $alignment->get_intron_coords();
    
    foreach my $intron_coordset (@intron_coords) {
        my ($intron_lend, $intron_rend) = @$intron_coordset; #already sorted
        my $node = $self->_add_intron($intron_lend, $intron_rend, $spliced_orient);
        push (@path_nodes, $node);
    }

    ## add a path
    $self->_extract_and_add_path_from_node_list(@path_nodes);


    return;
}


####
sub _extract_and_add_path_from_node_list {
    my $self = shift;
    my @path_nodes = @_;

    ## sort by genome order:
    @path_nodes = sort {$a->{lend} <=> $b->{lend}} @path_nodes;

    my @ordered_nodeID_list;
    my @coords;
    foreach my $path_node (@path_nodes) {
        my $nodeID = $path_node->get_nodeID();
        push (@ordered_nodeID_list, $nodeID);
        push (@coords, $path_node->get_coords());
        
    }
    my $orient = $path_nodes[0]->get_orient();
    @coords = sort {$a<=>$b} @coords;
    my $lend = shift @coords;
    my $rend = pop @coords;
    my $splice_graph_path = Splice_graph_path->new($lend, $rend, $orient, \@ordered_nodeID_list);
    
    $self->_add_splice_graph_path($splice_graph_path);
    
    return;
}

####
sub _try_terminal_exon_addition {
    my $self = shift;
    my ($exon_type, $segment, $orient) = @_;
    
    my ($segment_lend, $segment_rend) = $segment->get_coords();
        
    print "Examining terminal exon: $exon_type, $segment_lend, $segment_rend\n" if $SEE;
    
    ## only consider this a genuine terminal exon if it doesn't appear to be part of an 
    ## existing internal exon from a more complete alignment
    if (my @internal_segments = $self->_find_internal_exons_encompassing_coords_and_share_boundary($exon_type, $segment_lend, $segment_rend, $orient)) {
        if ($SEE) { 
            print "$exon_type, $segment_lend-$segment_rend, $orient, found already represented by internal exons:\n";
            foreach my $internal_segment (@internal_segments) {
                print "\t" . $internal_segment->toString() . "\n";
            }
        }
        return;
    }
    
    my $exon = $self->_find_existing_terminal_exon ($exon_type, $segment_lend, $segment_rend, $orient);
    if ($exon) {
        print "-found existing terminal exon with shared boundary: " . $exon->toString() . "\n" if $SEE;
        my $nodeID = $exon->get_nodeID();
        ## check for boundary adjustment:
        my ($exon_lend, $exon_rend) = $exon->get_coords();
        my $nonsplice_coord = undef;
        if ($exon_type eq "terminal_left_exon") {
            $nonsplice_coord = $segment_lend;
            if ($segment_lend < $exon_lend) {
                print "-extending left boundary to $segment_lend\n" if $SEE;
                $exon->set_coords($segment_lend, $exon_rend); # extend left boundary
            }
        }
        elsif ($exon_type eq "terminal_right_exon") {
            $nonsplice_coord = $segment_rend;
            if ($segment_rend > $exon_rend) {
                print "-extending right boundary to $segment_rend\n" if $SEE;
                $exon->set_coords($exon_lend, $segment_rend);
            }
        }
        else {
            confess "Error, exon type $exon_type not accounted for. "; # should never get here anyway
        }
        $exon->increment_evidence_support(); ## account for extra evidence
        $self->_add_to_terminal_exon_nonsplice_position_list($nodeID, $nonsplice_coord);
        
    }
    else {
        ## add new terminal exon:
        $self->_add_graph_node($exon_type, $segment_lend, $segment_rend, $orient);
    }
    
    return;
}

####
sub _add_to_terminal_exon_nonsplice_position_list {
    my $self = shift;
    my ($nodeID, $nonsplice_coord) = @_;

    my $nodeID_to_pos_list_href = $self->{_terminal_exon_nodeID_to_nonsplice_position_list};
    
    my $list_aref = $nodeID_to_pos_list_href->{$nodeID};
    unless (ref $list_aref) {
        $list_aref = $nodeID_to_pos_list_href->{$nodeID} = [];
    }
    push (@$list_aref, $nonsplice_coord);

    return;
}


####
sub _find_existing_terminal_exon {
    my $self = shift;
    my ($exon_type, $segment_lend, $segment_rend, $orient) = @_;
    my $exon_list_aref = undef;
    if ($exon_type eq "terminal_left_exon") {
        $exon_list_aref = $self->{_terminal_left_exons};
    }
    elsif ($exon_type eq "terminal_right_exon") {
        $exon_list_aref = $self->{_terminal_right_exons};
    }
    else {
        confess "exon_type $exon_type not accepted for terminal exons";
    }

    foreach my $exon (@$exon_list_aref) {
        my ($lend, $rend) = $exon->get_coords();
        if ( 
             ($exon->get_orient() eq $orient) && 
               (
                    ($exon_type eq "terminal_left_exon" && $rend == $segment_rend)
                             ||
                    ($exon_type eq "terminal_right_exon" && $lend == $segment_lend) 
               )
           ) 
        {
            
            return ($exon); # found it!
        }
    }

    return (undef); # didn't find one.
}

####
sub _augment_internal_exon_scores_using_terminal_alignment_segment {
    my $self = shift;
    my ($exon_type, $segment_lend, $segment_rend, $segment_orient) = @_;
    
    my @relevant_internal_exons = $self->_find_internal_exons_encompassing_coords_and_share_boundary($exon_type, $segment_lend, $segment_rend, $segment_orient);
    foreach my $internal_exon (@relevant_internal_exons) {
        
        $internal_exon->increment_evidence_support();
    }
    
    return;
}

####
sub _merge_overlapping_right_to_left_terminal_exons {
    my $self = shift;
    
    ## looking for this situation:
    #   <-----------     right terminal exon
    #         --------->   left terminal exon
    # that can be merged into:
    #   <-------------->   an internal exon
    #
    
    my %nodeIDs_targeted_for_removal; # if they fully overlap, construct a nice internal exon w/o extensions beyond splice boundary

    foreach my $right_terminal_exon (@{$self->{_terminal_right_exons}}) {
        my $right_orient = $right_terminal_exon->get_orient();
        my ($right_lend, $right_rend) = $right_terminal_exon->get_coords();
        my $right_nodeID = $right_terminal_exon->get_nodeID();

        foreach my $left_terminal_exon (@{$self->{_terminal_left_exons}}) {
            my $left_orient = $left_terminal_exon->get_orient();
            my ($left_lend, $left_rend) = $left_terminal_exon->get_coords();
            my $left_nodeID = $left_terminal_exon->get_nodeID();
            
            unless ($right_orient eq $left_orient) { next; } # must be transcribed on same strand!

            ## check for overlap:
            unless ($left_lend <= $right_rend && $left_rend >= $right_lend) { next;}

            ## make sure right's left splice  is before left's right splice (doh! should have better names).
            unless ($right_lend < $left_rend) { next; }
            
            ## check for extensions:
            my $left_overhang = $right_lend - $left_lend;
            my $right_overhang = $right_rend - $left_rend;
            
            my $merge_flag = 0; 
            if ($left_overhang <= $FUZZ_DIST && $right_overhang <= $FUZZ_DIST) {
                ## nice merge, as in illustration
                $merge_flag = 1;
                ## also, target these for deletion now.
                $nodeIDs_targeted_for_removal{$right_nodeID} = 1;
                $nodeIDs_targeted_for_removal{$left_nodeID} = 1;
            }
            else {
                ## must check position lists to see if transcripts are contained that have 
                ## nicely overlapping boundaries
                if ($self->_right_left_terminal_exons_overlap_via_position_lists($right_nodeID, $left_nodeID, $right_lend, $left_rend)) {
                    $merge_flag = 1; 
                    if ($left_overhang <= $FUZZ_DIST) { ## check for insufficient overhang
                        $nodeIDs_targeted_for_removal{$left_nodeID} = 1; 
                    }
                    if ($right_overhang <= $FUZZ_DIST) {
                        $nodeIDs_targeted_for_removal{$right_nodeID} = 1;
                    }
                }
            }
            if ($merge_flag) {
                $self->_merge_terminal_exons($right_terminal_exon, $left_terminal_exon);
            }
        }
    }
    
    ## process deletions:
    if (%nodeIDs_targeted_for_removal) {
        $self->_purge_nodes(%nodeIDs_targeted_for_removal);
    }
    
    return;
}

####
sub _merge_terminal_exons {
    my $self = shift;
    my ($right_terminal_exon, $left_terminal_exon) = @_;
    
    my ($right_lend, $right_rend) = $right_terminal_exon->get_coords();
    my ($left_lend, $left_rend) = $left_terminal_exon->get_coords();
    
    my ($exon_lend, $exon_rend) = ($right_lend, $left_rend); ## splice junctions for new internal exon
    my $orient = $right_terminal_exon->get_orient();
    if ($orient ne $left_terminal_exon->get_orient()) {
        confess "Error, trying to merge two terminal exons with opposite transcriptional orientations!";
    }
    
    if ($self->_graph_node_exists("internal_exon", $exon_lend, $exon_rend, $orient)) {
        confess "Error, trying to merge two terminal exons into an internal exon that already exists!";
    }
    
    my $node = $self->_add_graph_node("internal_exon", $exon_lend, $exon_rend, $orient);
    ## add to it the evidence from the other nodes being merged:
    my $evidence_support_to_add = $left_terminal_exon->get_num_evidence_support() + $right_terminal_exon->get_num_evidence_support();
    $evidence_support_to_add -= 1; # node already has value of one.
    if ($evidence_support_to_add) {
        $node->increment_evidence_support($evidence_support_to_add);
    }
    return;
}


####
sub _add_splice_graph_path {
    my $self = shift;
    my ($splice_graph_path) = @_;
    
    ## add it as long as:
    #    -it's not a subpath of an already existing path
    #    -if an existing path is a subpath of this, remove it and replace it with this path
    
    my @current_valid_splice_paths = $self->_get_valid_splice_paths();

    ## check to see if splice_graph_path is already represented in the current path list:
    foreach my $current_valid_splice_path (@current_valid_splice_paths) {
        if ($splice_graph_path->is_subpath_of($current_valid_splice_path)) {
            return; # nothing to do; path is already included as a subset of the current path set
        }
    }
    
    # if got this far, our new path is not already fully represented.
    #  add it, and any other existing splice paths that are not a subpath of it.
    my @new_valid_splice_paths = ($splice_graph_path);
    foreach my $current_valid_splice_path (@current_valid_splice_paths) {
        if (! $current_valid_splice_path->is_subpath_of($splice_graph_path)) {
            push (@new_valid_splice_paths, $current_valid_splice_path);
        }
    }

    $self->_set_valid_splice_paths(@new_valid_splice_paths);
    
    return;
}

####
sub _get_valid_splice_paths {
    my $self = shift;
    return (@{$self->{_valid_splice_paths}});
}

####
sub _set_valid_splice_paths {
    my $self = shift;
    my @valid_splice_paths = @_;
    ## completely stomps the existing contents!!!

    @{$self->{_valid_splice_paths}} = @valid_splice_paths;
    return;
}

sub _graph_node_exists {
    my $self = shift;
    my ($type, $lend, $rend, $orient) = @_;
    my $hashkey = $self->_get_hash_key($type, $lend, $rend, $orient);
    return (exists $self->{_graph_node_hashkey_lookup}->{$hashkey});
}

sub _get_graph_node_via_coords_n_type {
    my $self = shift;
    my ($type, $lend, $rend, $orient) = @_;
    my $hashkey = $self->_get_hash_key($type, $lend, $rend, $orient);
    my $node = $self->{_graph_node_hashkey_lookup}->{$hashkey};
    unless ($node) {
        confess "Error, no graph node retrieved based on data ($type, $lend, $rend, $orient)";
    }
    return ($node);
}



####
sub _add_graph_node {
    my $self = shift;
    my ($type, $lend, $rend, $orient) = @_;
    my $hash_key = $self->_get_hash_key($type, $lend, $rend, $orient);
    
    ## only add the node if it doesn't already exist!
    my $graph_node_hashkey_lookup_href = $self->{_graph_node_hashkey_lookup};
    if (my $existing_node = $graph_node_hashkey_lookup_href->{$hash_key}) {
        ## increment evidence for existing node:
        $existing_node->increment_evidence_support();
        return ($existing_node);
    } 
    else {
        # add it
        my $node = Splice_graph_node->new($type, $lend, $rend, $orient);
        push (@{$self->{_graph_nodes}}, $node); # add to complete node list
        ## helpers for node access:
        $graph_node_hashkey_lookup_href->{$hash_key} = $node; # store in lookup table
        my $nodeID = $node->get_nodeID();
        $self->{_graph_node_via_nodeID}->{$nodeID} = $node;
        
        ## splay based on type:
        my $type = $node->get_type();
        
        if ($type eq "internal_exon") {
            push (@{$self->{_internal_exons}}, $node);
        }
        elsif ($type eq "intron") {
            push (@{$self->{_introns}}, $node);
        }
        elsif ($type eq "terminal_left_exon") {
            push (@{$self->{_terminal_left_exons}}, $node);
        }
        elsif ($type eq "terminal_right_exon") {
            push (@{$self->{_terminal_right_exons}}, $node);
        }
        elsif ($type eq "singleton_exon") {
            push (@{$self->{_singleton_exons}}, $node);
        }
        else {
            confess "Error, do not recognize type: $type for node";
        }
        return ($node);
    }
    
}

####
sub _add_intron {
    my $self = shift;
    my ($intron_lend, $intron_rend, $spliced_orient) = @_;
    
    my $node = $self->_add_graph_node("intron", $intron_lend, $intron_rend, $spliced_orient);
    return ($node);
}

####
sub _add_internal_exon {
    my $self = shift;
    my ($exon_lend, $exon_rend, $spliced_orient) = @_;
    my $node = $self->_add_graph_node("internal_exon", $exon_lend, $exon_rend, $spliced_orient);
    return ($node);
}


####
sub get_graph_nodes {
    my $self = shift;
    return (sort {$a->{lend}<=>$b->{lend}} @{$self->{_graph_nodes}});
}




####
sub get_graph_node_via_nodeID {
    my $self = shift;
    my $nodeID = shift;
    my $node_obj = $self->{_graph_node_via_nodeID}->{$nodeID} or confess "Error, no node found based on nodeID: $nodeID\n" . $self->toString();
    return ($node_obj);
}

####
sub toString {
    my $self = shift;
    my @graph_nodes = $self->get_graph_nodes();
    my $num_graph_nodes = scalar (@graph_nodes);
    my $text = "Splice_graph_assembler instance with $num_graph_nodes graph nodes:\n";
    foreach my $graph_node (@graph_nodes) {
        $text .= $graph_node->toString() . "\n";
    }
    
    $text .= "\tvalid paths thru nodes:\n";
    foreach my $splice_path ($self->_get_valid_splice_paths()) {
        $text .= "\t" . $splice_path->toString() . "\n";
    }
    
    $text .= "\tassembled splice paths:\n";
    foreach my $splice_path (@{$self->{_assembled_splice_paths}}) {
        $text .= "\t" . $splice_path->toString() . "\n";
    }
    
    $text .= "\tFinal assemblies with termini\n";
    foreach my $assembly ($self->get_assemblies()) {
        my @node_list = @{$assembly->{__Splice_graph_assembler_nodeID_list}};
        $text .= "\t" . join (",", @node_list) . "\n";
    }

    $text .= $self->toAlignIllustration(60);
    
    return ($text);
    
}

=item toAlignIllustration()

=over 4

B<Description:> illustrates the individual cDNAs to be assembled along with the final products.

B<Parameters:> $max_line_chars(optional)

$max_line_chars is an integer representing the maximum number of characters in a single line of output to the terminal.  The default is  100.

B<Returns:> $alignment_illustration_text

$alignment_illustration_text is a string containing a paragraph of text which illustrates the alignments and assemblies. An example is below:

 --->    <-->  <----->     <--->    <----------------   (+)gi|1199466

 --->    <-->  <----->     <--->    <------------   (+)gi|1209702
                        
---->    <-->  <----    (+)AV827070

---->    <-->  <---     (+)AV828861

---->    <-->  <---      (+)AV830936

 --->    <-->  <-       (+)H36350

ASSEMBLIES: (1)

---->    <-->  <----->     <--->    <---------------- (+) gi|1199466, gi|1209702, AV827070, AV828861, AV830936, H36350




=back

=cut

    ;    

sub toAlignIllustration () {
    my $self = shift;
    my $max_line_chars = shift;
    $max_line_chars = ($max_line_chars) ? $max_line_chars : 100; #if not specified, 100  chars / line is default.
    
    ## Get minimum coord for relative positioning.
    my @coords;
    my @alignments = @{$self->{_incoming_alignments}};
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
    
    my @assemblies = @{$self->{_assemblies}};
    my $num_assemblies = $#assemblies + 1;
    $alignment_text .= "\n\nASSEMBLIES: ($num_assemblies)\n";
    foreach my $assembly (@assemblies) {
        $alignment_text .= "    " . $assembly->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
    }

    if (my @unincorporated_alignments = $self->get_unincorporated_alignments()) {
        my $num_unincorporated = scalar @unincorporated_alignments;
        $alignment_text .= "\n\nUNINCORPORATED_ALIGNMENTS($num_unincorporated)\n";
        foreach my $alignment (@unincorporated_alignments) {
            $alignment_text .= "    " . $alignment->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
        }
    }                                                   
        
    return ($alignment_text);
}

####
sub _get_hash_key {
    my $self = shift;
    my ($type, $lend, $rend, $orient) = @_;
    return ("$type,$lend,$rend,$orient");
}

####
sub _purge_nodes {
    my $self = shift;
    my (%nodeIDs) = @_;
    
    ## perform deletions based on hashkeys
    foreach my $nodeID (keys %nodeIDs) {
        print "PURGING node: $nodeID\n" if $SEE;
        my $node = $self->{_graph_node_via_nodeID}->{$nodeID};
        my $type = $node->get_type();
        my ($lend, $rend) = $node->get_coords();
        my $orient = $node->get_orient();
        
        my $hashkey = $self->_get_hash_key($type, $lend, $rend, $orient);
        delete $self->{_graph_node_hashkey_lookup}->{$hashkey};
        delete $self->{_graph_node_via_nodeID}->{$nodeID};
    }
    
    ## now, do array replacements:
    foreach my $node_list_aref ( $self->{_graph_nodes},
                                 ## only deleting the terminal exons, even though this could be more generic
                                 $self->{_terminal_left_exons},
                                 $self->{_terminal_right_exons},
                                 ) {
        my @replacments;
        my $need_replacement_flag = 0;
        foreach my $node (@$node_list_aref) {
            my $nodeID = $node->get_nodeID();
            if ($nodeIDs{$nodeID}) {
                ## must delete!
                $need_replacement_flag = 1;
            }
            else {
                push (@replacments, $node);
            }
        }
        if ($need_replacement_flag) {
            @$node_list_aref = @replacments; ## Doing Replacment
        }
    }
    
    return;
}

####
sub _right_left_terminal_exons_overlap_via_position_lists {
    my $self = shift;
    my ($right_nodeID, $left_nodeID, $left_splice_coord, $right_splice_coord) = @_;
    
    ## looking for the following:
    #                                 
    #                                                              ## termini of transcripts in right terminal
    #              <---------------------X----------------------   ## right terminal exon
    #        -------------------X--------------------------->      ## left terminal exon
    #                                                              ## termini of transcripts in left terminal
    #                    
    #
    #   The X termini show the endpoints of other transcripts incorporated that would define a proper merging situation.

    # look for boundary pair such that both are included within the splice junctions
    # and right boundary >= left boundary

    my $right_pos_list_aref = $self->{_terminal_exon_nodeID_to_nonsplice_position_list}->{$right_nodeID};
    my $left_pos_list_aref = $self->{_terminal_exon_nodeID_to_nonsplice_position_list}->{$left_nodeID};
    
    foreach my $right_pos (@$right_pos_list_aref) {
        unless ($right_pos > $left_splice_coord && $right_pos < $right_splice_coord) { next; }
        
        foreach my $left_pos (@$left_pos_list_aref) {
            unless ($left_pos > $left_splice_coord && $left_pos < $right_splice_coord) { next; }
            
            if ($right_pos >= $left_pos) {
                return (1); # found suitable case
            }
        }
    }

    return (0); # no such example found.
}

####
sub _analyze_intronless_alignments {
    my $self = shift;
    my @intronless_segments = @_;  ## actually alignment objects.


    print "method: _analyze_intronless_alignments()\n" if $SEE;

    my %applied_segment_indices;
    
    {
        ## hack in a hidden attribute that uniquely identifies each of these segments.
        ## perl allows this, but I'm not so happy with doing it.  for now, it'll be fine.
        ## This __intronless_segment_index will be used to identify those segments that are 
        ## applied to existing structures (internal/terminal exons), and those left over
        ## and still need to be accounted for.
        
        my $id = 0;
        foreach my $intronless_segment (@intronless_segments) {
            $id++;
            $intronless_segment->{__intronless_segment_ID} = $id;
        }
    }
    
    ## increment evidence for internal exons containing intronless segment
    ## and extend terminal exons overlapping intronless segments
    
    ## first, examine the internal_segments:
    foreach my $internal_exon (@{$self->{_internal_exons}}) {
        my $internal_exon_orient = $internal_exon->get_orient();
        my ($internal_exon_lend, $internal_exon_rend) = $internal_exon->get_coords();
        
        

        foreach my $intronless_segment (@intronless_segments) {
            my $intronless_segment_orient = $intronless_segment->get_spliced_orientation();
            my ($intronless_lend, $intronless_rend) = $intronless_segment->get_coords();
            
            if ($intronless_segment_orient eq '?' || $intronless_segment_orient eq $internal_exon_orient) {

                ## look for encapsulation
                if ( ($intronless_lend + $FUZZ_DIST) >= $internal_exon_lend &&
                     ($intronless_rend - $FUZZ_DIST) <= $internal_exon_rend) {
                    
                    $internal_exon->increment_evidence_support();
                    $applied_segment_indices{ $intronless_segment->{__intronless_segment_ID} } = 1; # incorporated already
                    
                    print "-intronless segment $intronless_lend, $intronless_rend, $intronless_segment_orient found incorporated in internal exon: " . $internal_exon->toString() . "\n" if $SEE;
                }
            }
        }
    }
    
    ## examine left terminal exons
    ## first sort so we examine the right boundaries in order from right to left to faciliate extensions
    @intronless_segments = reverse sort {$a->{rend}<=>$b->{rend}} @intronless_segments;

    foreach my $left_terminal_exon (@{$self->{_terminal_left_exons}}) {
        my $left_terminal_exon_orient = $left_terminal_exon->get_orient();
        
        #   ----------------->  # left terminal exon

        foreach my $intronless_segment (@intronless_segments) {
            
            my ($left_terminal_lend, $left_terminal_rend) = $left_terminal_exon->get_coords(); # do it here because it might change below
            my $intronless_segment_orient = $intronless_segment->get_spliced_orientation();
            my ($intronless_lend, $intronless_rend) = $intronless_segment->get_coords();
            
            unless ($intronless_segment_orient eq '?' || $intronless_segment_orient eq $left_terminal_exon_orient) { next; }

            ## must overlap
            unless ($left_terminal_lend <= $intronless_rend && $left_terminal_rend >= $intronless_lend) { next; } 

            unless ($intronless_rend - $FUZZ_DIST <= $left_terminal_rend) { next; }
            ## must be :
            #    ----------------------->  # left terminal exon
            #  ----------------            # intronless segment, overlaps but not passed the splice boundary
            
            # extend lend if intronless segment passes it
            if ($intronless_lend < $left_terminal_lend) {
                # update coords:
                $left_terminal_exon->set_coords($intronless_lend, $left_terminal_rend);
            }
            
            ## if got this far, evidence is incorporated.
            $left_terminal_exon->increment_evidence_support();
            $applied_segment_indices{ $intronless_segment->{__intronless_segment_ID} } = 1; # incorporated
            print "-intronless segment $intronless_lend, $intronless_rend, $intronless_segment_orient incorporated in left terminal exon: " . $left_terminal_exon->toString() . "\n" if $SEE;
        }
    }
            
    ## examine right terminal exons
    ## first sort so we examine the left boundaries in order from left to right to faciliate extensions
    @intronless_segments = sort {$a->{lend}<=>$b->{lend}} @intronless_segments;
    
    foreach my $right_terminal_exon (@{$self->{_terminal_right_exons}}) {
        my $right_terminal_exon_orient = $right_terminal_exon->get_orient();
        
        #   <-----------------  # right terminal exon
        
        foreach my $intronless_segment (@intronless_segments) {
            
            my ($right_terminal_lend, $right_terminal_rend) = $right_terminal_exon->get_coords(); # do it here because it might change below
            my $intronless_segment_orient = $intronless_segment->get_spliced_orientation();
            my ($intronless_lend, $intronless_rend) = $intronless_segment->get_coords();
            
            unless ($intronless_segment_orient eq '?' || $intronless_segment_orient eq $right_terminal_exon_orient) { next; }
            
            ## must overlap
            unless ($right_terminal_lend <= $intronless_rend && $right_terminal_rend >= $intronless_lend) { next; } 
            
            ## must be :
            #    <----------------                    # right terminal exon
            #             ----------------            # intronless segment, overlaps but not passed the splice boundary
            
            unless ($intronless_lend + $FUZZ_DIST >= $right_terminal_lend) { next; }
            
            # extend lend if intronless segment passes it
            if ($intronless_rend > $right_terminal_rend) {
                # update coords:
                $right_terminal_exon->set_coords($right_terminal_lend, $intronless_rend);
            }
            
            ## if got this far, evidence is incorporated.
            $right_terminal_exon->increment_evidence_support();
            $applied_segment_indices{ $intronless_segment->{__intronless_segment_ID} } = 1; # incorporated
            print "-intronless segment $intronless_lend, $intronless_rend, $intronless_segment_orient incorporated in right terminal exon: " . $right_terminal_exon->toString() . "\n" if $SEE;
        }
    }

    ## check to see if there are any unincorporated single segments, and if so, instantiate maximal single segments that contain them.
    my $have_leftover_single_segments_flag = 0;
    foreach my $intronless_segment (@intronless_segments) {
        unless ($applied_segment_indices{ $intronless_segment->{__intronless_segment_ID} }) {
            $have_leftover_single_segments_flag = 1;
            last;
        }
    }
    unless ($have_leftover_single_segments_flag) {
        ## all done; they've all been accounted for.
        print "-All intronless segments are accounted for.\n" if $SEE;
        return;
    }

    ## if got here, we have some single segments that haven't been accounted for.
    print "-some intronless segments are as of yet unincorporated. Need to instantiate singleton exons\n" if $SEE;
    $self->_instantiate_maximal_single_exons(\@intronless_segments, \%applied_segment_indices);
    
    return;
}

####
sub _instantiate_maximal_single_exons {
    my $self = shift;
    my ($intronless_alignments_aref, $applied_segment_indices_href) = @_;
    
    print "\nmethod: _instantiate_maximal_single_exons\n" if $SEE;

    ## get all coords for segments:
    my %segment_id_to_coords;
    my %segment_id_to_alignment;
    my $plus_strand_overlap_assembler = new Overlap_piler();
    my $minus_strand_overlap_assembler = new Overlap_piler();
    my $got_plus_strand_flag = 0;
    my $got_minus_strand_flag = 0;
    
    foreach my $alignment (@$intronless_alignments_aref) {
        print "Got intronless alignment: " . $alignment->toToken() . "\n" if $SEE;
        my ($align_lend, $align_rend) = $alignment->get_coords();
        my $align_id = $alignment->{__intronless_segment_ID};
        $segment_id_to_coords{$align_id} = [$align_lend, $align_rend];
        $segment_id_to_alignment{$align_id} = $alignment;
        
        my $spliced_orient = $alignment->get_spliced_orientation();
        if ($spliced_orient eq '+' || $spliced_orient eq '?') {
            $plus_strand_overlap_assembler->add_coordSet($align_id, $align_lend, $align_rend);
            print "Adding $align_lend-$align_rend [s$spliced_orient] to plus strand assembler\n" if $SEE;
            if ($spliced_orient eq '+') {
                $got_plus_strand_flag = 1;
            }
        }
        if ($spliced_orient eq '-' || $spliced_orient eq '?') {
            $minus_strand_overlap_assembler->add_coordSet($align_id, $align_lend, $align_rend);
            print "Adding $align_lend-$align_rend [s$spliced_orient] to minus strand assembler\n" if $SEE;
            if ($spliced_orient eq '-') {
                $got_minus_strand_flag = 1;
            }
        }
        
    }
    
    ## unless we have some known spliced orientation, then just assemble everything using the plus strand assembler; either could be used w/ no difference.
    unless ($got_plus_strand_flag || $got_minus_strand_flag) {
        $got_plus_strand_flag = 1;
    }
    

    my @singleton_clusters;
    if ($got_plus_strand_flag) {
        print "\nBuilding Plus strand clusters of singletons:\n" if $SEE;
        my @plus_strand_clusters = $plus_strand_overlap_assembler->build_clusters();
        if ($SEE) {
            print "\nSingleton assemblies, PLUS strand assembler:\n";
            @plus_strand_clusters = sort {$a->[0]<=>$b->[0]} @plus_strand_clusters;
            foreach my $cluster (@plus_strand_clusters) {
                print "clustered IDs [+]: " . join ("-", @$cluster) . "\n";
            }
        }
        push (@singleton_clusters, @plus_strand_clusters);
    }
    if ($got_minus_strand_flag) {
        print "\nBuilding Minus strand clusters of singletons:\n" if $SEE;
        my @minus_strand_clusters = $minus_strand_overlap_assembler->build_clusters();
        if ($SEE) {
            print "\nSingleton assemblies, MINUS strand assembler:\n";
            @minus_strand_clusters = sort {$a->[0]<=>$b->[0]} @minus_strand_clusters;
            foreach my $cluster (@minus_strand_clusters) {
                print "clustered IDs [-]: " . join ("-", @$cluster) . "\n";
            }
        }
        push (@singleton_clusters, @minus_strand_clusters);
    }
    
    ## convert to singleton assemblies:
    my @singleton_assemblies;
    foreach my $singleton_cluster (@singleton_clusters) {
        my @ids = @$singleton_cluster;
        my ($min_lend, $max_rend);
        my %aligned_orient_counts;
        my %spliced_orient_counts;
        foreach my $id (@ids) {
            my ($lend, $rend) = @{$segment_id_to_coords{$id}};
            if (! defined ($min_lend)) {
                ($min_lend, $max_rend) = ($lend, $rend);
            }
            else {
                if ($lend < $min_lend) { $min_lend = $lend; }
                if ($rend > $max_rend) { $max_rend = $rend; }
            }
            my $alignment = $segment_id_to_alignment{$id};
            my $spliced_orient = $alignment->get_spliced_orientation();
            my $aligned_orient = $alignment->get_aligned_orientation();
            if ($spliced_orient ne '?') {
                $spliced_orient_counts{$spliced_orient}++;
            }
            $aligned_orient_counts{$aligned_orient}++;
        }
        
        ## get orientation value with most support:
        my @spliced_orients = reverse sort {$spliced_orient_counts{$a}<=>$spliced_orient_counts{$b}} keys %spliced_orient_counts;
        if (scalar @spliced_orients > 1) {
            # should only be one spliced orient, or no spliced orient if everything was '?'
            confess "Error, captured a singleton assembly with multiple spliced orientations.";
        }
        my $assembly_spliced_orient = shift @spliced_orients;

        unless ($assembly_spliced_orient) {
            ## no evidence for assembly spliced orientation.
            ## take the aligned orientation that has the greatest support.
            my @aligned_orients = reverse sort {$aligned_orient_counts{$a}<=>$aligned_orient_counts{$b}} keys %aligned_orient_counts;
            $assembly_spliced_orient = shift @aligned_orients;
        }
        
        unless ($assembly_spliced_orient =~ /^[\+\-]$/) {
            confess "Error, couldn't decide upon assembly orientation for singleton assembly";
        }
                
        push (@singleton_assemblies, { lend => $min_lend,
                                       rend => $max_rend,
                                       ids => [@ids],
                                       length => $max_rend - $min_lend + 1,
                                       num_alignments_included => scalar @ids,
                                       orient => $assembly_spliced_orient,
                                   } 
              );
    }
    
    ## sort so that we maximize for number of alignments included and alignment span
    @singleton_assemblies = reverse sort {$a->{num_alignments_included}<=>$b->{num_alignments_included}
                                          ||
                                              $a->{length}<=>$b->{length}} @singleton_assemblies;
    
    my %ids_accounted_for;
    
  SINGLETON_ASSEMBLY_ANALYSIS:
    foreach my $singleton_assembly (@singleton_assemblies) {
        my ($lend, $rend, $orient, $ids_aref) = ($singleton_assembly->{lend},
                                                 $singleton_assembly->{rend},
                                                 $singleton_assembly->{orient},
                                                 $singleton_assembly->{ids});
        

        ## make sure singleton contains an alignment that wasn't already applied elsewhere:
        my $has_alignment_not_applied_elsewhere_flag = 0;
        foreach my $id (@$ids_aref) {
            unless ($applied_segment_indices_href->{$id}) {
                $has_alignment_not_applied_elsewhere_flag = 1;
                last;
            }
        }
        unless ($has_alignment_not_applied_elsewhere_flag) {
            ## no reason to bother pursuing it.  Everythings accounted for in other already instantiated exons.
            next SINGLETON_ASSEMBLY_ANALYSIS;
        }
        

        my $found_id_not_accounted_for_flag = 0; ## track those IDs that are built into other already instantiated singleton assemblies
        foreach my $id (@$ids_aref) {
            unless ($ids_accounted_for{$id}) {
                $found_id_not_accounted_for_flag = 1;
                last;
            }
        }
        if ($found_id_not_accounted_for_flag) {
            ## add a singleton assembly:
            $self->_add_graph_node("singleton_exon", $lend, $rend, $orient);
            foreach my $id (@$ids_aref) {
                $ids_accounted_for{$id} = 1;
            }
        }
    }


    ## Ensure that all are accounted for now:
    foreach my $intronless_alignment (@$intronless_alignments_aref) {
        my $id = $intronless_alignment->{__intronless_segment_ID};
        if (! $applied_segment_indices_href) {
            ## if not applied before entering this method, then should be applied now
            if (! $ids_accounted_for{$id}) {
                my $unaccounted_for_intronless_alignment = $segment_id_to_alignment{$id};
                confess "Error, intronless alignment " . $unaccounted_for_intronless_alignment->toToken() . "\n"
                    . "is not accounted for after instantiating maximal single exons.\n";
            }
        }
    }
    print "-all intronless alignments should now be accounted for.\n" if $SEE;

    return;
}


####
sub _chain_compatible_splice_paths {
    my $self = shift;
    
    ## wrap splice path objs into structs for scoring and dynamic programming to build chains:
    my @path_structs;
    foreach my $splice_path ($self->_get_valid_splice_paths()) {
        my $pathID = $splice_path->get_pathID();
        my @nodeIDs = $splice_path->get_ordered_nodeIDs();
        my $evidence_support_count = $self->_get_evidence_support_from_nodeIDs(@nodeIDs);
        my $orient = $splice_path->get_orient();
        my ($lend, $rend) = $splice_path->get_coords();
        
        my $struct = {  splice_path_obj => $splice_path,
                        path_score => $evidence_support_count,
                        struct_score => $evidence_support_count,
                        nodeIDs => [@nodeIDs],
                        lend => $lend,
                        rend => $rend,
                        orient => $orient,
                        pathID => $pathID,
                        prev => undef,   # previous link in a chain of compatible structs
                    };
        push (@path_structs, $struct);
    }
    
    @path_structs = sort {$a->{lend}<=>$b->{lend}} @path_structs;

    ## score them:
    for (my $i = 0; $i < $#path_structs; $i++) {

        my $i_path_struct = $path_structs[$i];
        my ($i_lend, $i_rend, $i_orient, $i_splice_path, $i_nodeIDs_aref) = ($i_path_struct->{lend},
                                                                             $i_path_struct->{rend},
                                                                             $i_path_struct->{orient},
                                                                             $i_path_struct->{splice_path_obj},
                                                                             $i_path_struct->{nodeIDs});
        

        for (my $j = $i+1; $j <= $#path_structs; $j++) {
            
            ## No self comparisons:
            if ($i == $j) { next; }

            ## struct j comes after struct i
            my $j_path_struct = $path_structs[$j];
            my ($j_lend, $j_rend, $j_orient, $j_splice_path, $j_nodeIDs_aref) = ($j_path_struct->{lend},
                                                                                 $j_path_struct->{rend},
                                                                                 $j_path_struct->{orient},
                                                                                 $j_path_struct->{splice_path_obj},
                                                                                 $j_path_struct->{nodeIDs});
            
            ## must have same orient:
            unless ($i_orient eq $j_orient) { next; }

            ## check for overlap:
            unless ($i_lend <= $j_rend && $i_rend >= $j_lend) { next; }
            
            ## must be compatible:
            unless ($i_splice_path->is_compatible($j_splice_path)) { next; }
            
            ## double check that neither is a subset of the other:
            if ($i_splice_path->is_subpath_of($j_splice_path) || $j_splice_path->is_subpath_of($i_splice_path)) {
                confess "Error, trying to chain two splice paths where one is a subpath of the other:\n" . $self->toString() 
                    . "\nOffending paths:\n"
                        . $i_splice_path->toString() . "\n"
                            . $j_splice_path->toString() . "\n";
            }
            
            ## calculate path score ending at struct j
            ## don't count the same nodes twice, though:
            my @nodes_found_in_both_paths = $self->_intersection_of_lists($i_nodeIDs_aref, $j_nodeIDs_aref);
            my $score_decrement = $self->_get_evidence_support_from_nodeIDs(@nodes_found_in_both_paths);
            my $candidate_path_score_j = $i_path_struct->{path_score} + $j_path_struct->{struct_score} - $score_decrement;
            if ($candidate_path_score_j > $j_path_struct->{path_score}) {
                ## found best path:
                $j_path_struct->{path_score} = $candidate_path_score_j;
                $j_path_struct->{prev} = $i_path_struct->{pathID};
            }
        }
    }

    ## get the highest scoring paths of compatible subpaths (structs)
    my %pathIDs_left;
    my %pathID_to_struct;
    foreach my $struct (@path_structs) {
        my $pathID = $struct->{pathID};
        $pathIDs_left{$pathID} = 1;
        $pathID_to_struct{$pathID} = $struct;
    }
    
    my @final_assembled_path_node_lists;
    while (%pathIDs_left) {
        ## find the highest scoring chain of paths:
        my $highest_path_score = 0;
        my $highest_scoring_pathID = undef;
        
        foreach my $pathID (keys %pathIDs_left) {
            my $struct = $pathID_to_struct{$pathID};
            my $path_score = $struct->{path_score};
            if ($path_score > $highest_path_score) {
                $highest_path_score = $path_score;
                $highest_scoring_pathID = $pathID;
            }
        }
        
        my %nodes_along_chain;
        ## walk along chain links and collect the nodes:
        my $pathID = $highest_scoring_pathID;
        while (defined $pathID) {
            delete $pathIDs_left{$pathID}; # remove it so we know it's accounted for as a traversed chain link.
            my $struct = $pathID_to_struct{$pathID};
            my $nodes_aref = $struct->{nodeIDs};
            foreach my $node (@$nodes_aref) {
                $nodes_along_chain{$node} = 1;
            }
            $pathID = $struct->{prev};
        }
        
        my @nodes_along_chain_list = keys %nodes_along_chain;
        push (@final_assembled_path_node_lists, [@nodes_along_chain_list]);
    }

    ## Create new path objects for these chains of subpaths:
    foreach my $node_list (@final_assembled_path_node_lists) {
        my @nodes = @$node_list;
        ## sort them according to position:
        @nodes = sort { $self->get_graph_node_via_nodeID($a)->{lend} 
                                       <=>
                        $self->get_graph_node_via_nodeID($b)->{lend}} @nodes;
        my $path_lend = $self->get_graph_node_via_nodeID($nodes[0])->{lend};
        my $path_rend = $self->get_graph_node_via_nodeID($nodes[$#nodes])->{rend};
        my $orient = $self->get_graph_node_via_nodeID($nodes[0])->get_orient();
        
        my $splice_path = Splice_graph_path->new($path_lend, $path_rend, $orient, [@nodes]);
        push (@{$self->{_assembled_splice_paths}}, $splice_path);
    }
    return;
}

####
sub _intersection_of_lists {
    my $self = shift;
    my ($list_A_aref, $list_B_aref) = @_;
    
    my %eles_in_A;
    foreach my $ele (@$list_A_aref) {
        $eles_in_A{$ele} = 1;
    }
    my @intersection;
    foreach my $ele (@$list_B_aref) {
        if ($eles_in_A{$ele}) {
            push (@intersection, $ele);
        }
    }

    return (@intersection);
}


####
sub _get_evidence_support_from_nodeIDs {
    my $self = shift;
    my @nodeIDs = @_;
    my $evidence_sum = 0;
    foreach my $nodeID (@nodeIDs) {
        my $node = $self->get_graph_node_via_nodeID($nodeID);
        my $ev_support = $node->get_num_evidence_support();
        $evidence_sum += $ev_support;
    }
    return ($evidence_sum);
}

####
sub _build_splice_graph {
    my $self = shift;
    
    print "## building splice graph:\n" if $SEE;

    ## do n^2 comparison of nodes:
    ## number of nodes is relatively small, so this isn't too much of a bottleneck
    my @nodes = $self->get_graph_nodes(); # already sorted by lend
    
    ## init base scores for nodes 
    foreach my $node (@nodes) {
        my $ev_support = $node->get_num_evidence_support();
        # init base score
        $node->{forward_base_score} = $ev_support;
        $node->{reverse_base_score} = $ev_support;
        # init path score
        $node->{forward_path_score} = $ev_support;
        $node->{reverse_path_score} = $ev_support;
    }
    
    ## build the splice graph 
    ## at the same time, do the forward dynamic programming calculation so we can backtrack to the best scoring path from any node
    ## nodeB must come after nodeA, be of same orientation, adjacent, and an acceptable linkage (intron to exon)
    ## perform calculation from left to right with backtracking from right to left.
    for (my $i = 1; $i <= $#nodes; $i++) {
        my $nodeB = $nodes[$i];
        my ($nodeB_lend, $nodeB_rend) = $nodeB->get_coords();
        my $nodeB_orient = $nodeB->get_orient();
        my $nodeB_type = $nodeB->get_type();
        my $nodeB_ID = $nodeB->get_nodeID();
        
        for (my $j = $i - 1; $j >= 0; $j--) {
            my $nodeA = $nodes[$j];
            my ($nodeA_lend, $nodeA_rend) = $nodeA->get_coords();
            my $nodeA_orient = $nodeA->get_orient();
            my $nodeA_type = $nodeA->get_type();
            my $nodeA_ID = $nodeA->get_nodeID();
            
            print "-comparing " . $nodeA->toString() . " <=> " . $nodeB->toString() if $SEE;

            unless ($nodeA_orient eq $nodeB_orient) { 
                print "-opposite orients, cannot link.\n" if $SEE;
                next;
            }
            
            unless ($nodeB_lend == $nodeA_rend + 1) { ## next base coordinate for next feature
                print "-not adjacent, cannot link.\n" if $SEE;
                next;
            }

            ## ensure proper connection:
            unless ($ACCEPTABLE_CONNECTIONS{$nodeA_type}->{$nodeB_type}) {
                print "-unacceptable linkage types, cannot link ($nodeA_type,$nodeB_type).\n" if $SEE;
                next;
            }

            print "* linking nodes.\n" if $SEE;
            
            ## if got this far, connections are perfect!
            $nodeA->connect_this_next_nodeID($nodeB_ID);
            $nodeB->connect_this_prev_nodeID($nodeA_ID);

            ## check scores for best scoring path solution
            my $path_score = $nodeA->{forward_path_score} + $nodeB->{forward_base_score};
            if ($path_score > $nodeB->{forward_path_score}) {
                ## link the nodes as current best path linkage
                $nodeB->{forward_path_score} = $path_score;
                $nodeB->{forward_backtrack_nodeID} = $nodeA_ID;
            }
        }
    }

    ## Do the reverse dynamic programming calculation so we can find the best scoring path in a left to right traversal from any node
    ## calculation from right to left with backtracking from left to right
    for (my $i = ($#nodes - 1); $i >= 0; $i--) {
        my $nodeA = $nodes[$i];
        my $nodeA_ID = $nodeA->get_nodeID();

        for (my $j = $i+1; $j <= $#nodes; $j++) {
            my $nodeB = $nodes[$j];
            my $nodeB_ID = $nodeB->get_nodeID();
            print "reverse comparison " . $nodeA->toString() . " <=> " . $nodeB->toString() . "\t" if $SEE;
            ## already stored the prev/next links, so rely on that to know if an acceptable linkage exists:
            if ($nodeA->has_next_nodeID($nodeB_ID)) {
                print "connectable.\n" if $SEE;
                ## check scores for DP
                my $base_score = $nodeA->{reverse_base_score};
                my $path_score = $nodeB->{reverse_path_score} + $base_score;
                if ($path_score > $nodeA->{reverse_path_score}) {
                    ## make it the best path and link the node
                    $nodeA->{reverse_path_score} = $path_score;
                    $nodeA->{reverse_backtrack_nodeID} = $nodeB_ID;
                }
            }
            else {
                print "non-connectable nodes.\n" if $SEE;
            }
        }
    }

    ## defensive programming:
    foreach my $node (@nodes) {
        
        ## ignore singletons
        if ($node->get_type() eq "singleton_exon") { next; }
        
        ## ensure that only terminal left exons have forward_backtrack_nodeID's as undef
        if ( (! defined $node->{forward_backtrack_nodeID}) && $node->get_type() ne "terminal_left_exon") {
            print $self->toString();
            print Dumper ($node);
            confess "Error, node has undefined forward backtrack ID but is not a terinal left exon: " . $node->toString() . "\n";
        }
        ## ensure that only terminal right exons have the reverse_backtrack_nodeID's as undef
        if ( (! defined $node->{reverse_backtrack_nodeID}) && $node->get_type() ne "terminal_right_exon") {
            print $self->toString();
            print Dumper ($node);
            confess "Error, node has undefined reverse backtrack ID but is not a terminal right exon: " . $node->toString() . "\n";
        } 
    }
    
    return;
}


####
sub _extend_splice_paths_to_termini {
    my $self = shift;
    
    my @splice_paths = $self->_get_valid_splice_paths();
    foreach my $splice_path (@splice_paths) {
        my @ordered_nodeIDs = $splice_path->get_ordered_nodeIDs();
        
        my $left_nodeID = $ordered_nodeIDs[0];
        my $right_nodeID = $ordered_nodeIDs[$#ordered_nodeIDs];

        print "-extending splice_path [" . join (",", @ordered_nodeIDs) . "] maximally to the left and right:\n" if $SEE;
        my @left_extension_nodes = $self->_extend_node_maximally_left($left_nodeID);
        @left_extension_nodes = grep { $_ != $left_nodeID } @left_extension_nodes; # remove the nucleating node

        my @right_extension_nodes = $self->_extend_node_maximally_right($right_nodeID);
        @right_extension_nodes = grep { $_ != $right_nodeID } @right_extension_nodes;
        
        unless (@left_extension_nodes && @right_extension_nodes) {
            confess "Error, couldn't extend splice path to the left or right to include termini\n"
                . "left: @left_extension_nodes\n"
                    . "right: @right_extension_nodes\n "
                        . $splice_path->toString();
        }

        ## create alignment assembly:
        my @all_nodes_in_maximal_path = (@left_extension_nodes, @ordered_nodeIDs, @right_extension_nodes);
        my $cdna_alignment = $self->_instantiate_cdna_alignment_from_nodeIDs(@all_nodes_in_maximal_path);

        $self->_add_alignment_assembly($cdna_alignment);
    }
    return;
}

####
sub _extend_node_maximally_left {
    my $self = shift;
    my (@nodeIDs) = @_;
    my $highest_scoring_path = $self->_perform_traceback("forward_backtrack", \@nodeIDs);
    my $highest_scoring_nodelist_aref = $highest_scoring_path->{path_nodeIDs_aref};
    return (@$highest_scoring_nodelist_aref);
}


####
sub _extend_node_maximally_right {
    my $self = shift;
    my (@nodeIDs) = @_;
    my $highest_scoring_path = $self->_perform_traceback("reverse_backtrack", \@nodeIDs); 
    my $highest_scoring_nodelist_aref = $highest_scoring_path->{path_nodeIDs_aref};
    return (@$highest_scoring_nodelist_aref);
}


####
sub _perform_traceback {
    my $self = shift;
    my ($traceback_direction, $nodeIDs_aref) = @_;

    my @paths_and_score_structs;
    
    foreach my $nodeID (@$nodeIDs_aref) {
        ## do backtracking from forward DP calculation:
        my @nodes_in_path;
        my $next_path_node = $nodeID;
        my $path_score = 0;
        while (defined $next_path_node) {
            my $node_obj = $self->get_graph_node_via_nodeID($next_path_node);
            my $ev_support = $node_obj->get_num_evidence_support();
            $path_score += $ev_support;
            push (@nodes_in_path, $next_path_node);
            if ($traceback_direction eq "forward_backtrack") {
                $next_path_node = $node_obj->{forward_backtrack_nodeID};
            }
            elsif ($traceback_direction eq "reverse_backtrack") {
                $next_path_node = $node_obj->{reverse_backtrack_nodeID};
            }
            else {
                confess "Don't understand traceback direction: $traceback_direction ";
            }
        }
        my $path_and_score_struct = { score => $path_score,
                                      path_nodeIDs_aref => [@nodes_in_path],
                                  };
        push (@paths_and_score_structs, $path_and_score_struct);
        
    }
    
    @paths_and_score_structs = reverse sort {$a->{score}<=>$b->{score}} @paths_and_score_structs;
    
    my $highest_scoring_path = shift @paths_and_score_structs;
    
    return ($highest_scoring_path);
    
}

=tree_traversal_replaced_by_DP_method


####
my $path_no = 0;
sub _tree_traversal {
    my $self = shift;
    my ($node_obj, $traversal_direction, $current_path_score, $current_path_list_aref, 
        $max_score_sref, $max_path_href) = @_;
    
    #print $self->toString();
    #print $self->toAlignIllustration(60);
    
    # print "_tree_traversal() current_path: ". join (",", @$current_path_list_aref) . "\n" if $SEE;
    
    $self->_ensure_unique_path_nodes(@$current_path_list_aref); # defensive programming

    my $edge_key = "_" . $traversal_direction . "_nodeID_href";
        
    my %edges = %{$node_obj->{$edge_key}};
    
    ## add current node evidence score to the sum score
    
    
    if (%edges) {
        ## not a terminal node
        ## explore each of the linked nodes:
        my @linked_nodes = keys %edges;
        foreach my $linked_node (@linked_nodes) {
            my $linked_node_obj = $self->get_graph_node_via_nodeID($linked_node);
            ## add linked node to path list:
            my @current_path_list = @$current_path_list_aref;
            push (@current_path_list, $linked_node);        
            my $num_evidence_support = $linked_node_obj->get_num_evidence_support();
            my $linked_node_path_score = $current_path_score + $num_evidence_support;    
                  
            $self->_tree_traversal($linked_node_obj, $traversal_direction, $linked_node_path_score, [@current_path_list],  
                                   $max_score_sref, $max_path_href);
        }
        
    }
    else {
        ## a terminal node
        ## base case!
        if ($node_obj->get_type() !~ /terminal/) {
            ## supposed to be a terminal node.
            confess "Reached base case but not a terminal node! " . $node_obj->toString() . "\n" . $self->toString();
        }
        
        $path_no++;
        print "$path_no _tree_traversal($traversal_direction) terminal_reached. PATH: " . join (",", @$current_path_list_aref) . "\n" if $SEE;
        if ($current_path_score >= $$max_score_sref) { #includes case where incoming (starting) node is in fact a terminal node.
            ## make this maximum scoring path:
            $max_path_href->{score} = $current_path_score;
            $max_path_href->{path_nodeIDs_aref} = [@$current_path_list_aref];
            $$max_score_sref = $current_path_score;
        }
    }
    return;
    
}

=cut

####
sub _ensure_unique_path_nodes {
    my $self = shift;
    my @nodes = @_;
    my %check;
    foreach my $node (@nodes) {
        if ($check{$node}) {
            confess $self->toString() . "Error, path list: " .  join (",", @nodes) . " has redundant node: $node\n";
        }
        $check{$node} = 1; # log it!
    }
    return; # all good!
}


####
sub _append_singletons_to_assembly_list {
    my $self = shift;
    foreach my $node_obj (@{$self->{_singleton_exons}}) {
        my $nodeID = $node_obj->get_nodeID();
        my $cdna_alignment = $self->_instantiate_cdna_alignment_from_nodeIDs($nodeID);
        $self->_add_alignment_assembly($cdna_alignment);
    }
    return;
}



####
sub _instantiate_cdna_alignment_from_nodeIDs {
    my $self = shift;
    my @nodeIDs = @_;

    my $spliced_orientation = undef;

    ## add alignment segments for all non-intron nodes:
    my @coordsets;
    foreach my $nodeID (@nodeIDs) {
        my $node_obj = $self->get_graph_node_via_nodeID($nodeID);
        my $node_type = $node_obj->get_type();
        if ($node_type eq "intron") { next; } ## no introns!

        my ($lend, $rend) = $node_obj->get_coords();
        my $orient = $node_obj->get_orient();
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        push (@coordsets, [$end5, $end3]);
        
        unless (defined $spliced_orientation) {
            $spliced_orientation = $orient;
        }

        if ($spliced_orientation ne $orient) {
            confess "Error, alignment segments of opposite orientation are in the same extended splice path! @nodeIDs\n" . $self->toString();
        }
    }
    
    @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;
    my $match_coord_sum = 0;
    my @alignment_segments;
    foreach my $coordset (@coordsets) {
        my ($end5, $end3) = @$coordset;
        my $seg_length = abs ($end3 - $end5) + 1;
        my $match_lend = $match_coord_sum + 1;
        my $match_rend = $match_coord_sum + $seg_length;
        
        my $alignment_segment = new CDNA::Alignment_segment($end5, $end3, $match_lend, $match_rend, 100); # 100% identity
        
        push (@alignment_segments, $alignment_segment);

        $match_coord_sum += $seg_length;
    }

    my $cdna_alignment = new CDNA::CDNA_alignment($match_coord_sum, \@alignment_segments);

    $cdna_alignment->set_spliced_orientation($spliced_orientation);
    $cdna_alignment->force_spliced_validation($spliced_orientation);

    ## hack in the node list as a hidden object attribute.
    $cdna_alignment->{__Splice_graph_assembler_nodeID_list} = [@nodeIDs];
    
    return ($cdna_alignment);
}

####
sub _correlate_assemblies_with_incoming_alignments {
    my $self = shift;
    my @assemblies = $self->get_assemblies();
    
    my %alignment_accs; # keys to alignment objects
    my %incorporated_alignments;
    foreach my $assembly (@assemblies) {

        my $num_segments_in_assembly = $assembly->get_num_segments();
        my @incorporated_alignments;
        
        foreach my $alignment (@{$self->{_incoming_alignments}}) {
            my $alignment_acc = $alignment->get_acc();
            
            ## store alignment object in lookup table for later retrieval.
            unless (exists $alignment_accs{$alignment_acc}) {
                $alignment_accs{$alignment_acc} = $alignment;
            }
        
            ## check for alignment incorporation:
            if ($assembly->encapsulates($alignment, $FUZZ_DIST) &&
                ( 
                  ## since assemblies have fixed spliced orient, compat will fail when encapsulation 
                  ## of single exon alignments of ambiguous orientation exists
                  ($num_segments_in_assembly == 1 && $alignment->get_spliced_orientation() eq '?')  
                  ||
                  ## do rigorous compatibility check:
                  ($assembly->is_compatible($alignment, $FUZZ_DIST)) 
                  )
                ) 
            {
                push (@incorporated_alignments, $alignment_acc);
                $incorporated_alignments{$alignment_acc} = 1;
            }
        }
        $assembly->set_acc( join ("/", @incorporated_alignments)); # use / as DELIMETER for cdna accessions
    }
    
    my @unincorporated_alignments;
    foreach my $alignment_acc (keys %alignment_accs) {
        unless ($incorporated_alignments{$alignment_acc}) {
            push (@unincorporated_alignments, $alignment_accs{$alignment_acc});
        }
    }
    

    @{$self->{_unincorporated_alignments}} = (); #clear
    
    if (@unincorporated_alignments) {
        @{$self->{_unincorporated_alignments}} = @unincorporated_alignments;
        return (1); # still have unincorporated alignments!
    }
    
    return(0); ## all alignments are accounted for. 
}

####
sub _explore_assemblies_from_unincorporated_alignments {
    my $self = shift;
    print "method _explore_assemblies_from_unincorporated_alignments()\n" if $SEE;
    ## decompose unincorporated alignments and create maximal paths that contain them.
    
    my @added_alignment_assemblies;


    my @unincorporated_alignments = $self->get_unincorporated_alignments();

    ## examine in order of descending numbers of segments and length:
    @unincorporated_alignments = reverse sort { $a->{num_segments} <=> $b->{num_segments}
                                         ||
                                             $a->{cdna_length} <=> $b->{cdna_length} } @unincorporated_alignments;
    
    

    foreach my $unincorporated_alignment (@unincorporated_alignments) {
        print "\nUnincorporated alignment under scrutiny: " . $unincorporated_alignment->toToken() . "\n" if $SEE;
        
        my $orient = $unincorporated_alignment->get_spliced_orientation();
        my ($alignment_lend, $alignment_rend) = $unincorporated_alignment->get_coords();

        my $found_in_new_assembly_flag = 0;
        foreach my $new_assembly (@added_alignment_assemblies) {
            if ($new_assembly->is_compatible($unincorporated_alignment, $FUZZ_DIST) && $new_assembly->encapsulates($unincorporated_alignment, $FUZZ_DIST)) {
                $found_in_new_assembly_flag = 1;
                last;
            }
        }
        if ($found_in_new_assembly_flag) {
            print "-found in a newly created assembly based on an earlier unincorporated alignment.  Now accounted for.\n" if $SEE;
            next;
        }

        my $num_segments = $unincorporated_alignment->get_num_segments();
        
        my @all_nodes;
        
        if ($num_segments == 1) {
            ## must be included in some terminal exon not already represented by maximal splice paths (only explanation)
            ## find a maximal path that includes this segment stemming from a terminal exon
            
            my @terminal_nodes = $self->_find_terminal_exons_encompassing_segment($alignment_lend, $alignment_rend, $orient);
            unless (@terminal_nodes) {
                confess "Error, searching for terminal nodes that encompass unincorproated singleton segment and none found.\n" 
                    . $self->toString();
            }
            ## find maximal path:
            my @paths_from_terminal_node;
            foreach my $terminal_node (@terminal_nodes) {
                my $nodeID = $terminal_node->get_nodeID();
                my $type = $terminal_node->get_type();
                my $max_path_href;
                if ($type eq "terminal_left_exon") {
                    ## must extend to the right
                    $max_path_href = $self->_perform_traceback("forward_backtrack", [$nodeID]); #_extend_node_maximally_right($nodeID);
                }
                elsif ($type eq "terminal_right_exon") {
                    $max_path_href = $self->_perform_traceback("reverse_backtrack", [$nodeID]); #_extend_node_maximally_left($nodeID);
                }
                else {
                    confess "Error, terminal node is supposed to be left or right, but type($type) is not recognized.";
                }
                
                push (@paths_from_terminal_node, $max_path_href);
            }
            ## get highest scoring path:
            @paths_from_terminal_node = reverse sort { $a->{score}<=>$b->{score} } @paths_from_terminal_node;
            my $highest_scoring_path = shift @paths_from_terminal_node;
            my $node_IDs_in_path_aref = $highest_scoring_path->{path_nodeIDs_aref};
            @all_nodes = @$node_IDs_in_path_aref;
        }
        
        else {
            ## has an intron! Some alternate termini needed that were not found in a maximal splice path
            
            my @segments = $unincorporated_alignment->get_alignment_segments();
            my $first_segment = shift @segments;
            my $last_segment = pop @segments;
            
            ## create node list for internal segments and introns:
            my @nodes_in_central_path;
            my @introns = $unincorporated_alignment->get_intron_coords();
            foreach my $intron (@introns) {
                my ($lend, $rend) = @$intron;
                my $node_obj = $self->_get_graph_node_via_coords_n_type("intron", $lend, $rend, $orient);
                my $nodeID = $node_obj->get_nodeID();
                push (@nodes_in_central_path, $nodeID);
            }
            foreach my $internal_segment (@segments) {
                my ($lend, $rend) = $internal_segment->get_coords();
                my $node_obj = $self->_get_graph_node_via_coords_n_type("internal_exon", $lend, $rend, $orient);
                my $nodeID = $node_obj->get_nodeID();
                push (@nodes_in_central_path, $nodeID);
            }
            my ($first_segment_lend, $first_segment_rend) = $first_segment->get_coords();
            my @candidate_preceding_nodes = $self->_find_candidate_preceding_exons($first_segment_lend, $first_segment_rend, $orient);
            unless (@candidate_preceding_nodes) {
                confess "Error, no candidate preceding nodes for ($first_segment_lend, $first_segment_rend, $orient)";
            }
            
            my ($last_segment_lend, $last_segment_rend) = $last_segment->get_coords();
            my @candidate_following_nodes = $self->_find_candidate_following_exons($last_segment_lend, $last_segment_rend, $orient);
            unless (@candidate_following_nodes) {
                confess "Error, no candidate following nodes for ($last_segment_lend, $last_segment_rend, $orient)";
            }
            
            ## find maximal paths left and right:
            my @nodeIDs_left;
            foreach my $candidate_preceding_node (@candidate_preceding_nodes) {
                print "Candidate preceding node: " . $candidate_preceding_node->toString() . "\n" if $SEE;
                my $nodeID = $candidate_preceding_node->get_nodeID();
                push (@nodeIDs_left, $nodeID);
            }
            my @path_left = $self->_extend_node_maximally_left(@nodeIDs_left);
            
            my @nodeIDs_right;
            foreach my $candidate_following_node (@candidate_following_nodes) {
                print "Candidate following node: " . $candidate_following_node->toString() . "\n" if $SEE;
                my $nodeID = $candidate_following_node->get_nodeID();
                push (@nodeIDs_right, $nodeID);
            }
            my @path_right = $self->_extend_node_maximally_right(@nodeIDs_right);
            
            @all_nodes = (@path_left, @nodes_in_central_path, @path_right);
        }

        
        my $assembly = $self->_instantiate_cdna_alignment_from_nodeIDs(@all_nodes);

        ## verify that this assembly accounts for the unincorporated alignment:
        unless ( (my $compatible_result = $assembly->is_compatible($unincorporated_alignment, $FUZZ_DIST))
                 && 
                 (my $encapsulated_result = $assembly->encapsulates($unincorporated_alignment, $FUZZ_DIST)) ) {
            confess "Error, created new assembly to house missing alignment, but it doesn't afterall!\n"
                . "unincorp: " . $unincorporated_alignment->toToken() . "\n"
                . "newAsmb: " . $assembly->toToken() . "\n"
                . "compatible: $compatible_result\n"
                . "encapsulated: $encapsulated_result\n";
        }
        print "-verified incorporation in new assembly.\n" if $SEE;
        push (@added_alignment_assemblies, $assembly);
    }
    
    ## add to assembly list
    $self->_add_alignment_assembly(@added_alignment_assemblies);
    
    return;
}

####
sub _find_candidate_preceding_exons {
    my $self = shift;
    my ($lend, $rend, $orient) = @_;
    
    ## must have same orient, share rend boundary, and be internal or terminal exons that contain it
    my @candidate_exons;
    foreach my $exon (@{$self->{_internal_exons}}, @{$self->{_terminal_left_exons}}) {
        my ($exon_lend, $exon_rend) = $exon->get_coords();
        my $exon_orient = $exon->get_orient();
        if ($exon_orient eq $orient && $lend + $FUZZ_DIST >= $exon_lend && $exon_rend == $rend) {
            push (@candidate_exons, $exon);
        }
    }
    
    return (@candidate_exons);
}

####
####
sub _find_candidate_following_exons {
    my $self = shift;
    my ($lend, $rend, $orient) = @_;
    
    ## must have same orient, share rend boundary, and be internal or terminal exons that contain it
    my @candidate_exons;
    foreach my $exon (@{$self->{_internal_exons}}, @{$self->{_terminal_right_exons}}) {
        my ($exon_lend, $exon_rend) = $exon->get_coords();
        my $exon_orient = $exon->get_orient();
        if ($exon_orient eq $orient && $rend - $FUZZ_DIST <= $exon_rend && $exon_lend == $lend) {
            push (@candidate_exons, $exon);
        }
    }
    
    return (@candidate_exons);
}







######################################################
######################################################

package Splice_graph_node;

use strict;
use warnings;
use Carp;

my $nodeID = 0; ## class attribute

####
sub new {
    my $packagename = shift;
    
    my ($type, $lend, $rend, $orient) = @_;
    
    ## type can be:  intron | internal_exon | terminal_left_exon | terminal_right_exon | singleton_exon
    $nodeID++;
    
    # object atts:
    my $self = {
        nodeID => $nodeID,
        type => undef,
        lend => undef,
        rend => undef,
        orient => undef,
        length => undef,
        num_evidence_support => 1,  # indicate number of transcripts supporting this feature
        
        ## graph connections:
        _prev_nodeID_href => {}, ## keys nodeIDs connected to before this node
        _next_nodeID_href => {}, ## keys nodeIDs connected to after this node
    
        ## attributes for dynamic programming calculations to find longest path
        forward_base_score => 0, ## calculation performed from left to right, with right to left backtracking
        forward_path_score => 0,
        forward_backtrack_nodeID => undef,
        
        reverse_base_score => 0, ## calculation performed from right to left with backtracking from left to right
        reverse_path_score => 0,
        reverse_backtrack_nodeID => undef,
        
    };
    
    bless ($self, $packagename);

    ## use set methods to set attribute values and validate settings:
    $self->set_coords($lend, $rend);
    $self->set_type($type);
    $self->set_orient($orient);
    
    return ($self);
}

####
sub set_num_evidence_support {
    my $self = shift;
    my ($evidence_support) = @_;
    
    $self->{num_evidence_support} = $evidence_support;
    return;
}

####
sub increment_evidence_support {
    my $self = shift;
    my ($inc_val) = @_;
    if (defined $inc_val) {
        $self->{num_evidence_support} += $inc_val;
    }
    else {
        ## just add one.
        $self->{num_evidence_support}++;
    }
    return;
}

####
sub get_num_evidence_support {
    my $self = shift;
    return ($self->{num_evidence_support});
}

####
sub toString {
    my $self = shift;
    my $text = "node: " . $self->get_nodeID()
        . "\t" . join ("-", $self->get_coords())
            . "\torient(" . $self->get_orient() . ")"
                . "\t" . $self->get_type()
                    . "\tEvSupport: " . $self->get_num_evidence_support();
    
    return ($text);
}

####
sub get_nodeID {
    my $self = shift;
    return ($self->{nodeID});
}


####
sub set_coords {
    my $self = shift;
    my ($lend, $rend) = @_;

    unless ($lend =~ /^\d+$/ && $rend =~ /^\d+$/) {
        confess "Error, coordinates [$lend] or [$rend] are not integers";
    }
    
    unless ($lend <= $rend) {
        confess "Error, coordinates out of order";
    }
    
    $self->{lend} = $lend;
    $self->{rend} = $rend;
    
    $self->{length} = $rend - $lend + 1;
    
    return;
}

####
sub get_coords {
    my $self = shift;
    return ($self->{lend}, $self->{rend});
}

####
sub set_type {
    my $self = shift;
    my ($type) = @_;
    unless ($type =~ /^(intron|internal_exon|terminal_left_exon|terminal_right_exon|singleton_exon)$/) {
        confess "type $type is not acceptible";
    }
    
    $self->{type} = $type;
    return;
}

####
sub get_type {
    my $self = shift;
    return ($self->{type});
}

####
sub set_orient {
    my $self = shift;
    my ($orient) = @_;
    unless ($orient =~ /^[\+\-]$/) {
        confess "Error, orient($orient) is not allowed here";
    }
    
    $self->{orient} = $orient;
    return;
}

####
sub get_orient {
    my $self = shift;
    return ($self->{orient});
}


####
sub connect_this_prev_nodeID {
    my $self = shift;
    my ($nodeID) = @_;
    $self->{_prev_nodeID_href}->{$nodeID} = 1;
    return;
}

####
sub connect_this_next_nodeID {
    my $self = shift;
    my ($nodeID) = @_;
    $self->{_next_nodeID_href}->{$nodeID} = 1;
    return;
}

####
sub get_connected_prev_nodeIDs {
    my $self = shift;
    my @prev_nodeIDs = keys %{$self->{_prev_nodeID_href}};
    return (@prev_nodeIDs);
}

####
sub get_connected_next_nodeIDs {
    my $self = shift;
    my @next_nodeIDs = keys %{$self->{_next_nodeID_href}};
    return (@next_nodeIDs);
}

####
sub has_next_nodeID {
    my $self = shift;
    my ($nodeID) = @_;
    if ($self->{_next_nodeID_href}->{$nodeID}) {
        return (1); #yes
    }
    else {
        return (0); # no
    }
}

####
sub has_prev_nodeID {
    my $self = shift;
    my ($nodeID) = @_;
    if ($self->{_prev_nodeID_href}->{$nodeID}) {
        return (1); #yes
    }
    else {
        return (0); #no
    }
}

#############################################
#############################################

package Splice_graph_path;

use strict;
use warnings;
use Carp;

my $PATH_ID = 0;  ## class variable

####
sub new {
    my $packagename = shift;
    
    my ($lend, $rend, $orient, $ordered_nodeIDs_aref) = @_;
    
    unless ($lend =~ /^\d+$/ && $rend =~ /^\d+$/) {
        confess "Error, lend $lend and rend $rend must be integers";
    }
    unless ($orient =~ /^[\+\-]$/) {
        confess "Error, orient($orient) is invalid";
    }
    
    unless (@$ordered_nodeIDs_aref) {
        confess "Error, need list of ordered nodeIDs as constructor param";
    }
    
    $PATH_ID++;

    my $self = {
        _lend => $lend,
        _rend => $rend,
        _orient => $orient,
        _ordered_nodeIDs => [@$ordered_nodeIDs_aref], # a path is an ordered list of intron and exon nodeIDs traversed  
        _pathID => $PATH_ID,
        
    };
    
    bless ($self, $packagename);

    return ($self);
}

####
sub get_ordered_nodeIDs {
    my $self = shift;
    return (@{$self->{_ordered_nodeIDs}});
}

####
sub toString {
    my $self = shift;
    my @ordered_nodeIDs = $self->get_ordered_nodeIDs();
    return (join (",", @ordered_nodeIDs));
}

####
sub is_subpath_of {
    my $self = shift;
    my ($other_path) = @_;

    my $self_path_string = "," . $self->toString() . ","; ## add comma terminal delimiters

    my $other_path_string = "," . $other_path->toString() . ","; 

    if ($other_path_string =~ /$self_path_string/) {
        return (1); #true
    }
    else {
        return (0); #false
    }
}

####
sub is_compatible () {
    my $self = shift;
    my ($other_splice_path_obj) = @_;
    
    ## make sure they have at least one node in common:
    my @node_IDs_A = $self->get_ordered_nodeIDs();
    my @node_IDs_B = $other_splice_path_obj->get_ordered_nodeIDs();

    my $found_nodeID_in_common = 0;
    my $common_nodeID = undef;
    { # do this in a more restricted scope
        my %nodeIDs;
        foreach my $nodeID (@node_IDs_A) {
            $nodeIDs{$nodeID} = 1;
        }
        foreach my $nodeID (@node_IDs_B) {
            if ($nodeIDs{$nodeID}) {
                $found_nodeID_in_common = 1;
                $common_nodeID = $nodeID;
                last;
            }
        }
    }
    unless ($found_nodeID_in_common) {
        return (0); # false; no common node so definitely incompatible
    }

    ## examine all nodes adjacent to common node; they should be completely identical
    my $common_node_pos_A = $self->_find_ele_index_pos_via_value(\@node_IDs_A, $common_nodeID);
    my $common_node_pos_B = $self->_find_ele_index_pos_via_value(\@node_IDs_B, $common_nodeID);

    ## search left:
    my ($i, $j);
    for ($i=$common_node_pos_A, $j=$common_node_pos_B; $i >= 0 && $j >= 0; $i--, $j--) {
        my $A_node = $node_IDs_A[$i];
        my $B_node = $node_IDs_B[$j];
        if ($A_node ne $B_node) {
            # incompatible
            return (0); 
        }
    }
    
    ## search right:
    for ($i = $common_node_pos_A, $j = $common_node_pos_B; $i <= $#node_IDs_A && $j <= $#node_IDs_B; $i++,$j++) {
        my $A_node = $node_IDs_A[$i];
        my $B_node = $node_IDs_B[$j];
        if ($A_node ne $B_node) {
            # incompatible
            return (0); 
        }
    }

    ## if got here, then passed compatibility tests.
    return (1); # Yes, compatible
    
}


sub _find_ele_index_pos_via_value {
    my $self = shift;
    my ($list_aref, $value) = @_;
    for (my $i = 0; $i <= $#$list_aref; $i++) {
        if ($list_aref->[$i] eq $value) {
            return ($i);
        }
    }
    confess "Error finding index position of $value in list";
}
                             
####
sub get_pathID {
    my $self = shift;
    return ($self->{_pathID});
}

####
sub get_coords {
    my $self = shift;
    return ($self->{_lend}, $self->{_rend});
}

####
sub get_orient {
    my $self = shift;
    return ($self->{_orient});
}


1; #EOM
