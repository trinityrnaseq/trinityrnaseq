package SegmentGraph;

use strict;
use warnings;
use Carp;
use Data::Dumper;

sub new {
    my $packagename = shift;
    

    my $self = {  id_to_node => {}, # primary container for all graph elements
                  
                  _next_edges => {}, # id => { id => 1 }
                  _prev_edges => {}, 
    };

    bless ($self, $packagename);

    return($self);
}

####
sub delete_node {
    my $self = shift;
    my ($node) = @_;
    
    my $node_ID = $node->get_ID();
   
    my @next_nodes = $self->get_next_nodes($node);
    foreach my $next_node (@next_nodes) {
        $self->_delete_prev_edge_to($next_node, $node_ID);
    }
    my @prev_nodes = $self->get_prev_nodes($node);
    foreach my $prev_node (@prev_nodes) {
        $self->_delete_next_edge_to($prev_node, $node_ID);
    }
    
    delete( $self->{id_to_node}->{$node_ID} );

    return;
}

sub get_next_nodes {
    my $self = shift;
    my ($node) = @_;

    my $node_ID = $node->get_ID();

    my @next_nodes;

    if (my $href = $self->{_next_edges}->{$node_ID}) {
        my @ids = keys %$href;
        foreach my $id (@ids) {
            my $next_node = $self->{id_to_node}->{$id} or confess "Error, no node found for ID: $id";
            push (@next_nodes, $next_node);
        }
    }
    
    return(@next_nodes);
}

sub get_prev_nodes {
    my $self = shift;
    my ($node) = @_;

    my $node_ID = $node->get_ID();

    my @prev_nodes;

    if (my $href = $self->{_prev_edges}->{$node_ID}) {
        my @ids = keys %$href;
        foreach my $id (@ids) {
            my $next_node = $self->{id_to_node}->{$id} or confess "Error, no node found for ID: $id";
            push (@prev_nodes, $next_node);
        }
    }
    
    return(@prev_nodes);
}



sub add_segment {
    my $self = shift;
    my ($lend, $rend, $parent_feature_name) = @_;



    unless ($lend && $rend && $parent_feature_name) {
        confess "Error, need params (lend, rend, parent_feature_name)";
    }


    #print "ADDING SEGMENT: $lend,$rend\n";

    my $overlapping_node = $self->find_overlapping_segment($lend, $rend);
    
    if ($overlapping_node) {
        
        my ($overlapping_node_lend, $overlapping_node_rend) = $overlapping_node->get_coords();
        
        #print "Found overlapping nodes: ($lend,$rend) to ($overlapping_node_lend, $overlapping_node_rend)\n";

        ## number of things could happen here...

        if ($overlapping_node->has_same_coordinates($lend, $rend)) {
            $overlapping_node->add_owners($parent_feature_name);
        }
        elsif ($overlapping_node->overlaps_coordinates($lend, $rend)) {
            
            ## fracture into segments based on overlaps
            
            my @coords = sort {$a<=>$b} ($lend, $rend, $overlapping_node->get_coords());
            
            ## bounds stay the same, but internal set needs to be adjusted.
            my @fractured_segments = ( [$coords[0], $coords[1]-1],
                                       [$coords[1], $coords[2]],
                                       [$coords[2]+1, $coords[3]] );
            
            my @remaining_input_segments;
            my @contained_by_both;
            my @contained_by_overlapping_segment_only;
            foreach my $fractured_segment (@fractured_segments) {
                my ($seg_lend, $seg_rend) = @$fractured_segment;
                
                if ($seg_lend > $seg_rend) { next; }
                
                #print "Searching fragment: $seg_lend,$seg_rend\n";

                # corresponds to the input and the overlapping segments
                if ($overlapping_node->envelops_coordinates($seg_lend, $seg_rend)
                    &&
                    ($seg_lend >= $lend && $seg_rend <= $rend)
                    ) {
                    #print "fragment contained by segment: $seg_lend, $seg_rend\n";
                    push (@contained_by_both, $fractured_segment);
                }
                
                elsif ($overlapping_node->envelops_coordinates($seg_lend, $seg_rend)) {
                    ## a part of the original overlapping segment
                    push (@contained_by_overlapping_segment_only, $fractured_segment);
                }


                # just the input segment
                elsif ($seg_lend >= $lend && $seg_rend <= $rend) {
                    
                    #print "adding to remaining segment: $seg_lend, $seg_rend\n";
                    push (@remaining_input_segments, $fractured_segment);
                }
                else {
                    confess "Error, ended up with coordinate segment that isnt placed: $seg_lend,$seg_rend ";
                }
            }
            
            my @owners = $overlapping_node->get_owners();
            
            $self->delete_node($overlapping_node);
            
            #print "Contained by both: " . Dumper(\@contained_by_both) . "\n";
            #print "Contained by overlapping segment only: " . Dumper(\@contained_by_overlapping_segment_only). "\n";
            #print "Remaining input segments: " . Dumper(\@remaining_input_segments) . "\n";
            

            foreach my $seg_coords_aref (@contained_by_both) {
                my ($seg_lend, $seg_rend) = @$seg_coords_aref;
                $self->_add_segment_node($seg_lend, $seg_rend, [@owners, $parent_feature_name]);
            }
            
            foreach my $seg_coords_aref (@contained_by_overlapping_segment_only) {
                my ($seg_lend, $seg_rend)= @$seg_coords_aref;
                $self->_add_segment_node($seg_lend, $seg_rend, \@owners);
            }
            

            foreach my $remaining_input_seg_aref (@remaining_input_segments) {
                my ($seg_lend, $seg_rend) = @$remaining_input_seg_aref;
                $self->add_segment($seg_lend, $seg_rend, $parent_feature_name); # recursive call
            }
        }
    }
    
    else {
        ## add it
        $self->_add_segment_node($lend, $rend, $parent_feature_name);
    }
    return;
}


####
sub get_all_nodes {
    my $self = shift;
   
    my $nodes_href = $self->{id_to_node};
    my @nodes = values %$nodes_href;

    @nodes = sort {$a->{lend} <=> $b->{lend}
                   ||
                       $a->{rend} <=> $b->{rend} } @nodes;;
    

    return(@nodes);
}


####
sub toString {
    my $self = shift;

    
    my $text = "";

    my @nodes = $self->get_all_nodes(); # already nicely sorted
    
    foreach my $node (@nodes) {

        my ($seg_lend, $seg_rend) = $node->get_coords();
        my @owners = sort $node->get_owners();
        
        $text .= "$seg_lend-$seg_rend\t@owners\n";
    }

    return($text);
}

####
sub find_overlapping_segment {
    my $self = shift;
    my ($lend, $rend) = @_;

    my @nodes = $self->get_all_nodes();
    foreach my $node (@nodes) {
        unless (ref $node) {
            
            confess "Error, got node thats not a ref" . Dumper($node);
        }
        
        if ($node->overlaps_coordinates($lend, $rend)) {
            return($node);
        }
    }

}


####
sub identify_all_owners {
    my $self = shift;
    

    my %all_owners;
    my @nodes = $self->get_all_nodes();
    foreach my $node (@nodes) {
        
        my @owners = $node->get_owners();
        foreach my $owner (@owners) {
            $all_owners{$owner} = 1;
        }
    }

    return(keys %all_owners);
}



###################
## Private methods
##################

####
sub _add_segment_node {
    my $self = shift;
    
    my ($lend, $rend, $parent_feature_name) = @_;

    my $segment_node = SegmentNode->new($lend, $rend, $parent_feature_name);

    my $segment_node_ID = $segment_node->get_ID();
    
    $self->{id_to_node}->{$segment_node_ID} = $segment_node;
    
    return;
}

####
sub _delete_prev_edge_to {
    my $self = shift;
    my ($node_obj, $prev_node_ID) = @_;
    unless (ref $node_obj) {
        confess "Error, need node_obj as parameter";
    }
    
    my $node_ID = $node_obj->get_ID();
    
    my $prev_edges_href = $self->{_prev_edges};
    delete($prev_edges_href->{$node_ID}->{$prev_node_ID});
    unless (%{$prev_edges_href->{$node_id}}) {
        delete($prev_edges_href->{$node_id});
    }

    return;
}

####
sub _delete_next_edge_to {
    my $self = shift;
    my ($node_obj, $next_node_ID) = @_;
    unless (ref $node_obj) {
        confess "Error, need node_obj as parameter";
    }
    
    my $node_ID = $node_obj->get_ID();
    
    my $next_edges_href = $self->{_next_edges};
    delete($next_edges_href->{$node_ID}->{$next_node_ID});
    unless (%{$next_edges_href->{$node_id}}) {
        delete($next_edges_href->{$node_id});
    }

    return;
}


#####################################################################################################
package SegmentNode;

use strict;
use warnings;
use Carp;

my $NODE_COUNTER = 0;

sub new {
    my $packagename = shift;
    my ($lend, $rend, $parent_feature_name) = @_;

    unless ($lend && $rend && $parent_feature_name) {
        die "Error, need params(lend, rend, parent_feature_name)";
    }

    my $self = { lend => $lend,
                 rend => $rend,
                 owners => {},
                 ID => ++$NODE_COUNTER,
    };

    bless ($self, $packagename);

    my @owners;
    
    if (ref $parent_feature_name eq "ARRAY") {
        @owners = @$parent_feature_name;
    }
    else {
        @owners = ($parent_feature_name);
    }

    $self->add_owners(@owners);
        
    return($self);

}


####
sub get_ID {
    my $self = shift;
    
    return($self->{ID});
}

####
sub has_same_coordinates {
    my $self = shift;
    my ($lend, $rend) = @_;

    if ($self->{lend} == $lend  && $self->{rend} == $rend) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub add_owners {
    my $self = shift;
    my (@parent_feature_names) = @_;

    foreach my $parent_feature_name (@parent_feature_names) {
        $self->{owners}->{$parent_feature_name} = 1;
    }
    
    return;
}

####
sub get_owners {
    my $self = shift;
    return(keys %{$self->{owners}});
}

####
sub is_contained_by_coordinates {
    my $self = shift;
    my ($lend, $rend) = @_;

    if ($self->{lend} >= $lend && $self->{rend} <= $rend) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub envelops_coordinates {
    my $self = shift;
    my ($lend, $rend) = @_;

    if ($self->{lend} <= $lend && $self->{rend} >= $rend) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub overlaps_coordinates {
    my $self = shift;
    my ($lend, $rend) = @_;

    if ($lend <= $self->{rend} && $rend >= $self->{lend}) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub get_coords {
    my $self = shift;
    
    return($self->{lend}, $self->{rend});

}


####
sub has_owners {
    my $self = shift;
    my @accs = @_;

    foreach my $acc (@accs) {

        if (! exists $self->{owners}->{$acc}) {
            return(0);
        }
    }

    return(1); # must have had them all
}



1; #EOM






