package main;
our $SEE;

package KmerGraph;

use strict;
use warnings;
use Carp;
use KmerNode;
use ReadTracker;

use base qw (ReadCoverageGraph);

no warnings qw (recursion);


my $MIN_DISPLAY_SEQ_TAG_LENGTH = 10; # lower than this, and shows up in the dot file as a label.

sub new {
	my $packagename = shift;
	my $KmerLength = shift;
	
	unless ($KmerLength > 1) { croak "Error, need KmerLength > 1"; }
	
	my $self = { 
		_nodes => {},  #kmerSeq => nodeAddress until after compaction... becomes nodeaddress => $nodeaddress
		KmerLength => $KmerLength,
		compacted => 0,
        color_to_accs => {},
        
	};
	
	bless ($self, $packagename);

	return($self);
}


sub add_sequence_to_graph {
	my $self = shift;
	my ($acc, $sequence, $weight, $color) = @_;  ## Sequence type can be undef, ReferenceTranscript, or GreedyPath
	
	unless (defined($weight) && $weight =~ /^\d+/) {
		confess "Error, params: (seq, weight, color), and weight must be numeric.";
	}
	
	if ($self->{compacted}) {
		confess "Cannot add more sequences to the graph after it has been compacted";
	}
	
    push (@{$self->{color_to_accs}->{$color}}, $acc); # add this accession to color indicator. # prefer 1-1 relationship here.
    
	$sequence = uc $sequence;
	
	my $KmerLength = $self->{KmerLength};

	if (length($sequence) < $KmerLength+1) {
		warn "sequence $sequence is less than KmerLength $KmerLength+1; skipping it\n";
		return;
	}

	my $prevKmer = substr($sequence, 0, $KmerLength);
	my $prevKmerNode = $self->get_or_create_node($prevKmer);
	$prevKmerNode->{_count}++;
	$prevKmerNode->track_reads($acc);
   
	if ($color) {
		$prevKmerNode->add_color($color);
	}
	
	for (my $i = 1; $i <= length($sequence)-$KmerLength; $i++) {
		my $nextKmer = substr($sequence, $i, $KmerLength);
		my $nextKmerNode = $self->get_or_create_node($nextKmer);
		$nextKmerNode->track_reads($acc);
		
		if ($color) {
			$nextKmerNode->add_color($color);
		}
				
		$nextKmerNode->{_count}++;

		my $prev_edge_weight = $self->get_edge_count($prevKmerNode, $nextKmerNode);

		## only link together kmers that contain only recognizable nucleotides
		$self->link_adjacent_nodes($prevKmerNode, $nextKmerNode, $weight);
		
		my $edge_weight = $self->get_edge_count($prevKmerNode, $nextKmerNode);
		if ($edge_weight < $weight) {
			confess "Error, added weight $weight to edge $prevKmerNode -> $nextKmerNode and didn't stick";
		}

		if ($prev_edge_weight > $edge_weight) {
			confess "Error, after adding new edge, prev weight of $prev_edge_weight became: $edge_weight ";
		}
		#print "$prevKmerNode -> $nextKmerNode : $prev_edge_weight -> $edge_weight\n";
		
		$prevKmerNode = $nextKmerNode;
		$prevKmer = $nextKmer;
	}
	
	return;
}

sub get_or_create_node {
        my $self = shift;
        my ($node_name, $read_accession) = @_;
        
        if ($self->node_exists($node_name)) {
                my $node = $self->get_node($node_name);
				return($node);
                
        }
        else {
                # instantiate it, add it to the graph
                return($self->create_node($node_name));
        }
}
        

sub create_node {
        my $self = shift;

        my ($node_name) = @_;
        
        if ($self->node_exists($node_name)) {
                confess "Error, node $node_name already exists in the graph";
        }

        my $node = new KmerNode($node_name);
        
        $self->{_nodes}->{$node_name} = $node;
        
        return($node);
}




####
sub compact_graph {
	my $self = shift;

	$self->validate_graph();

	my $kmer_length = $self->{KmerLength};

	unless ($self->{compacted}) {
		$self->_prep_graph_for_compaction();
		
		$self->{compacted} = 1;
	}

	my $COMPACTED = 0; # flag

	
	## identify all nodes that are branched on the left and not branched on the right
	## then join together the unbranched nodes starting from these and walking right.

	my @nodes = $self->get_all_nodes();

	## merge unbranched nodes in runs.




	####
	####  \
	####    (-)---
	#### /

	# or

	#  (-)---


	# or

	####      (-)------------
	####    /
	####  -
	####    \ 
	####     (-) ------------


	# or, single line with incompatible decorations:
	
	#    aaaaaaaaaa b bbbbbbbbbb
	#    ----------(-)----------


	my @init_nodes;
	foreach my $node (@nodes) {
		my @prev_nodes = $node->get_all_prev_nodes();

		my @next_nodes = $node->get_all_next_nodes();
		
		if 
			
			(  (scalar(@next_nodes) == 1)
			  &&
			
			   (
				



				( 
				  (scalar(@prev_nodes) != 1)  # branched before or no node before. 
				  ||
				  (scalar (@prev_nodes) == 1 && scalar($prev_nodes[0]->get_all_next_nodes() > 1) ) # previous node branches off.
				  
				)

				||
				
				## decoration switch
				(scalar(@prev_nodes) == 1 && (! &_compatible_decorations($node, $prev_nodes[0])) && (&_compatible_decorations($node, $next_nodes[0])) )
				
			   )

			)
		{
			# got one
			push (@init_nodes, $node);
		}
		
	}
	
	print "Got " . scalar(@init_nodes) . " init nodes.\n" if $main::SEE;
	
	#foreach my $node (@init_nodes) {
	#	print "Init node $node:\n" . $node->toString() . "\n";# if $main::SEE;
	#}
	

	## Collapse neighboring nodes before branching
	foreach my $node (@init_nodes) {
		
		#print "Init node $node:\n" . $node->toString() . "\n";# if $main::SEE;
				
		my @unbranched_path;
		
		my ($next_node) = $node->get_all_next_nodes();
		
		unless ($next_node) {
			print "warning, init node now lacks next nodes...  skipping.\n"; # potential bug?  or circular dealt with already?
			next;
		}

		my %seen = ($node => 1);
		
	    while (1) {
			# avoid circularity
			if ($seen{$next_node}) {
				print "Seen $next_node already... avoiding circularity.\n" if $main::SEE;
				last;
			}
			$seen{$next_node} = 1;
			
			unless (&_compatible_decorations($node, $next_node)) {
				print "Incompatible decorations.\n" if $main::SEE;
				last;
			}
			
			# check for left branching on next node
			my $num_prev_nodes = scalar ($next_node->get_all_prev_nodes());
			if (scalar $num_prev_nodes > 1) {
				print "Got $num_prev_nodes prev nodes for $next_node; terminating extension.\n" if $main::SEE;
				last;
			}
			
			push (@unbranched_path, $next_node);
			my @next_nodes = $next_node->get_all_next_nodes();
		
			if (scalar @next_nodes != 1) {
				last;
			}
			$next_node = shift @next_nodes;
		
		}
		
		if ($main::SEE) {
			print "Single unbranched path:\n";
			$self->print_path($node, @unbranched_path);
		}
		
		foreach my $next_node (@unbranched_path) {
			# merge pairs
			$self->_append_node($node, $next_node);
			$COMPACTED = 1;
						
		}
		
		if ($main::SEE) {
			print "Collapsed:\n";
			$self->print_path($node);
			print "\n\n\n";
		}

	}

	
	$self->validate_graph();

	
	return ($COMPACTED);
}

####
sub _append_node {
	my $self = shift;
	my ($node, $next_node) = @_;

	## In compacted mode already.

	if (! _compatible_decorations($node, $next_node)) {
		confess "Error, trying to join nodes with incompatible decorations.";
	}
	
	my $KmerLength = $self->{KmerLength};

	
	# add counts
	$node->set_count($node->get_count() + $next_node->get_count());
	
	my $orig_node_seq = $node->get_sequence();
	my $orig_next_node_seq = $next_node->get_sequence();

	my $new_node_sequence = $orig_node_seq . substr($orig_next_node_seq, $KmerLength-1); # exclude K-1 prefix.
	
	# add terminal base
	$node->set_sequence($new_node_sequence);

	# add read evidence
	$node->{_ReadTracker}->append_to_ReadTracker($next_node->{_ReadTracker});
	
	# sever connection between node and next
	$node->delete_next_node($next_node);
	$next_node->delete_prev_node($node);
	
	# sever connection between next and next-next nodes, and reconnect next-next's with node.
	
	my @next_node_next_nodes = $next_node->get_all_next_nodes();
	
	my %next_next_node_edge_count;

	## update edge counts.
	foreach my $next_next_node (@next_node_next_nodes) {
		my $edge_count = $self->get_edge_count($next_node, $next_next_node);
		
		$next_next_node_edge_count{$next_next_node} = $edge_count;
	}
	
	$self->delete_node_from_graph($next_node);

	foreach my $next_next_node (@next_node_next_nodes) {
		$self->link_adjacent_nodes($node, $next_next_node, $next_next_node_edge_count{$next_next_node});
		
	}


	
	return;
}


####
sub toGraphViz {
	my $self = shift;

	my (%settings) = @_;

	my @nodes = $self->get_all_nodes();

	my $text = "digraph G {\n";

	$text .= 
		#"node [width=0.1,height=0.1,fontsize=10,shape=point];\n"
		"node [width=0.1,height=0.1,fontsize=10];\n"
		. "edge [fontsize=12];\n"
		. "margin=1.0;\n"
		. "rankdir=LR;\n"
		. "labeljust=l;\n";
	
	
	foreach my $node (sort {$a->get_ID() cmp $b->get_ID()}  @nodes) {
		
		my @prev_nodes = $node->get_all_prev_nodes();
		my @next_nodes = $node->get_all_next_nodes();
		
		my $sequence = $node->get_sequence();

		if (my $min_len = $settings{no_short_singletons}) {
			
			if ( scalar(@prev_nodes) == 0
				 &&
				 scalar(@next_nodes) == 0
				 &&
				 length($sequence) < $min_len) {

				next;
			}
		}
		
		

		my $seqLen = length($sequence);
		
		my $count = $node->get_count();
		
		my $len = $seqLen;
		
		if (scalar(@prev_nodes)==0) { # no left connecting node:
			$len = $seqLen - ($self->{KmerLength} - 1);
		}
		
		my $avg_count = $count; #int($count/$len + 0.5);
		
		my $id = hex($node->get_ID());
		
		my $atts = "";
		
		if (length($sequence) < $MIN_DISPLAY_SEQ_TAG_LENGTH) {
			$atts = "-$sequence";
		}
		else {
			$atts = "-" . substr($sequence, 0, 3) . "..." . substr($sequence, -3);
		}
		
		my $depth ="";
		if (defined $node->{depth}) {
			$depth = $node->{depth};
		}
		
		#$text .= "\t$id \[label=\"$id-L$seqLen-C$avg_count$atts\"$color];\n";
		$text .= "\t$id \[label=\"D:$depth-L$seqLen-C$avg_count$atts\"];\n";
		
		
		foreach my $next_node (@next_nodes) {
			
			my $next_id = hex($next_node->get_ID());
						
			my @colors_in_common = &get_colors_in_common($node, $next_node);
			
			my $edge_weight = $self->get_edge_count($node, $next_node);
			
			if (@colors_in_common) {
				foreach my $color (@colors_in_common) {
					
                    my $label = "$edge_weight";
                    unless ($color eq 'black') { # reserved for general data
                        my @accs = @{$self->{color_to_accs}->{$color}};
                        $label = join(",", @accs);
                    }
					#$text .= "\t$id->$next_id [label=$edge_weight, color=\"$color\"];\n";
					$text .= "\t$id->$next_id [label=\"$label\", color=\"$color\"];\n";
                    
				}
			}

		}
	}
		
	$text .= "}\n";
	
	return($text);
}


sub get_colors_in_common {
	my ($node, $next_node) = @_;

	my @colorsA = $node->get_colors();
	my @colorsB = $next_node->get_colors();
	
	my %counts;
	
	foreach my $color (@colorsA, @colorsB) {
		$counts{$color}++;
	}
	my @colors = grep { $counts{$_} > 1 } keys %counts;

	return(@colors);
}



####
sub _prep_graph_for_compaction {
	my $self = shift;

	my @nodes = $self->get_all_nodes();

	my %new_node_set;  ## change lookup so it's not based on sequence anymore.

	foreach my $node (@nodes) {
		$new_node_set{"$node"} = $node;
	}

	$self->{_nodes} = \%new_node_set;
}


####
sub delete_node_from_graph {
	my $self = shift;
	my ($node) = @_;

	$self->prune_nodes_from_graph($node);
	
	if ($self->{compacted}) {
		
		delete $self->{_nodes}->{$node};
		
	}
	else {
		my $kmer = $node->get_sequence();
		delete $self->{_nodes}->{"$kmer"};
	}
	
	return;
}


####
sub purge_nodes_below_count {
	my $self = shift;

	my ($min_count) = @_;

	my @nodes = $self->get_all_nodes();

    my $num_deleted = 0;

	foreach my $node (@nodes) {

        if ($node->get_count() < $min_count) {
			
			my @colors = $node->get_colors();
            unless (scalar(@colors) == 1 && $colors[0] eq "black") { 
                next;
                ## skipping base data (reads) below coverage limit
            }
            
            $self->delete_node_from_graph($node);
            
            $num_deleted++;
        }
	}
	
    print STDERR "  -purged nodes below count($min_count), deleted $num_deleted nodes.\n";

	return;
}


####
sub validate_graph {
	my $self = shift;

	print STDERR "Validating graph.\n";
	
	my @nodes = $self->get_all_nodes();

	my %node_ids_in_graph;

	foreach my $node (@nodes) {
	   
		my $node_id = $node->get_ID();
		$node_ids_in_graph{$node_id} = $node;
	}

	foreach my $node (@nodes) {

		my $id = $node->get_ID();
		
		my @prev_nodes = $node->get_all_prev_nodes();
		foreach my $prev_node (@prev_nodes) {
			my $prev_id = $prev_node->get_ID();
			
			unless (exists ($node_ids_in_graph{$prev_id})) {
				print STDERR "Error, $prev_id exists as prev_node to $id, but not in graph\n";
				print STDERR $prev_node->get_sequence() . "\n" . $prev_node->toString();
				confess;
			}
			
			unless ($prev_node->has_next_node($node)) {
				confess "Error, prev_node: $prev_id lacks current next node $id  as link.\n";
			}
			
			
		}

		my @next_nodes = $node->get_all_next_nodes();
		foreach my $next_node (@next_nodes) {
			my $next_id = $next_node->get_ID();
			
			unless (exists ($node_ids_in_graph{$next_id})) {
				print STDERR "Error, $next_id exists as next_node to $id, but not in graph\n";
				print STDERR $next_node->get_sequence() . "\n" . $next_node->toString();
				confess;
			}
		
			unless ($next_node->has_prev_node($node)) {
			    confess "Error, next node: $next_id lacks current node $id as prev link.\n";
			}
			
			
		}
		
	}
	
	return;
}

####
sub print_path {
	my $self = shift;
	my @nodes = @_;

	my $counter = 0;
	my $prev_node;
	
	foreach my $node (@nodes) {
		$counter++;
		printf("%4s", $counter);
		
		my $edge_count = "";
		if ($prev_node) {
			$edge_count = " edge_weight: " . $self->get_edge_count($prev_node, $node);
		}
		
		print " " . $node->get_value() . " C:" . $node->get_count() . " " . join(",", $node->get_colors()) . "$edge_count\n";
		
		$prev_node = $node;

	}
	print "\n";
	
	return;
}


sub _compatible_decorations {
	my ($nodeA, $nodeB) = @_;


	my @colorsA = $nodeA->get_colors();
	my @colorsB = $nodeB->get_colors();

	unless (@colorsA || @colorsB) {
		# no decorations.
		return(1);
	}

	my %counts;
	foreach my $color (@colorsA, @colorsB) {
		$counts{$color}++;
	}
	
	my @colors_unique = grep { $counts{$_} == 1 } keys %counts; 
		
	if (@colors_unique) {
		#print STDERR "unique colors: @colors_unique\n";
		return(0);
	}

	else {
		return(1);
	}
}




####
sub score_nodes {
	my $self = shift;
	
	$self->assign_depth_to_nodes();

	##  score = sum prev counts for those reads that are consistent, excluding those that are from bifurcating reads.

	
	my @nodes = sort {$a->{depth}<=>$b->{depth}} $self->get_all_nodes();

	## assign base scores:
	foreach my $node (@nodes) {
		
		## define base scores: contributions by reads that are not in prev or next nodes.
		my %neighboring_reads = $self->gather_read_indices($node->get_all_prev_nodes(), $node->get_all_next_nodes());
		
		my $base_score = 0;

		my $read_tracker = $node->get_ReadTracker();
		my @reads = $read_tracker->get_tracked_read_indices();
		foreach my $read (@reads) {
			if ($node->{depth} == 0 || ! $neighboring_reads{$read}) {
				my $count = $read_tracker->get_read_base_count($read);
				
				$base_score += $count;
				#print "$node has read [$read] with count [$count]\n";

			}
		}

		$node->{_base_score} = $base_score;
		$node->{_sum_score} = $base_score; # init
		#print "$node : base score = $base_score\n";
	}
	
	
	## assign sum scores
	foreach my $node (@nodes) {
		my %indices_exclude;
		
		my $read_tracker = $node->get_ReadTracker();
		
		my $highest_prev_score = $node->{_sum_score};
		my $highest_prev_node = undef;

		foreach my $prev_node ($node->get_all_prev_nodes()) {
			## want reads that are in current node and prev node, but not other next nodes of prev node (exclusive to current node).
			my $sum_score = $prev_node->{_sum_score};

			my @other_prev_next_nodes = grep { $_ ne $node } $prev_node->get_all_next_nodes();
			
			my %reads_ignore = $self->gather_read_indices(@other_prev_next_nodes);
			
			foreach my $read ($read_tracker->get_tracked_read_indices()) {
				
				if (! $reads_ignore{$read}) {
					
					$sum_score += $read_tracker->get_read_base_count($read);
				}
			}
			
			push (@{$prev_node->{_forward_scores}}, { next_node => $node, 
													  score => $sum_score, 
				  } );
			

			if ($sum_score > $highest_prev_score) {
				$highest_prev_score = $sum_score;
				$highest_prev_node = $prev_node;
			}
		}
		
		$node->{_sum_score} = $highest_prev_score;
		$node->{_best_prev} = $highest_prev_node;
	}
			

	## find highest sum score
	my $highest_sum_score = 0;
	my $best_node = undef;
	
	foreach my $node (@nodes) {
		#print "Node: $node, has sum score: " . $node->{_sum_score} . "\n";
		if ($node->{_sum_score} > $highest_sum_score) {
			$highest_sum_score = $node->{_sum_score};
			$best_node = $node;
		}
	}

	
	my @path_nodes;
	my $prev_node = $best_node;
	while ($prev_node) {
		#print "Backtracked: $prev_node " . $prev_node->{_sum_score} . "\n";
		#print $prev_node->toString();
		push (@path_nodes, $prev_node);
		$prev_node = $prev_node->{_best_prev};
	
	}

	@path_nodes = reverse @path_nodes;
	
	return (@path_nodes);
	

}


####
sub extract_sequence_from_path {
	my $self = shift;
	my (@ordered_nodes) = @_; # ordered left to right
	
	my @seqs;
	
	foreach my $node (@ordered_nodes) {
		
		$node->{_visited} = 1;

		push (@seqs, $node->get_sequence());
		
	}

	
	my $kmer_length = $self->{KmerLength};
	
	my $final_seq = shift @seqs; ## first one gets full kmer treatment.
	
	foreach my $seq (@seqs) {
		$seq = substr($seq, $kmer_length-1);
			
		$final_seq .= $seq;
	}
   

	return($final_seq);
}



####
sub gather_read_indices {
	my $self = shift;
	my (@nodes) = @_;

	my %indices;

	foreach my $node (@nodes) {
		
		my $read_tracker = $node->get_ReadTracker();
		
		my @read_indices = $read_tracker->get_tracked_read_indices();
		
		foreach my $read (@read_indices) {
			$indices{$read}=1;
		}
	}

	return(%indices);
}


####
sub assign_depth_to_nodes {
	my $self = shift;
	print STDERR "Assigning depth to nodes\n";	
	
	my @root_nodes = $self->get_root_nodes();

	## base cases.
	foreach my $root (@root_nodes) {
		$root->{depth} = 0; 
	}
	
	foreach my $node ($self->get_all_nodes()) {
		if (! defined ($node->{depth})) {
			$self->assign_depth($node, {});
		}
	}
	
	
	return;
}


####
sub assign_depth {
	my $self = shift;
	my $node = shift;
	my $seen_href = shift;

	my @ancestors = $node->get_all_prev_nodes();
	
	my $max_depth = 0;
	
	foreach my $ancestor (@ancestors) {
				
		my $depth;
		if (defined($ancestor->{depth})) {
			$depth = 1 + $ancestor->{depth};
		}
		else {
			
			if ($seen_href->{$ancestor}) { 
				# print "Breaking cycle: $ancestor\n" . join("\n", keys %$seen_href) . "\n";
				next; 
			}
			$seen_href->{$ancestor} = 1;
			
			$depth = 1 + $self->assign_depth($ancestor, $seen_href);
		}
		if ($depth > $max_depth) {
			$max_depth = $depth;
		}
	}
	
	# print "$node\tdepth: $max_depth\n";
	$node->{depth} = $max_depth;
	return($max_depth);
}


####
sub all_nodes_visited {
	my $self = shift;
	
	my @nodes = $self->get_all_nodes();

	foreach my $node (@nodes) {

		if (! $node->{_visited}) {
			return(0);
		}
	}

	return(1); # all visited.
}

####
sub extract_path_from_unvisited_node {
	my $self = shift;
	
	my @nodes = $self->get_all_nodes();
	
	my $highest_scoring_unvisited_node = undef;
	my $highest_score = -1;
	
	foreach my $node (@nodes) {
		if (! $node->{_visited}) {
			if ($node->{_base_score} > $highest_score) {
				$highest_score = $node->{_base_score};
				$highest_scoring_unvisited_node = $node;
			}
		}
	}

	unless ($highest_scoring_unvisited_node) {
		confess "Error, no unvisited node found";
	}

	my @best_path;

	my %seen; # don't follow cycles.  FIXME: cycles shouldn't be part of the scoring, though... bug somewhere upstream?

	my $prev_node = $highest_scoring_unvisited_node;
	while ($prev_node && ! $prev_node->{_visited} && ! $seen{$prev_node}) {
		print "walking prev_node: $prev_node\n" if $main::SEE;
		$seen{$prev_node} = 1;
		unshift(@best_path, $prev_node);
		$prev_node = $prev_node->{_best_prev};
	}

	## now track forward.
	my $next_node = $highest_scoring_unvisited_node;
	delete $seen{$next_node}; # already used above.
	while ($next_node && ! $next_node->{_visited} && ! $seen{$next_node}) {
		print "walking next_node: $next_node\n" if $main::SEE;
		$seen{$next_node} = 1;
		my @structs = @{$next_node->{_forward_scores}};
		if (@structs) {
			@structs = sort {$a->{score}<=>$b->{score}} @structs;
			my $best_struct = pop @structs;
			push (@best_path, $next_node);
			$next_node = $best_struct->{next_node};
		}
		else {
			$next_node = undef;
		}
	}

	return(@best_path);
}



####
sub prune_low_weight_edges {
	my $self = shift;
	my %params = @_;

	my $edge_weight_threshold = $params{edge_weight_threshold};


	my @nodes = $self->get_all_nodes();
	
	my @edges_to_prune;

	foreach my $node (@nodes) {
		my @prev_nodes = $node->get_all_prev_nodes();
		
		foreach my $prev_node (@prev_nodes) {
			
			my $edge_ratio = $self->compute_edge_support_ratio($prev_node, $node);
			
			#print "EdgeRatio: $edge_ratio\n";
			
			if ($edge_ratio < $edge_weight_threshold) {
				push (@edges_to_prune, [$prev_node, $node]);
				
			}
		}
	}
	

	if (@edges_to_prune) {
		
		foreach my $edges_to_prune (@edges_to_prune) {
			
			my ($prev_node, $node) = @$edges_to_prune;
			$self->prune_edge($prev_node, $node);
		}
		
		return(scalar(@edges_to_prune));
	}
	else {
		return(0);
	}

}

####
sub prune_singletons {
	my $self = shift;

	my ($min_seq_length) = @_;

	my @nodes = $self->get_all_nodes();

	foreach my $node (@nodes) {

		if (  (! $node->get_all_prev_nodes())
			  &&
			  (! $node->get_all_next_nodes()) ) {
			
			if (length ($node->get_sequence()) < $min_seq_length) {
			
				$self->delete_node_from_graph($node);
			}
		}
	}
			
	return;
}


####
sub compute_edge_support_ratio {
	my $self = shift;
	my ($prev_node, $node) = @_;

	## ratio = reads in common / (all reads out of prev UNION all reads into node)
	
	my @reads_exiting_prev_node = $prev_node->get_reads_exiting_node();
	my %exiting_reads = map { + $_ => 1 } @reads_exiting_prev_node;

	my @reads_entering_node = $node->get_reads_entering_node();
	my %entering_reads = map { + $_ => 1 } @reads_entering_node;

	my %all_reads = map { + $_ => 1 } (@reads_exiting_prev_node, @reads_entering_node);

	my @all = keys %all_reads;
	my @shared;
	
	foreach my $read (@all) {
		if ($entering_reads{$read} && $exiting_reads{$read}) {
			push (@shared, $read);
		}
	}

	my $ratio = scalar(@shared) / scalar(@all);

	return($ratio);
}

	

####
sub prune_dangling_nodes {
	my $self = shift;
	my %params = @_;

	my $min_leaf_node_length = $params{min_leaf_node_length};
	my $min_leaf_node_avg_cov = $params{min_leaf_node_avg_cov};

	unless (defined ($min_leaf_node_length)
					 && 
			defined ($min_leaf_node_avg_cov) ) {

		confess "Error, must define values for both: min_leaf_node_length && min_leaf_node_avg_cov";
	}

	my @nodes = $self->get_all_nodes();
	
	my @nodes_to_prune;

	my $kmer_length = $self->{KmerLength};

	foreach my $node (@nodes) {

		## check to see if it's a dangling node.
		if (  (! $node->get_all_prev_nodes())
			  ||
			  (! $node->get_all_next_nodes())
			) {

			my $length = length($node->get_sequence());
			
			my $cov = $node->get_count();

			my $avg_cov = $cov / ($length - ($kmer_length - 1));

			if ($avg_cov < $min_leaf_node_avg_cov && $length < $min_leaf_node_length) {
				push (@nodes_to_prune, $node);
			}

		}
	}

	if (@nodes_to_prune) {

		foreach my $node (@nodes_to_prune) {
			
			$self->prune_nodes_from_graph($node);
		}

		return(scalar(@nodes_to_prune)); # nodes pruned.
	}

	else {
		return(0); # no nodes pruned.
	}

}



1; #EOM

