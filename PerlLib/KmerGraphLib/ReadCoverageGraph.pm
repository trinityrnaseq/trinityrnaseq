package main;
our $SEE;

package ReadCoverageGraph;

use strict;
use warnings;
use Carp;

use base qw(GenericGraph);
use ReadCoverageNode;

no warnings qw (recursion);

sub new {
	my $packagename = shift;
	
	my $self = $packagename->SUPER::new();

	bless ($self, $packagename);

	return($self);
}



sub get_or_create_node {
	my $self = shift;
	my ($node_name, $read_accession) = @_;
	
	if ($self->node_exists($node_name)) {
		my $node = $self->get_node($node_name);
		$node->add_read($read_accession);
		return($node);
		
	}
	else {
		# instantiate it, add it to the graph
		return($self->create_node($node_name, $read_accession));
	}
}
	

sub create_node {
	my $self = shift;

	my ($node_name, $read_accession) = @_;
	
	if ($self->node_exists($node_name)) {
		confess "Error, node $node_name already exists in the graph";
	}

	my $node = new ReadCoverageNode($node_name, $read_accession);
	
	$self->{_nodes}->{$node_name} = $node;
	
	return($node);
}


####
sub get_nodes_sorted_by_count_desc {
	my $self = shift;
	my @nodes = $self->get_all_nodes();

	@nodes = reverse sort {$a->{_count}<=>$b->{_count}} @nodes;
	
	return(@nodes);
}


sub find_maximal_path_including_node {
	my $self = shift;
	
	my ($nucleating_node, $max_recurse_depth) = @_;

	my $path_forward_aref = [$nucleating_node];
	my $depth_forward = 0;
	my $sum_forward_count = 0;
	do {
		my $start_node = $path_forward_aref->[-1];
		($path_forward_aref, $sum_forward_count, $depth_forward) = $self->extend_path("next", $start_node, $path_forward_aref, $max_recurse_depth, 0, 0);
		
		if ($main::SEE) {
			print "Forward.\n";
			&print_path(@$path_forward_aref);
		}

	} while ($depth_forward > 0);

	if ($main::SEE) {
		print "Forward, done.\n";
		&print_path(@$path_forward_aref);
	}
	
	
	my $path_reverse_aref = [$nucleating_node];
	my $depth_reverse = 0;
	my $sum_reverse_count = 0;
	do {
		my $start_node = $path_reverse_aref->[-1];
		($path_reverse_aref, $sum_reverse_count, $depth_reverse) = $self->extend_path("prev", $start_node, $path_reverse_aref, $max_recurse_depth, 0, 0);
		
		if ($main::SEE) {
			print "Reverse:\n";
			&print_path(@$path_reverse_aref);
		}

	} while ($depth_reverse > 0);
	

	if ($main::SEE) {
		print "Reverse, done.\n";
		&print_path(@$path_reverse_aref);
	}
	
	## unwrap path

	# pull out the nucleating node, should be first one in each path.
	shift @$path_forward_aref;
	shift @$path_reverse_aref;
	
	my @path = ( (reverse @$path_reverse_aref), $nucleating_node, @$path_forward_aref);
	
	if ($main::SEE) {
		print "Done.\n";
		&print_path(@path);
	}
	
	return(@path);
	
}


####
sub extend_path {
	my $self = shift;
	
	my ($direction, 
		$node, 
		$curr_path_list_aref,
		$max_recurse_depth, 
		$sum_counts, 
		$curr_depth) = @_;

	
	my $path_length = scalar (@$curr_path_list_aref);
	
	print "extending $direction from " . $node->get_value() . ", K:$path_length S:$sum_counts, D:$curr_depth\n" if $main::SEE;
	
	#print join("\t", @_) . "\n";
	
	## curr_path_list_aref should include the incoming node already

	if ($curr_depth >= $max_recurse_depth) {
		return($curr_path_list_aref, $sum_counts, $curr_depth);
	}

	## explore connected nodes
	my @other_nodes;
	if ($direction eq 'next') {
		@other_nodes = $node->get_all_next_nodes();
	}
	else {
		@other_nodes = $node->get_all_prev_nodes();
	}

	## only explore those other_node's that are not already seen along the current path
	my @unseen_nodes;

	foreach my $other_node (@other_nodes) {
		unless (grep {$_ == $other_node} @$curr_path_list_aref) {
			push (@unseen_nodes, $other_node);
		}
	}
	
	#print "Curr path nodes: " . join(" ", @$curr_path_list_aref) . "\n";

	@other_nodes = @unseen_nodes; # reset to those that haven't been seen already.

	#print "\tOther nodes: " . join(" ", @other_nodes) . "\n\n";

	if (@other_nodes) {
		## examine possible paths:
		
		if ($main::SEE) {
			print "Extending from:\n" . $node->get_value() . " (" . $node->get_count() . ") Depth:$curr_depth to\n";
			foreach my $other_node (@other_nodes) {
				print $other_node->get_value() . " (" . $other_node->get_count() . ")\n";
			}
			print "\n";
		}
		
		
		my @alt_paths;
		foreach my $other_node (@other_nodes) {
			#print $other_node->toString() . "\n";
			
			my ($path_list_aref, $counts, $depth) = $self->extend_path($direction,
																	   $other_node,
																	   [@$curr_path_list_aref, $other_node], # tack it on to the list
																	   $max_recurse_depth, 
																	   $sum_counts + $other_node->get_count(), 
																	   $curr_depth+1);
			
			push (@alt_paths, [$path_list_aref, $counts, $depth]);
		}
		
		## select the greatest one
		@alt_paths = sort {  #$a->[2] <=> $b->[2]   ## Perhaps include extension length
							 #	 ||
								 $a->[1] <=>$b->[1] } @alt_paths;
		
		my $top_path = pop @alt_paths;
		
		my ($path_list_aref, $counts, $depth) = @$top_path;
		return($path_list_aref, $counts, $depth);
	}
	else {
		## no extensions possible
		return($curr_path_list_aref, $sum_counts, $curr_depth);
	}
}



####
sub print_path {
	my @nodes = @_;

	my $counter = 0;
	foreach my $node (@nodes) {
		$counter++;
		printf("%4s", $counter);
		print " " . $node->get_value() . " C:" . $node->get_count() . "\n";
	}
	print "\n";
	
	return;
}



1; #EOM
