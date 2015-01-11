package GenericNode;

use strict;
use warnings;

use Carp;

my $ID_counter = 0;

sub new {
	my $packagename = shift;
	
	my ($value) = @_;
	
	my $self = {

		_value => $value,
		_prev => {},
		_next => {},
		
		_ID => ++$ID_counter,

		_colors => {},

	};

	bless ($self, $packagename);

	return($self);
	
}

sub get_hexID {
	my $self = shift;

	my $id = sprintf("%x", $self->{_ID});
	
	return("$id");
}

sub get_ID {
    my $self = shift;
    return($self->{_ID});
}


sub get_value {
	my $self = shift;

	return($self->{_value});
}

sub set_value {
	my ($self) = shift;
	my ($value) = @_;
	
	$self->{_value} = $value;
}

sub add_prev_node {
	my $self = shift;
	my ($prev_node) = @_;

	unless (ref $prev_node) {
		croak "Error, prev_node should be a node object";
	}
	$self->{_prev}->{$prev_node} = $prev_node;
	return;
}

sub add_next_node {
	my $self = shift;
	my ($next_node) = @_;

	unless (ref $next_node) {
		croak "Error, next_node should be a node object";
	}
	
	$self->{_next}->{$next_node} = $next_node;
}

sub has_next_node {
	my $self = shift;
	my $next_node = shift;

	if (exists $self->{_next}->{$next_node}) {
		return(1);
	}
	else {
		return(0);
	}
}

sub has_prev_node {
	my $self = shift;
	my $prev_node = shift;

	if (exists $self->{_prev}->{$prev_node}) {
		return(1);
	}
	else {
		return(0);
	}
}



sub get_all_prev_nodes {
	my $self = shift;
	my @prev_nodes = grep { ref($_) } values %{$self->{_prev}};
	
	#print "Prev: " . join(",", @prev_nodes) . "\n";
	
	return(@prev_nodes);
	
}

sub get_all_next_nodes {
	my $self = shift;
	my @next_nodes = grep { ref($_) } values %{$self->{_next}};
	
	#print "Next: " . join(",", @next_nodes) . "\n";
	
	return(@next_nodes);
	
}

sub delete_prev_node {
	my $self = shift;
	my $node = shift;

	my $prev_nodes_href = $self->{_prev};
	delete $prev_nodes_href->{$node};

	return;
}

sub delete_next_node {
	my $self = shift;
	my $node = shift;

	my $next_nodes_href = $self->{_next};
	delete $next_nodes_href->{$node};
	
	return;
}


sub add_color {
	my $self = shift;

	my $color = shift;

	$self->{_colors}->{$color}++;

	return;
}

sub get_colors {
	my $self = shift;

	return(keys %{$self->{_colors}});
}




1; #EOM


