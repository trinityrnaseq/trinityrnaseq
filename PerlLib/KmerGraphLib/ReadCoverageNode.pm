package ReadCoverageNode;

use strict;
use warnings;

use Carp;
use base qw (GenericNode);


## static vars
my %READ_TRACKER;
my $READ_COUNTER = 0;


sub new {
	my $packagename = shift;
	my ($node_name, $read_acc) = @_;

	my $self = $packagename->SUPER::new($node_name); # sets _value
	


	$self->{_reads} = {};
	$self->{_count} = 0;


	bless ($self, $packagename);

	$self->add_read($read_acc);


	return($self);
}

####
sub get_count {
	my $self = shift;
	return($self->{_count});
}


sub set_count {
	my $self = shift;
	my $count = shift;

	$self->{_count} = $count;
	return;
}


####
sub add_read {
	my $self = shift;
	my $read_acc = shift;

	
	my $node_name = $self->get_value();
	unless ((defined $read_acc) && $read_acc =~ /\w/) {
		confess "Error, need read accession";
	}

	my $read_tracking_no = $self->_get_or_create_tracking_number($read_acc);

	if (exists $self->{_reads}->{$read_tracking_no}) {
				
		#print "Already got read $read_acc for $node_name, count:" . $self->get_count() . "\n";
		
	}
	else {
		$self->{_reads}->{$read_tracking_no} = 1;
		$self->{_count}++;
	
		#print "-incrementing count for $read_acc for $node_name => $self->{_count}\n";
		
	}
	
	return;
}
		
####
sub get_reads {
	my $self = shift;
	
	my @reads = keys %{$self->{_reads}};

	return(@reads);
}

###
sub has_read {
	my $self = shift;
	my $read_acc = shift;
	
	if (exists $self->{_reads}->{$read_acc}) {
		return(1);
	}
	else {
		return(0);
	}
}



####
sub count_reads_in_common {
	my $self = shift;
	my $other_node = shift;

	my $count = 0;
	
	foreach my $read ($self->get_reads()) {
		if ($other_node->has_read($read)) {
			$count++;
		}
	}
	
	return($count);
}

####
sub toString {
	my $self = shift;
	
	my @prev_nodes = $self->get_all_prev_nodes();
	
	my @next_nodes = $self->get_all_next_nodes();
	
	my $text = "";
	
	foreach my $prev_node (@prev_nodes) {
	    $text .= "P " . $prev_node->get_value() . "(" . $prev_node->get_count() . ")\n";
	}
	$text .= "X  " . $self->get_value() . "(" . $self->get_count() . ")\n";
	
	foreach my $next_node (@next_nodes) {
		$text .= "N   " . $next_node->get_value() . "(" . $self->get_count() . ")\n";
	}
	
	return($text);
}

#### Private read tracking   ## this should probably be a separate class at some point, with a singleton class object.

sub _read_is_tracked {
	my $self = shift;
	my $read_acc = shift;

	if (exists $READ_TRACKER{$read_acc}) {
		return(1);
	}
	else {
		return(0);
	}
}

sub _get_read_tracking_number {
	my $self = shift;
	my $read_acc = shift;

	if (! $self->_read_is_tracked($read_acc)) {
		die "Error, read $read_acc is not tracked";
	}
	
	my $tracking_number = $READ_TRACKER{$read_acc};
	return($tracking_number);
}

sub _get_or_create_tracking_number {
	my $self = shift;
	my $read_acc = shift;

	if ($self->_read_is_tracked($read_acc)) {
		return($self->_get_read_tracking_number($read_acc));
	}
	else {
		## track this new read:
		$READ_COUNTER++;
		$READ_TRACKER{$read_acc} = $READ_COUNTER;
		return($self->_get_read_tracking_number($read_acc));
	}
}


   



1; # EOM
