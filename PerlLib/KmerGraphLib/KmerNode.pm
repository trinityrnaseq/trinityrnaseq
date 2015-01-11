package KmerNode;

use strict;
use warnings;
use Carp;
use ReadTracker;

use base qw (GenericNode);

sub new {
	my $packagename = shift;
	my ($kmer_seq, $accession) = @_;
	
	my $self = $packagename->SUPER::new($kmer_seq);
	
	$self->{_count} = 0;
	$self->{_ReadTracker} = new ReadTracker();
	
	$self->{depth} = undef;
	
	## for DP scans:
	$self->{_visited} = 0;
	$self->{_base_score} = 0;
	$self->{_sum_score} = 0;
	$self->{_best_prev} = undef; # prev_node in highest scoring path.
	$self->{_forward_scores} = []; # stores structs of { node => ref, score => $score}
	
	
	bless ($self, $packagename);

	return($self);
}


####
sub get_ReadTracker {
	my $self = shift;
	return($self->{_ReadTracker});
}


sub get_sequence {
	my $self = shift;

	return($self->get_value());
}

sub set_sequence {
	my $self = shift;
	my $sequence = shift;
	
	unless ($sequence =~ /\w/) {
		confess "Error, need sequence";
	}

	$self->set_value($sequence);

	return;
}

sub track_reads {
	my $self = shift;
	my @accs = @_;

	$self->{_ReadTracker}->track_reads(@accs);
	
	return;
}

sub get_count {
	my $self = shift;

	return($self->{_count});
}

sub set_count {
	my $self = shift;

	my ($count) = @_;

	$self->{_count} = $count;
	return;
}



####
sub toString {
        my $self = shift;
        
        my @prev_nodes = $self->get_all_prev_nodes();
        
        my @next_nodes = $self->get_all_next_nodes();
        
        my $text = "";
        
        foreach my $prev_node (@prev_nodes) {
            $text .= "P " . $prev_node->get_value() . "(" . $prev_node->get_count() . ") $prev_node " . join(",", $prev_node->get_colors()) . "\n";
        }
        $text .= "X  " . $self->get_value() . "(" . $self->get_count() . ") $self " . join(",", $self->get_colors()) . "\n";
        
        foreach my $next_node (@next_nodes) {
			$text .= "N   " . $next_node->get_value() . "(" . $self->get_count() . ") $next_node " . join(",", $next_node->get_colors()) . "\n";
		}
		
		return($text);
}


####
sub get_reads_exiting_node {
	my $self = shift;

	my @next_nodes = $self->get_all_next_nodes();
	unless (@next_nodes) {
		return();
	}

	my $read_tracker = $self->get_ReadTracker();

	my @reads_in_node = $read_tracker->get_tracked_read_indices();

	my @reads_in_next_nodes;
	foreach my $next_node (@next_nodes) {
		push (@reads_in_next_nodes, $next_node->get_ReadTracker()->get_tracked_read_indices());
	}

	my %next_node_reads = map { + $_ => 1 } @reads_in_next_nodes;

	my @exiting_reads;
	foreach my $read (@reads_in_node) {
		if ($next_node_reads{$read}) {
			push (@exiting_reads, $read);
		}
	}

	return(@exiting_reads);
}

####
sub get_reads_entering_node {
	my $self = shift;
	
	my @prev_nodes = $self->get_all_prev_nodes();
	unless (@prev_nodes) {
		return();
	}

	my $read_tracker = $self->get_ReadTracker();

	my @reads_in_node = $read_tracker->get_tracked_read_indices();

	my @reads_in_prev_nodes;
	foreach my $prev_node (@prev_nodes) {
		push (@reads_in_prev_nodes, $prev_node->get_ReadTracker()->get_tracked_read_indices());
	}

	my %prev_node_reads = map { + $_ => 1 } @reads_in_prev_nodes;

	my @entering_reads;
	foreach my $read (@reads_in_node) {
		if ($prev_node_reads{$read}) {
			push (@entering_reads, $read);
		}
	}

	return(@entering_reads);
}



1; #EOM

