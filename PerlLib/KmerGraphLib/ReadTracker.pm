package ReadTracker;

use strict;
use warnings;
use Carp;
use ReadManager;


sub new {
	my ($packagename) = shift;
	
	
	my $self = { reads => {}, # read indices, use ReadManager to hold full acc strings.

	};

	bless ($self, $packagename);

	return($self);
}

sub track_reads {
	my $self = shift;
	my @reads = @_;

	foreach my $read (@reads) {
		my $index = &ReadManager::get_read_index($read);
		$self->{reads}->{$index}++;
	}

	return;
}

sub append_to_ReadTracker {
	my $self = shift;
	my $to_add_ReadTracker = shift;

	foreach my $read_index (keys %{$to_add_ReadTracker->{reads}}) {
		$self->{reads}->{$read_index} += $to_add_ReadTracker->{reads}->{$read_index};
	}
	
	return;
}

sub get_tracked_read_indices {
	my $self = shift;

	return(keys %{$self->{reads}});
}

sub get_read_base_count {
	my $self = shift;
	my $read = shift;

	unless (defined $self->{reads}->{$read}) {
		confess "Error, no read index stored [$read] ";
	}

	return($self->{reads}->{$read});
}

	

1; #EOM

	
					 
