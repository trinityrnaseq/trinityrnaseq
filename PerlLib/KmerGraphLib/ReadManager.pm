package ReadManager;

use strict;
use warnings;
use Carp;


my $counter = 0;
my %read_acc_to_counter;


####
sub get_read_index {
	my ($read_acc) = @_;

	if (my $read_index = $read_acc_to_counter{$read_acc}) {
		
		return($read_index);
	}

	else {

		$counter++;
		$read_acc_to_counter{$read_acc} = $counter;

		return($counter);
	}
}

1; #EOM

