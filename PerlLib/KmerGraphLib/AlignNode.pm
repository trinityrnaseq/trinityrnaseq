package AlignNode;

use base qw (ReadCoverageNode);

## Instead of storing Kmers, will store base number and strand positions.

sub new {
	my $packagename = shift;
	my ($stranded_base, $read_accession) = @_;
	
	my $self = $packagename->SUPER::new($stranded_base, $read_accession);
	
	bless ($self, $packagename);

	return($self);
}


1; #EOM
