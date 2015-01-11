package Process_cmd;

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(process_cmd);


sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";

	my $ret = system($cmd);
	if ($ret) {
		confess "Error, cmd:\n$cmd\n died with ret ($ret)";
	}

	return;
}


1; #EOM

