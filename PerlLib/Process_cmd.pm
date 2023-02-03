package Process_cmd;

use strict;
use warnings;
use Carp;
use Cwd;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(process_cmd ensure_full_path);


sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";

	my $ret = system($cmd);
	if ($ret) {
		confess "Error, cmd:\n$cmd\n died with ret ($ret)";
	}

	return;
}


sub ensure_full_path {
    my ($path) = @_;

    my @ret_paths;

    foreach my $p (split(/,\s*/, $path)) {
        push (@ret_paths, &ensure_full_single_path($p));
    }

    my $ret_path = join(",", @ret_paths);

    return($ret_path);
    
}

sub ensure_full_single_path {
    my ($path) = @_;
    
    unless ($path =~ m|^/|) {
        $path = cwd() . "/$path";
    }

    return($path);
}



1; #EOM

