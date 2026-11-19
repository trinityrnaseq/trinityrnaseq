package Process_cmd;

use strict;
use warnings;
use Carp;
use Cwd;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(process_cmd ensure_full_path);
our @EXPORT_OK = qw(describe_exit_status);


sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";

	my $ret = system($cmd);
	if ($ret) {
		confess "Error, cmd:\n$cmd\n died with ret ($ret)" . &describe_exit_status($ret);
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


####
# Translate a raw system() return value into a human-readable explanation.
sub describe_exit_status {
    my ($ret) = @_;

    return "" unless $ret;
    
    if ($ret == -1) {
        return " [failed to execute command: $!]";
    }

    my $signal = $ret & 127;
    my $exit_code = $ret >> 8;

    if ($signal) {
        my %sig_names = (2 => 'SIGINT', 6 => 'SIGABRT', 9 => 'SIGKILL', 11 => 'SIGSEGV', 15 => 'SIGTERM');
        my $name = $sig_names{$signal} || "signal $signal";
        my $msg = " [process was killed by $name";
        if ($signal == 9) {
            $msg .= ": this usually means it ran out of memory (killed by the OS or job scheduler). Check available RAM vs. --max_memory / --CPU settings";
        }
        elsif ($signal == 11) {
            $msg .= ": segmentation fault";
        }
        return "$msg]";
    }

    # shell-reported exit codes
    if ($exit_code == 127) {
        return " [exit code 127: command not found - check that all required tools are installed and in your PATH]";
    }
    elsif ($exit_code == 126) {
        return " [exit code 126: command found but not executable]";
    }
    elsif ($exit_code == 137) {
        return " [exit code 137: process was killed (SIGKILL), usually because it ran out of memory. Check available RAM vs. --max_memory / --CPU settings]";
    }
    elsif ($exit_code == 134) {
        return " [exit code 134: process aborted (SIGABRT), often due to running out of memory or a crash]";
    }
    elsif ($exit_code == 139) {
        return " [exit code 139: segmentation fault]";
    }
    
    return " [exit code $exit_code]";
}


1; #EOM

