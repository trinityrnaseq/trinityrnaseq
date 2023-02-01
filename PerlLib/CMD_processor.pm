package CMD_processor;

use strict;
use warnings;

use Carp;

use vars qw (@ISA @EXPORT); ## set in script using this module for verbose output.
@ISA = qw(Exporter);
@EXPORT = qw(process_cmd process_parallel_cmds);

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret ($ret)";
    }

    return; # all good.
}

####
sub process_parallel_cmds {
    my (@cmds) = @_;

    print "Processing commands:\n" . join("\n", @cmds) . "\n\n";
    

    foreach my $cmd (@cmds) {
        my $pid = fork();
        unless ($pid) {
            # child
            &process_cmd($cmd);
            exit(0);
        }
    }

    my $process_failed = 0;
    while (my $pid = wait()) {
        if ($pid < 0) {
            last;
        }
        else {
            if ($?) {
                print STDERR "-process $pid exited ($?) indicating error.\n";
                $process_failed = 1;
            }
        }
    }

    if ($process_failed) {
        die "Error, at least one command failed.  See above errors for more details. ";
    }

    return; # all good.
}

1; #EOM

