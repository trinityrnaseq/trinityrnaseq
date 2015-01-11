#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 dir/ \n\n";

my $dirname = $ARGV[0] or die $usage;


main: {

    my $cmd = "find $dirname -type f -exec rm -f {} \\;";
    &process_cmd($cmd);


    $cmd = "rm -rf $dirname";

    &process_cmd($cmd);
    
    exit(0);
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}
