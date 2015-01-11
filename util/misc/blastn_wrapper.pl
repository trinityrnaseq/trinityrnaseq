#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 db query [opts]\n\n";

unless (@ARGV) {
    die $usage;
}
my $db = $ARGV[0] or die $usage;
my $query = $ARGV[1] or die $usage;

shift @ARGV;
shift @ARGV;

main: {

    my $cmd = "makeblastdb -in $db -dbtype nucl";
    &process_cmd($cmd) unless (-s "$db.nin"); # only build it once

    $cmd = "blastn -db $db -query $query -dust no @ARGV";
    &process_cmd($cmd);

    exit(0);
}

####
sub process_cmd {
    my ($cmd) = @_;

    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }
    
    return;
}

    
    
