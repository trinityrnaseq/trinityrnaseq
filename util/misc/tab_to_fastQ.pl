#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
    chomp;
    my @x = split(/\t/);
    unless (scalar(@x) == 3) {
        die "Error, 3 fields not encountered for: $_";
    }

    print "\@$x[0]\n"
        . "$x[1]\n"
        . "+\n"
        . "$x[2]\n";
}

exit(0);

