#!/usr/bin/env perl

use strict;
use warnings;

my @files = @ARGV;

foreach my $file (@files) {
    open (my $fh, $file) or die $!;
    while (<$fh>) {
        chomp;
        print "$_\n";
    }
    close $fh;
}

exit(0);

